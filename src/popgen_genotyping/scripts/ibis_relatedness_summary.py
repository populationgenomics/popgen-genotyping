"""Pure functions and job-side steps for the IBIS relatedness table.

jobs/ibis_relatedness_job.py queues these steps as Hail Batch jobs, and the unit tests
(test/test_ibis_relatedness_summary.py) import the functions directly. This module needs
only pandas, so it runs in any image that has it, as a single file staged into the job: it must
not import from the rest of this package.

Two steps run inside the batch, each as `python3 ibis_relatedness_summary.py <step> ...`:

prepare-maps
    Takes the .bim files plink2 wrote with centimorgan positions (autosomes and X), makes cM
    non-negative, checks the maps are usable and writes the corrected .bim files plus the
    autosomal map length. plink2 --cm-map extrapolates variants before a chromosome's first map
    row, which gives them a slightly negative cM. Setting those to 0 would make the first two
    markers of chromosome 1 both 0 cM, which IBIS reads as "no genetic map" and then rejects the
    file, so each chromosome is instead shifted up until its first marker is at 0. That keeps
    every distance between markers, and so every segment length and the map length, unchanged.

summarise
    Turns the IBIS segment files and the KING .seg into one row per related pair, with the
    relatedness re-estimated from IBIS runs at several length cutoffs. See
    build_relatedness_table.

IBIS reports one continuous stretch of IBD as several abutting segments: it reports IBD2 stretches
as separate segments, so every IBD1/IBD2 boundary starts a new one, and now and then two IBD1
segments abut with no gap. A parent-offspring pair, which cannot be IBD2, gets a few false short
IBD2 calls that each split an IBD1 stretch in three. The summary therefore joins the segments of a
pair on a chromosome into runs (join_into_runs): consecutive segments, by genetic start, belong to
one run when the gap between the previous segment's end and this one's start is under
--join-gap-cm, whatever their IBD types. The gaps between segments are bimodal (many are exactly
0, very few are 0.5 to 2 cM), so the default of 0.5 cM joins the abutting segments and leaves
separate ones apart. A run is as long as its genetic span, its kinship contribution is the sum over
its segments of length x kinship weight, and the --min-seg-cm floor and every length cutoff apply
to the run, not to its segments.

IBIS .seg layout (README at commit 0e1d9aa, "Segment file format"), tab-separated, no header:

    sample1 sample2 chrom phys_start_pos phys_end_pos IBD_type genetic_start_pos
    genetic_end_pos genetic_seg_length marker_count error_count error_density

IBD_type is IBD1 or IBD2, and sample IDs are `<FID>:<IID>` unless IBIS runs with -noFamID.
"""

import argparse
import gzip
import json
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Final

import numpy as np
import pandas as pd

AUTOSOMES: Final[tuple[str, ...]] = tuple(str(chromosome) for chromosome in range(1, 23))
# IBIS recognises the X chromosome only under this exact name (IBIS.cc passes XchrName "X" to
# its reader), so the X .bim must name it X, not chrX or 23.
X_CHROMOSOME: Final = 'X'

BIM_COLUMNS: Final[tuple[str, ...]] = ('CHROM', 'ID', 'CM', 'POS', 'A1', 'A2')
FAM_COLUMNS: Final[tuple[str, ...]] = ('FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO')
IBIS_SEG_COLUMNS: Final[tuple[str, ...]] = (
    'ID1',
    'ID2',
    'CHROM',
    'PHYS_START',
    'PHYS_END',
    'IBD_TYPE',
    'CM_START',
    'CM_END',
    'LENGTH_CM',
    'N_MARKERS',
    'N_ERRORS',
    'ERROR_DENSITY',
)
IBIS_SEG_DTYPES: Final[dict[str, type]] = {
    'ID1': str,
    'ID2': str,
    'CHROM': str,
    'IBD_TYPE': str,
    'CM_START': float,
    'CM_END': float,
    'LENGTH_CM': float,
}
IBIS_SEG_USED_COLUMNS: Final[tuple[str, ...]] = tuple(IBIS_SEG_DTYPES)

# Kinship weight of one cM of each IBIS segment type: IBD1 shares one allele (kinship 1/4),
# IBD2 shares both (kinship 1/2). KING's IBD1Seg/4 + IBD2Seg/2 and IBIS's seg2coef use the same
# weights; IBIS additionally adds a supplemental 0.00138 to its own coefficients, which is
# not added here.
IBD1_KINSHIP_WEIGHT: Final = 0.25
IBD2_KINSHIP_WEIGHT: Final = 0.5
KINSHIP_WEIGHT_BY_IBD_TYPE: Final[dict[str, float]] = {'IBD1': IBD1_KINSHIP_WEIGHT, 'IBD2': IBD2_KINSHIP_WEIGHT}

# Lower kinship bound of each degree, highest first. Dup/MZ is strictly above its bound; each
# lower degree includes its bound. Below the last bound the pair is unrelated (NA).
# Manichaikul et al. (2010) Bioinformatics 26:2867, doi:10.1093/bioinformatics/btq559:
# kinship > 0.354 duplicate/MZ, [0.177, 0.354] 1st degree, [0.0884, 0.177] 2nd degree,
# [0.0442, 0.0884] 3rd degree.
KING_DEGREE_THRESHOLDS: Final[tuple[tuple[float, str], ...]] = (
    (0.354, 'Dup/MZ'),
    (0.177, '1st'),
    (0.0884, '2nd'),
    (0.0442, '3rd'),
)

# KING InfType values in a --ibdseg .seg file.
KING_INFTYPES: Final[frozenset[str]] = frozenset({'Dup/MZ', 'PO', 'FS', '2nd', '3rd', '4th'})

# Fixed segment-length thresholds for the IBIS_N_SEG_GE<N>CM columns.
N_SEG_THRESHOLDS_CM: Final[tuple[int, ...]] = (7, 10)

PAIR_KEYS: Final[list[str]] = ['ID1', 'ID2']

# Columns of the runs join_into_runs returns: KINSHIP_CM is the sum over the run's segments of
# length x kinship weight, so kinship is its sum divided by the map length.
RUN_COLUMNS: Final[tuple[str, ...]] = (*PAIR_KEYS, 'CHROM', 'CM_START', 'CM_END', 'LENGTH_CM', 'KINSHIP_CM')

# IBIS looks at this many leading markers to decide whether the .bim has a genetic map.
IBIS_MAP_PROBE_MARKERS: Final = 2

# IBIS prints float32 sums, so its reported map length differs from ours in the 6th digit.
IBIS_MAP_LENGTH_REL_TOL: Final = 1e-4
IBIS_MAP_LENGTH_REGEX: Final = re.compile(r'Total Genetic Length in use:\s*(\S+)')

FID_PREFIX_HINT: Final = 'IBIS prefixes sample IDs with <FID>:; check it ran without -noFamID'


def read_bim(path: str | Path) -> pd.DataFrame:
    """Read a PLINK .bim file into the columns CHROM, ID, CM, POS, A1, A2."""
    return pd.read_csv(
        path,
        sep=r'\s+',
        header=None,
        names=list(BIM_COLUMNS),
        dtype={'CHROM': str, 'ID': str, 'CM': float, 'POS': int, 'A1': str, 'A2': str},
        keep_default_na=False,
    )


def write_bim(bim: pd.DataFrame, path: str | Path) -> None:
    """Write a .bim in PLINK's tab-separated layout, with cM to nine significant digits."""
    bim.to_csv(path, sep='\t', header=False, index=False, float_format='%.9g')


def check_chromosome_blocks(bim: pd.DataFrame, expected_chromosomes: Sequence[str]) -> None:
    """Raise unless the .bim holds exactly the expected chromosomes, one contiguous block each, in order.

    IBIS treats a chromosome name that reappears later as a new chromosome, and recognises X only
    under the name X, so a .bim with chr-prefixed or numeric (23) names must fail here rather than
    produce a silently wrong run.
    """
    chromosome = bim['CHROM']
    observed = chromosome[chromosome != chromosome.shift()].tolist()
    if observed != list(expected_chromosomes):
        raise ValueError(
            f'.bim chromosome blocks are {observed}, expected {list(expected_chromosomes)} '
            '(each chromosome once, contiguous, in this order)',
        )


def cm_shift_by_chromosome(bim: pd.DataFrame) -> pd.Series:
    """How far up each chromosome must move so its first variant is at least 0 cM (0 if it already is)."""
    first = bim.groupby('CHROM', sort=False)['CM'].first()
    return -first.clip(upper=0.0)


def shift_to_non_negative_cm(bim: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with each chromosome shifted up so its first variant is at 0 cM, if it was below 0.

    plink2 --cm-map extrapolates variants before a chromosome's first map row, which gives them a
    slightly negative cM. Clamping those to 0 is not an option: IBIS (genetio marker.cc) treats the
    map as missing when the first two markers of the file are both at 0 cM, and then fails on any
    later marker with a non-zero cM, which is what clamping gives chromosome 1. Shifting keeps the
    cM distance between every pair of markers, so segment lengths and the map length (last minus
    first cM per chromosome) are exactly what the unshifted map gives.
    """
    return bim.assign(CM=bim['CM'] + bim['CHROM'].map(cm_shift_by_chromosome(bim)))


def check_cm_non_decreasing(bim: pd.DataFrame) -> None:
    """Raise if cM decreases between consecutive variants of a chromosome."""
    decreases = bim.groupby('CHROM', sort=False)['CM'].diff() < 0
    if decreases.any():
        first = int(decreases.idxmax())
        raise ValueError(
            f'.bim cM decreases at variant line {first + 1} on chromosome {bim.loc[first, "CHROM"]} '
            f'({decreases.sum()} decreases in total)',
        )


def check_map_applied(bim: pd.DataFrame) -> None:
    """Raise if any chromosome spans 0 cM, meaning the genetic map was not applied to it.

    PLINK 1.9 --cm-map with an '@' pattern silently leaves chromosome X at 0 cM because it looks
    for a file named for X's numeric code (23); this catches that class of failure.
    """
    cm = bim.groupby('CHROM', sort=False)['CM']
    spans = cm.last() - cm.first()
    flat = spans[spans <= 0].index.tolist()
    if flat:
        raise ValueError(f'.bim chromosome(s) {flat} span 0 cM: the genetic map was not applied')


def check_ibis_reads_map(bim: pd.DataFrame) -> None:
    """Raise if IBIS would take the genetic map for missing.

    IBIS treats the map as absent, and sets cM from physical position, when the first two markers in
    the file are both at 0 cM (float32); it then exits with an error at the first marker whose cM
    is not 0.
    """
    first_markers = bim['CM'].iloc[:IBIS_MAP_PROBE_MARKERS].to_numpy(dtype=np.float32)
    if len(first_markers) == IBIS_MAP_PROBE_MARKERS and (first_markers == 0).all():
        raise ValueError(
            'the first two markers are both at 0 cM, which IBIS reads as "no genetic map" and then rejects'
        )


def map_length_cm(bim: pd.DataFrame) -> float:
    """Total map length: the sum over chromosomes of the last minus the first variant's cM.

    This is what IBIS prints as "Total Genetic Length in use" and divides by to get kinship.
    """
    cm = bim.groupby('CHROM', sort=False)['CM']
    return float((cm.last() - cm.first()).sum())


def finalise_bim(bim: pd.DataFrame, expected_chromosomes: Sequence[str]) -> tuple[pd.DataFrame, float]:
    """Check the chromosome blocks, make cM non-negative and check the map; return the .bim and its map length."""
    check_chromosome_blocks(bim, expected_chromosomes)
    shifted = shift_to_non_negative_cm(bim)
    check_cm_non_decreasing(shifted)
    check_map_applied(shifted)
    check_ibis_reads_map(shifted)
    return shifted, map_length_cm(shifted)


def read_fam(path: str | Path) -> pd.DataFrame:
    """Read the FID, IID and SEX columns of a PLINK .fam file."""
    return pd.read_csv(
        path,
        sep=r'\s+',
        header=None,
        names=list(FAM_COLUMNS),
        usecols=['FID', 'IID', 'SEX'],
        dtype=str,
        keep_default_na=False,
    )


def single_fid(fam: pd.DataFrame) -> str:
    """The one family ID shared by every sample in the .fam; raise if there is more than one."""
    fids = sorted(fam['FID'].unique())
    if len(fids) != 1:
        raise ValueError(f'expected one FID in the .fam, found {len(fids)}: {fids[:5]}')
    return fids[0]


def strip_fid_prefix(ids: pd.Series, fid: str) -> pd.Series:
    """Remove the `<fid>:` prefix IBIS puts on sample IDs, raising if any ID lacks it."""
    prefix = f'{fid}:'
    missing_prefix = ~ids.str.startswith(prefix)
    if missing_prefix.any():
        raise ValueError(
            f'{missing_prefix.sum()} IBIS sample IDs do not start with "{prefix}" '
            f'(e.g. {ids[missing_prefix].iloc[0]!r}); {FID_PREFIX_HINT}',
        )
    return ids.str.slice(len(prefix))


def check_ids_known(ids: pd.Series, known_ids: set[str], label: str) -> None:
    """Raise if any ID is not in the known set (the IIDs of the samples IBIS and KING were run on)."""
    unknown = sorted(set(ids.unique()) - known_ids)
    if unknown:
        raise ValueError(f'{len(unknown)} {label} sample IDs are not in the .fam (e.g. {unknown[:3]})')


def _read_first_line(path: str | Path) -> str:
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        return handle.readline()


def read_ibis_segments(path: str | Path) -> pd.DataFrame:
    """Read an IBIS .seg or .seg.gz into ID1, ID2, CHROM, IBD_TYPE, CM_START, CM_END and LENGTH_CM.

    IDs are returned as IBIS wrote them (`<FID>:<IID>`); see strip_fid_prefix. An empty file, which
    IBIS writes when it finds no segments, gives an empty frame.

    Raises:
        ValueError: if the first row does not have the 12 columns of the IBIS layout, or a
            segment is neither IBD1 nor IBD2.
    """
    first_line = _read_first_line(path)
    if not first_line.strip():
        return pd.DataFrame({column: pd.Series(dtype=dtype) for column, dtype in IBIS_SEG_DTYPES.items()})
    n_fields = len(first_line.rstrip('\n').split('\t'))
    if n_fields != len(IBIS_SEG_COLUMNS):
        raise ValueError(f'{path}: IBIS .seg rows have {n_fields} columns, expected {len(IBIS_SEG_COLUMNS)}')

    segments = pd.read_csv(
        path,
        sep='\t',
        header=None,
        names=list(IBIS_SEG_COLUMNS),
        usecols=list(IBIS_SEG_USED_COLUMNS),
        dtype=IBIS_SEG_DTYPES,
        keep_default_na=False,
    )
    unknown_types = set(segments['IBD_TYPE']) - set(KINSHIP_WEIGHT_BY_IBD_TYPE)
    if unknown_types:
        raise ValueError(f'{path}: unexpected IBIS segment type(s) {sorted(unknown_types)}')
    return segments


def read_king_seg(path: str | Path) -> pd.DataFrame:
    """Read a KING --ibdseg .seg file into ID1, ID2, IBD1Seg, IBD2Seg and InfType.

    Raises:
        ValueError: if a required column is missing or an InfType is not one KING emits.
    """
    king = pd.read_csv(path, sep=r'\s+', dtype={'ID1': str, 'ID2': str}, keep_default_na=False)
    required = ['ID1', 'ID2', 'IBD1Seg', 'IBD2Seg', 'InfType']
    missing = [column for column in required if column not in king.columns]
    if missing:
        raise ValueError(f'{path}: KING .seg is missing column(s) {missing}')
    unknown = set(king['InfType']) - KING_INFTYPES
    if unknown:
        raise ValueError(f'{path}: unrecognised KING InfType value(s) {sorted(unknown)}')
    return king[required]


def parse_ibis_map_length(log_text: str) -> float:
    """The "Total Genetic Length in use" value from an IBIS log."""
    match = IBIS_MAP_LENGTH_REGEX.search(log_text)
    if match is None:
        raise ValueError('no "Total Genetic Length in use" line in the IBIS log')
    return float(match.group(1))


def check_map_length_matches_ibis(expected_cm: float, ibis_cm: float) -> None:
    """Raise unless the .bim map length matches IBIS's own, to within float32 rounding."""
    if abs(expected_cm - ibis_cm) > IBIS_MAP_LENGTH_REL_TOL * expected_cm:
        raise ValueError(
            f'map length from the corrected .bim ({expected_cm:.4f} cM) differs from the length IBIS '
            f'used ({ibis_cm:.4f} cM)',
        )


def king_kinship(ibd1_seg: pd.Series, ibd2_seg: pd.Series) -> pd.Series:
    """KING kinship from the .seg IBD proportions: IBD1Seg/4 + IBD2Seg/2 (as in merge_qc)."""
    return ibd1_seg.astype(float) * IBD1_KINSHIP_WEIGHT + ibd2_seg.astype(float) * IBD2_KINSHIP_WEIGHT


def kinship_to_degree(kinship: pd.Series) -> pd.Series:
    """Degree label for each kinship from KING's thresholds (KING_DEGREE_THRESHOLDS); NA if unrelated or NA."""
    degree = pd.Series(pd.NA, index=kinship.index, dtype='string')
    for position, (lower_bound, label) in enumerate(KING_DEGREE_THRESHOLDS):
        in_range = kinship > lower_bound if position == 0 else kinship >= lower_bound
        degree = degree.mask(in_range & degree.isna(), label)
    return degree


def normalise_pair_order(pairs: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with ID1 < ID2 on every row, so a pair has one key whichever way a tool listed it."""
    if (pairs['ID1'] == pairs['ID2']).any():
        raise ValueError('a pair lists the same sample twice')
    swap = pairs['ID1'] > pairs['ID2']
    normalised = pairs.copy()
    normalised.loc[swap, 'ID1'] = pairs.loc[swap, 'ID2']
    normalised.loc[swap, 'ID2'] = pairs.loc[swap, 'ID1']
    return normalised


def cutoff_label(cutoff_cm: float) -> str:
    """The cutoff as it appears in column names: 5 -> '5', 7.5 -> '7.5'."""
    return f'{cutoff_cm:g}'


def join_into_runs(segments: pd.DataFrame, join_gap_cm: float, allow_overlap: bool = False) -> pd.DataFrame:
    """Join the abutting segments of each pair on each chromosome into runs, and normalise pair order.

    Segments of a pair on a chromosome are taken in order of genetic start; a segment joins the
    previous one's run when its start is less than join_gap_cm after the previous end, whatever the
    IBD types (join_gap_cm 0 joins nothing). A run's CM_START and CM_END are its first start and
    last end, LENGTH_CM is the span between them (so it includes any gap inside the run), and
    KINSHIP_CM is the sum over its segments of LENGTH_CM x the kinship weight of the IBD type.
    Returns RUN_COLUMNS, one row per run, with no length floor.

    With allow_overlap, a segment that starts before the end of the run so far joins that run. Use it
    only where KINSHIP_CM is not used: IBIS can call the same X stretch of two males as both IBD1 and
    IBD2 (seen once in 1.9M X segments), which would count it twice.

    Raises:
        ValueError: without allow_overlap, if a segment starts before the previous segment of its pair
            and chromosome ends. IBIS emits no overlapping autosomal segments, so an overlap there
            means the input is malformed.
    """
    chromosome_keys = [*PAIR_KEYS, 'CHROM']
    ordered = (
        normalise_pair_order(segments).sort_values([*chromosome_keys, 'CM_START', 'CM_END']).reset_index(drop=True)
    )
    # The gap is measured from the furthest end so far, which is the previous end unless segments overlap.
    run_end = ordered['CM_END'].astype(float).groupby([ordered[key] for key in chromosome_keys]).cummax()
    gap_cm = ordered['CM_START'] - run_end.groupby([ordered[key] for key in chromosome_keys]).shift()
    overlaps = gap_cm < 0
    if overlaps.any() and not allow_overlap:
        first = int(overlaps.idxmax())
        raise ValueError(
            f'{overlaps.sum()} IBIS segments start before the previous segment of the same pair and chromosome '
            f'ends (first on chromosome {ordered.loc[first, "CHROM"]}, overlapping by {-gap_cm[first]:.6f} cM)',
        )
    # The first segment of each pair and chromosome has no gap (NaN) and so always starts a run.
    starts_run = gap_cm.isna() | (gap_cm >= join_gap_cm)
    ordered = ordered.assign(
        RUN=starts_run.cumsum(),
        KINSHIP_CM=ordered['LENGTH_CM'] * ordered['IBD_TYPE'].map(KINSHIP_WEIGHT_BY_IBD_TYPE),
    )
    runs = (
        ordered.groupby([*chromosome_keys, 'RUN'], sort=False)
        .agg(CM_START=('CM_START', 'min'), CM_END=('CM_END', 'max'), KINSHIP_CM=('KINSHIP_CM', 'sum'))
        .reset_index()
    )
    runs['LENGTH_CM'] = runs['CM_END'] - runs['CM_START']
    return runs[list(RUN_COLUMNS)]


def _prepare_runs(
    segments: pd.DataFrame,
    min_seg_cm: float,
    join_gap_cm: float,
    allow_overlap: bool = False,
) -> pd.DataFrame:
    """Join segments into runs (see join_into_runs) and keep only runs of at least min_seg_cm.

    The floor applies to the run, not to its segments: IBIS reports IBD2 segments down to its
    -min_l2 and keeps any IBD1 piece that abuts an IBD2 segment whatever its length, so a segment
    shorter than min_seg_cm can be part of a long run and still counts towards it.
    """
    required = set(IBIS_SEG_USED_COLUMNS)
    if not required <= set(segments.columns):
        raise ValueError(f'segments are missing column(s) {sorted(required - set(segments.columns))}')
    runs = join_into_runs(segments, join_gap_cm, allow_overlap)
    return runs[runs['LENGTH_CM'] >= min_seg_cm]


def _autosomal_pair_stats(runs: pd.DataFrame, map_length: float, cutoffs_cm: Sequence[float]) -> pd.DataFrame:
    by_pair = runs.groupby(PAIR_KEYS)['LENGTH_CM']
    stats = pd.DataFrame(
        {
            'IBIS_N_SEG': by_pair.size(),
            'IBIS_TOTAL_CM': by_pair.sum(),
            'IBIS_LONGEST_CM': by_pair.max(),
        },
    )
    for threshold in N_SEG_THRESHOLDS_CM:
        long_enough = runs[runs['LENGTH_CM'] >= threshold].groupby(PAIR_KEYS).size()
        stats[f'IBIS_N_SEG_GE{threshold}CM'] = long_enough.reindex(stats.index, fill_value=0)

    for cutoff in cutoffs_cm:
        kept = runs[runs['LENGTH_CM'] >= cutoff].groupby(PAIR_KEYS)['KINSHIP_CM'].sum()
        stats[f'IBIS_KINSHIP_GE{cutoff_label(cutoff)}CM'] = kept.reindex(stats.index, fill_value=0.0) / map_length
    return stats.reset_index()


def _x_pair_stats(runs: pd.DataFrame) -> pd.DataFrame:
    by_pair = runs.groupby(PAIR_KEYS)['LENGTH_CM']
    return pd.DataFrame(
        {
            'X_N_SEG': by_pair.size(),
            'X_TOTAL_CM': by_pair.sum(),
            'X_LONGEST_CM': by_pair.max(),
        },
    ).reset_index()


def _king_pairs(king_seg: pd.DataFrame) -> pd.DataFrame:
    king = normalise_pair_order(
        pd.DataFrame(
            {
                'ID1': king_seg['ID1'],
                'ID2': king_seg['ID2'],
                'KING_INFTYPE': king_seg['InfType'],
                'KING_KINSHIP': king_kinship(king_seg['IBD1Seg'], king_seg['IBD2Seg']),
            },
        ),
    )
    if king.duplicated(PAIR_KEYS).any():
        raise ValueError('the KING .seg lists some pair more than once')
    return king


def ibis_columns(cutoffs_cm: Sequence[float]) -> list[str]:
    """Output column names, in order, for the given cutoffs."""
    per_cutoff = [
        column
        for cutoff in cutoffs_cm
        for column in (f'IBIS_KINSHIP_GE{cutoff_label(cutoff)}CM', f'IBIS_DEGREE_GE{cutoff_label(cutoff)}CM')
    ]
    return [
        *PAIR_KEYS,
        'KING_INFTYPE',
        'KING_KINSHIP',
        *per_cutoff,
        'IBIS_N_SEG',
        'IBIS_TOTAL_CM',
        'IBIS_LONGEST_CM',
        *(f'IBIS_N_SEG_GE{threshold}CM' for threshold in N_SEG_THRESHOLDS_CM),
        'X_N_SEG',
        'X_TOTAL_CM',
        'X_LONGEST_CM',
    ]


def build_relatedness_table(
    autosomal_segments: pd.DataFrame,
    x_segments: pd.DataFrame,
    king_seg: pd.DataFrame,
    autosomal_map_length_cm: float,
    cutoffs_cm: Sequence[float],
    min_seg_cm: float,
    join_gap_cm: float,
) -> pd.DataFrame:
    """One row per pair related by KING or by an IBIS autosomal run, with IBIS kinship at length cutoffs.

    IBIS segments of a pair on a chromosome are joined into runs (see join_into_runs): a segment
    whose start is less than join_gap_cm after the previous segment's end joins that segment's
    run, whatever the IBD types. Every IBIS column describes runs, not segments.

    Args:
        autosomal_segments: IBIS autosomal segments with columns ID1, ID2, CHROM, IBD_TYPE
            (IBD1 or IBD2), CM_START, CM_END and LENGTH_CM, with sample IDs matching KING's (FID
            prefix stripped).
        x_segments: the same for the separate IBIS X chromosome run.
        king_seg: the KING .seg with ID1, ID2, IBD1Seg, IBD2Seg and InfType.
        autosomal_map_length_cm: total autosomal genetic map length the kinships are divided by.
        cutoffs_cm: run-length cutoffs (cM); each must be at least min_seg_cm.
        min_seg_cm: the -min_l IBIS ran with. Runs shorter than this are ignored throughout (see
            _prepare_runs).
        join_gap_cm: segments closer than this (cM) join the same run; 0 joins nothing.

    Returns:
        One row per pair, ID1 < ID2, sorted, with the columns of ibis_columns(cutoffs_cm). A run is
        as long as its genetic span (last end minus first start), and its kinship contribution is
        the sum over its segments of IBD1 cM/4 and IBD2 cM/2:

        - KING_INFTYPE, KING_KINSHIP (IBD1Seg/4 + IBD2Seg/2); NA for pairs KING did not call.
        - IBIS_KINSHIP_GE<L>CM: the kinship contributions of the autosomal runs of at least L cM,
          summed and divided by the autosomal map length. 0 if the pair has IBIS runs but none
          that long. IBIS_DEGREE_GE<L>CM: the degree for that kinship from KING's thresholds
          (KING_DEGREE_THRESHOLDS), NA if below the 3rd-degree bound.
        - IBIS_N_SEG, IBIS_TOTAL_CM (sum of run lengths), IBIS_LONGEST_CM, IBIS_N_SEG_GE7CM,
          IBIS_N_SEG_GE10CM over all autosomal runs of at least min_seg_cm.
        - X_N_SEG, X_TOTAL_CM, X_LONGEST_CM over the X runs, joined the same way. X is never mixed
          into the autosomal kinship, and a pair whose only IBIS runs are on X is not listed.

        Every IBIS column is NA for a pair with no IBIS autosomal run (X columns: no X run). IBIS
        compares every pair, so for a KING-called pair that means no run of at least min_seg_cm,
        not missing data.

    Raises:
        ValueError: if the map length is not positive, a cutoff is below min_seg_cm or repeated,
            join_gap_cm is negative, segments overlap, or the inputs are malformed.
    """
    if autosomal_map_length_cm <= 0:
        raise ValueError(f'autosomal map length must be positive, got {autosomal_map_length_cm}')
    cutoffs = list(cutoffs_cm)
    if not cutoffs:
        raise ValueError('at least one cutoff is required')
    below_floor = [cutoff for cutoff in cutoffs if cutoff < min_seg_cm]
    if below_floor:
        raise ValueError(f'cutoffs {below_floor} are below min_seg_cm={min_seg_cm}, which IBIS never reports')
    if len({cutoff_label(cutoff) for cutoff in cutoffs}) != len(cutoffs):
        raise ValueError(f'cutoffs {cutoffs} repeat a value')
    if join_gap_cm < 0:
        raise ValueError(f'join_gap_cm must not be negative, got {join_gap_cm}')

    king = _king_pairs(king_seg)
    autosomal = _autosomal_pair_stats(
        _prepare_runs(autosomal_segments, min_seg_cm, join_gap_cm),
        autosomal_map_length_cm,
        cutoffs,
    )
    # X runs carry no kinship, so an overlapping IBD1/IBD2 double call there is joined, not rejected.
    x_stats = _x_pair_stats(_prepare_runs(x_segments, min_seg_cm, join_gap_cm, allow_overlap=True))

    pairs = pd.concat([king[PAIR_KEYS], autosomal[PAIR_KEYS]]).drop_duplicates()
    table = pairs.merge(king, how='left', on=PAIR_KEYS)
    table = table.merge(autosomal, how='left', on=PAIR_KEYS)
    table = table.merge(x_stats, how='left', on=PAIR_KEYS)

    for cutoff in cutoffs:
        label = cutoff_label(cutoff)
        table[f'IBIS_DEGREE_GE{label}CM'] = kinship_to_degree(table[f'IBIS_KINSHIP_GE{label}CM'])

    count_columns = [column for column in table.columns if '_N_SEG' in column]
    table[count_columns] = table[count_columns].astype('Int64')
    return table[ibis_columns(cutoffs)].sort_values(PAIR_KEYS).reset_index(drop=True)


def table_counts(table: pd.DataFrame, cutoffs_cm: Sequence[float]) -> dict[str, Any]:
    """Aggregate counts describing a built table, for the sidecar JSON."""
    has_king = table['KING_INFTYPE'].notna()
    has_ibis = table['IBIS_N_SEG'].notna()
    return {
        'n_pairs': len(table),
        'n_king_pairs': int(has_king.sum()),
        'n_ibis_pairs': int(has_ibis.sum()),
        'n_king_pairs_without_ibis_segments': int((has_king & ~has_ibis).sum()),
        'n_pairs_with_x_segments': int(table['X_N_SEG'].notna().sum()),
        'n_pairs_by_ibis_degree': {
            f'GE{cutoff_label(cutoff)}CM': {
                str(degree): int(count)
                for degree, count in table[f'IBIS_DEGREE_GE{cutoff_label(cutoff)}CM'].value_counts().items()
            }
            for cutoff in cutoffs_cm
        },
    }


def parse_cutoffs(text: str) -> list[float]:
    """Parse a comma-separated list of cutoffs such as '5,10,15'."""
    return [float(part) for part in text.split(',')]


def prepare_maps(args: argparse.Namespace) -> None:
    """Step prepare-maps: make cM non-negative in the autosomal and X .bim files, check them and write the copies."""
    autosomal_in = read_bim(args.auto_bim)
    x_in = read_bim(args.x_bim)
    autosomal, autosomal_length = finalise_bim(autosomal_in, AUTOSOMES)
    x_bim, x_length = finalise_bim(x_in, [X_CHROMOSOME])
    write_bim(autosomal, args.out_auto_bim)
    write_bim(x_bim, args.out_x_bim)

    stats = {
        'autosomal_map_length_cm': autosomal_length,
        'x_map_length_cm': x_length,
        'n_autosomal_variants': len(autosomal),
        'n_x_variants': len(x_bim),
        'n_autosomal_variants_before_map_start': int((autosomal_in['CM'] < 0).sum()),
        'n_x_variants_before_map_start': int((x_in['CM'] < 0).sum()),
        'autosomal_max_cm_shift': float(cm_shift_by_chromosome(autosomal_in).max()),
        'x_cm_shift': float(cm_shift_by_chromosome(x_in).max()),
    }
    Path(args.out_stats).write_text(json.dumps(stats, indent=2) + '\n')
    print(json.dumps(stats, indent=2))


def summarise(args: argparse.Namespace) -> None:
    """Step summarise: build the pairs table and its sidecar JSON from the IBIS and KING outputs."""
    cutoffs = parse_cutoffs(args.cutoffs_cm)
    stats = json.loads(Path(args.map_stats).read_text())
    map_length = stats['autosomal_map_length_cm']

    ibis_reported_length = parse_ibis_map_length(Path(args.ibis_auto_log).read_text())
    check_map_length_matches_ibis(map_length, ibis_reported_length)

    fam = read_fam(args.fam)
    fam_x = read_fam(args.fam_x)
    fid = single_fid(fam)
    iids = set(fam['IID'])
    iids_x = set(fam_x['IID'])
    if not iids_x <= iids:
        raise ValueError('the X .fam holds samples that are not in the autosomal .fam')

    segments = {}
    for name, path, known_ids in (
        ('autosomal', args.ibis_auto_seg, iids),
        ('x', args.ibis_x_seg, iids_x),
    ):
        raw = read_ibis_segments(path)
        for column in ('ID1', 'ID2'):
            raw[column] = strip_fid_prefix(raw[column], fid)
            check_ids_known(raw[column], known_ids, f'IBIS {name}')
        segments[name] = raw

    king_seg = read_king_seg(args.king_seg)
    for column in ('ID1', 'ID2'):
        check_ids_known(king_seg[column], iids, 'KING')

    table = build_relatedness_table(
        segments['autosomal'],
        segments['x'],
        king_seg,
        map_length,
        cutoffs,
        args.min_seg_cm,
        args.join_gap_cm,
    )
    table.to_csv(args.out_table, sep='\t', index=False, na_rep='NA', float_format='%.6f', compression='gzip')

    sidecar = {
        'run': json.loads(args.run_info),
        'min_seg_cm': args.min_seg_cm,
        'join_gap_cm': args.join_gap_cm,
        'cutoffs_cm': cutoffs,
        'degree_thresholds': [
            {'kinship_lower_bound': bound, 'degree': label} for bound, label in KING_DEGREE_THRESHOLDS
        ],
        'autosomal_map_length_cm': map_length,
        'ibis_reported_autosomal_map_length_cm': ibis_reported_length,
        'x_map_length_cm': stats['x_map_length_cm'],
        'n_samples': len(fam),
        'n_samples_x': len(fam_x),
        'n_autosomal_ibis_segments': len(segments['autosomal']),
        'n_x_ibis_segments': len(segments['x']),
        'map_stats': stats,
        'counts': table_counts(table, cutoffs),
    }
    Path(args.out_sidecar).write_text(json.dumps(sidecar, indent=2) + '\n')
    print(json.dumps(sidecar['counts'], indent=2))


def build_parser() -> argparse.ArgumentParser:
    """The command line of the job-side steps."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    steps = parser.add_subparsers(dest='step', required=True)

    prepare = steps.add_parser('prepare-maps', help=prepare_maps.__doc__)
    prepare.add_argument('--auto-bim', required=True)
    prepare.add_argument('--x-bim', required=True)
    prepare.add_argument('--out-auto-bim', required=True)
    prepare.add_argument('--out-x-bim', required=True)
    prepare.add_argument('--out-stats', required=True)
    prepare.set_defaults(run=prepare_maps)

    summary = steps.add_parser('summarise', help=summarise.__doc__)
    summary.add_argument('--ibis-auto-seg', required=True)
    summary.add_argument('--ibis-x-seg', required=True)
    summary.add_argument('--ibis-auto-log', required=True)
    summary.add_argument('--king-seg', required=True)
    summary.add_argument('--fam', required=True, help='the autosomal .fam IBIS ran on')
    summary.add_argument('--fam-x', required=True, help='the X .fam IBIS ran on')
    summary.add_argument('--map-stats', required=True, help='the JSON written by prepare-maps')
    summary.add_argument('--min-seg-cm', required=True, type=float)
    summary.add_argument(
        '--join-gap-cm',
        required=True,
        type=float,
        help='segments of a pair on a chromosome with a gap below this (cM) are joined into one run',
    )
    summary.add_argument('--cutoffs-cm', required=True)
    summary.add_argument('--run-info', required=True, help='JSON object of run parameters for the sidecar')
    summary.add_argument('--out-table', required=True)
    summary.add_argument('--out-sidecar', required=True)
    summary.set_defaults(run=summarise)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    args.run(args)


if __name__ == '__main__':
    main()
