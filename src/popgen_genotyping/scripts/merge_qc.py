"""Merge PLINK2 QC data files, relatedness estimates, and bafregress results into a CSV summary.

Relatedness comes from three estimators, one column each: KING --ibdseg (``KING_IBDSEG_RELATIVES``),
plink2 KING-robust (``KING_ROBUST_RELATIVES``) and IBIS (``IBIS_RELATIVES_GE5CM``). Each lists a
sample's relatives at 3rd degree or closer as semicolon-separated ``REL_ID:KINSHIP:DEGREE`` strings.

Also folds in the set of samples excluded from all three estimators for contamination or high
missingness (flagged RELATEDNESS_EXCLUDED with the exclusion reason), so they keep their
QC columns without appearing to have relatives.
"""

import argparse
import sys

import pandas as pd

KING_IBDSEG_RELATIVES_COL: str = 'KING_IBDSEG_RELATIVES'
KING_ROBUST_RELATIVES_COL: str = 'KING_ROBUST_RELATIVES'
IBIS_RELATIVES_COL: str = 'IBIS_RELATIVES_GE5CM'
RELATIVES_COLS: list[str] = [KING_IBDSEG_RELATIVES_COL, KING_ROBUST_RELATIVES_COL, IBIS_RELATIVES_COL]

# KING --ibdseg InfType values kept in KING_IBDSEG_RELATIVES, which is also the DEGREE label.
# Ref: Manichaikul et al. (2010), doi:10.1093/bioinformatics/btq559;
# Chen et al. (2024) KING --ibdseg. KING emits InfType in {Dup/MZ, PO, FS, 2nd, 3rd}
# under --degree 3.
KING_IBDSEG_INFTYPES: frozenset[str] = frozenset({'Dup/MZ', 'PO', 'FS', '2nd', '3rd'})

# InfType values beyond the --degree 3 we ask KING for. KING still labels the odd pair
# on the 3rd/4th boundary (PropIBD ~0.088) as 4th; those rows are dropped, not binned.
INFTYPE_BEYOND_DEGREE: frozenset[str] = frozenset({'4th'})

# Lower kinship bound of each degree, highest first. Dup/MZ is strictly above its bound; each
# lower degree includes its bound. Below the last bound the pair is beyond 3rd degree.
# Manichaikul A et al. (2010) Bioinformatics 26:2867-2873, doi:10.1093/bioinformatics/btq559:
# kinship > 0.354 duplicate/MZ, [0.177, 0.354] 1st degree, [0.0884, 0.177] 2nd degree,
# [0.0442, 0.0884] 3rd degree. The same bounds as IBIS_DEGREE_GE5CM, so a pair at a bound is binned
# the same way by KING-robust and IBIS, and every pair kept by plink2's --king-table-filter 0.0442
# (kinship at least that) has a degree.
KING_ROBUST_DEGREE_THRESHOLDS: tuple[tuple[float, str], ...] = (
    (0.354, 'Dup/MZ'),
    (0.177, '1st'),
    (0.0884, '2nd'),
    (0.0442, '3rd'),
)

# Columns of the IBIS relatedness pairs table (ibis_relatedness_summary.py) read here: the kinship
# and degree from IBD runs of at least 5 cM. IBIS_NA is how the table writes a missing value. The
# degree labels are those of KING_ROBUST_DEGREE_THRESHOLDS; a pair below the 3rd-degree bound, or
# without any IBIS run, has no degree (IBIS_NA).
IBIS_KINSHIP_COL: str = 'IBIS_KINSHIP_GE5CM'
IBIS_DEGREE_COL: str = 'IBIS_DEGREE_GE5CM'
IBIS_NA: str = 'NA'
IBIS_DEGREES: frozenset[str] = frozenset({'Dup/MZ', '1st', '2nd', '3rd'})


def read_qc_file(filepath: str) -> pd.DataFrame:
    """Read a whitespace-delimited PLINK2 QC file.

    Strips leading ``#`` from column headers. Compression (``.gz``) is inferred from the
    file name.

    Args:
        filepath: Path to the whitespace-delimited file.

    Returns:
        DataFrame with cleaned column names.
    """
    df = pd.read_csv(filepath, sep=r'\s+', engine='python')
    df.columns = df.columns.str.lstrip('#')
    return df


def format_relatives(pairs: pd.DataFrame, column: str) -> pd.DataFrame:
    """Collapse related pairs into one ``REL_ID:KINSHIP:DEGREE`` string per sample.

    Every pair is listed under both of its samples. A sample's relatives are joined with ``;``
    in the order the pairs appear (the pairs where it is ``ID1`` first, then those where it is
    ``ID2``).

    Args:
        pairs: One row per related pair, with columns ``ID1``, ``ID2``, ``KINSHIP`` and ``DEGREE``.
        column: Name of the output column.

    Returns:
        DataFrame with ``IID`` and ``column``, one row per sample with at least one relative; empty
        with those columns if there are no pairs.
    """
    if pairs.empty:
        return pd.DataFrame(columns=pd.Index(['IID', column]))

    rel_1 = pairs[['ID1', 'ID2', 'KINSHIP', 'DEGREE']].copy()
    rel_1.columns = pd.Index(['IID', 'REL_ID', 'KINSHIP', 'DEGREE'])

    rel_2 = pairs[['ID2', 'ID1', 'KINSHIP', 'DEGREE']].copy()
    rel_2.columns = pd.Index(['IID', 'REL_ID', 'KINSHIP', 'DEGREE'])

    rel_stacked = pd.concat([rel_1, rel_2], ignore_index=True)

    # Format the string payload: REL_ID:KINSHIP:DEGREE
    rel_stacked[column] = (
        rel_stacked['REL_ID'].astype(str)
        + ':'
        + rel_stacked['KINSHIP'].astype(str)
        + ':'
        + rel_stacked['DEGREE'].astype(str)
    )

    return rel_stacked.groupby('IID')[column].apply(';'.join).reset_index()


def kinship_to_degree(kinship: pd.Series) -> pd.Series:
    """Degree label for each KING-robust kinship, from ``KING_ROBUST_DEGREE_THRESHOLDS``.

    Args:
        kinship: Kinship coefficients.

    Returns:
        Series of ``Dup/MZ``, ``1st``, ``2nd`` or ``3rd``, with a missing value for a kinship below
        the 3rd-degree bound.
    """
    degree = pd.Series(pd.NA, index=kinship.index, dtype='object')
    for position, (lower_bound, label) in enumerate(KING_ROBUST_DEGREE_THRESHOLDS):
        in_range = kinship > lower_bound if position == 0 else kinship >= lower_bound
        degree = degree.mask(in_range & degree.isna(), label)
    return degree


def process_seg(seg_path: str) -> pd.DataFrame:
    """Process a KING ``--ibdseg`` ``.seg`` file into the ``KING_IBDSEG_RELATIVES`` column.

    The column contains semicolon-separated ``REL_ID:KINSHIP:DEGREE`` strings. Kinship is the
    IBD-based KING kinship coefficient ``IBD1Seg/4 + IBD2Seg/2`` (equivalent to ``PropIBD/2``),
    to four decimals. DEGREE is KING's relationship inference (``PO``, ``FS``, ``2nd``,
    ``3rd``, ``Dup/MZ``). Pairs labelled beyond the requested degree (``4th``) are dropped,
    with the count printed to stderr; any other InfType raises ``ValueError``.

    Args:
        seg_path: Path to the autosomal ``.seg`` file emitted by KING
            ``--ibdseg``.

    Returns:
        DataFrame with IID and KING_IBDSEG_RELATIVES columns, or an empty
        DataFrame with those columns if no relationships are found.

    Raises:
        ValueError: If the ``.seg`` file is missing any of the required
            columns or contains an unrecognised ``InfType``.
    """
    empty = pd.DataFrame(columns=pd.Index(['IID', KING_IBDSEG_RELATIVES_COL]))

    seg_df = read_qc_file(seg_path)
    if seg_df.empty:
        return empty

    required = {'ID1', 'ID2', 'IBD1Seg', 'IBD2Seg', 'InfType'}
    missing = required - set(seg_df.columns)
    if missing:
        raise ValueError(
            f'.seg file {seg_path} is missing required column(s): {sorted(missing)}',
        )

    unknown = set(seg_df['InfType']) - KING_IBDSEG_INFTYPES - INFTYPE_BEYOND_DEGREE
    if unknown:
        raise ValueError(
            f'.seg file {seg_path} contains unrecognised InfType value(s): {sorted(unknown)}. '
            f'Expected one of {sorted(KING_IBDSEG_INFTYPES)}.',
        )

    beyond = seg_df['InfType'].isin(list(INFTYPE_BEYOND_DEGREE))
    n_beyond = int(beyond.sum())
    if n_beyond:
        print(
            f'Dropping {n_beyond} pair(s) from {seg_path} with InfType beyond the requested degree: '
            f'{sorted(INFTYPE_BEYOND_DEGREE)}',
            file=sys.stderr,
        )
        seg_df = seg_df.loc[~beyond]

    pairs = pd.DataFrame(
        {
            'ID1': seg_df['ID1'],
            'ID2': seg_df['ID2'],
            'KINSHIP': (seg_df['IBD1Seg'].astype(float) / 4.0 + seg_df['IBD2Seg'].astype(float) / 2.0).round(4),
            'DEGREE': seg_df['InfType'],
        },
    )
    return format_relatives(pairs, KING_IBDSEG_RELATIVES_COL)


def process_king_robust(kin0_path: str) -> pd.DataFrame:
    """Process a plink2 ``--make-king-table`` ``.kin0`` file into the ``KING_ROBUST_RELATIVES`` column.

    The column contains semicolon-separated ``REL_ID:KINSHIP:DEGREE`` strings. Kinship is the
    KING-robust estimate in the ``KINSHIP`` column, to four decimals. DEGREE is binned from the
    kinship by ``KING_ROBUST_DEGREE_THRESHOLDS`` (``Dup/MZ``, ``1st``, ``2nd``, ``3rd``). The
    table is written with a kinship floor (plink2 ``--king-table-filter``), so every pair is
    normally within 3rd degree; pairs below the 3rd-degree bound (a floor set lower than that) are
    dropped, with the count printed to stderr.

    Args:
        kin0_path: Path to the ``.kin0`` file (plain or ``.gz``).

    Returns:
        DataFrame with IID and KING_ROBUST_RELATIVES columns, or an empty
        DataFrame with those columns if no relationships are found.

    Raises:
        ValueError: If the ``.kin0`` file is missing any of the required columns or has a
            missing or non-numeric ``KINSHIP``.
    """
    kin0_df = read_qc_file(kin0_path)

    required = {'IID1', 'IID2', 'KINSHIP'}
    missing = required - set(kin0_df.columns)
    if missing:
        raise ValueError(
            f'.kin0 file {kin0_path} is missing required column(s): {sorted(missing)}',
        )
    if kin0_df.empty:
        return pd.DataFrame(columns=pd.Index(['IID', KING_ROBUST_RELATIVES_COL]))

    kinship = pd.Series(pd.to_numeric(kin0_df['KINSHIP'], errors='coerce'))
    if kinship.isna().any():
        raise ValueError(
            f'.kin0 file {kin0_path} has {int(kinship.isna().sum())} missing or non-numeric KINSHIP value(s)'
        )

    degree = kinship_to_degree(kinship)
    beyond = degree.isna()
    n_beyond = int(beyond.sum())
    if n_beyond:
        print(
            f'Dropping {n_beyond} pair(s) from {kin0_path} with kinship below the 3rd-degree bound '
            f'{KING_ROBUST_DEGREE_THRESHOLDS[-1][0]}',
            file=sys.stderr,
        )

    pairs = pd.DataFrame(
        {
            'ID1': kin0_df['IID1'],
            'ID2': kin0_df['IID2'],
            'KINSHIP': kinship.round(4),
            'DEGREE': degree,
        },
    ).loc[~beyond]
    return format_relatives(pairs, KING_ROBUST_RELATIVES_COL)


def process_ibis_pairs(pairs_path: str) -> pd.DataFrame:
    """Process the IBIS relatedness pairs table into the ``IBIS_RELATIVES_GE5CM`` column.

    The column contains semicolon-separated ``REL_ID:KINSHIP:DEGREE`` strings, from the table's
    ``IBIS_KINSHIP_GE5CM`` (kinship from IBD runs of at least 5 cM, to four decimals) and
    ``IBIS_DEGREE_GE5CM`` (``Dup/MZ``, ``1st``, ``2nd``, ``3rd``). A pair whose degree is ``NA``
    (below the 3rd-degree kinship bound, or a KING pair with no IBIS run of that length) is not a
    relative here and is dropped, with the count printed to stderr. Any other degree label raises
    ``ValueError``.

    Args:
        pairs_path: Path to ``ibis_relatedness_pairs.tsv.gz``.

    Returns:
        DataFrame with IID and IBIS_RELATIVES_GE5CM columns, or an empty
        DataFrame with those columns if no relationships are found.

    Raises:
        ValueError: If the table is missing any of the required columns, has an unrecognised
            degree label, or has a related pair with a missing or non-numeric kinship.
    """
    # Strings throughout, and no default NA parsing: a missing value is exactly IBIS_NA.
    table = pd.read_csv(pairs_path, sep='\t', dtype=str, keep_default_na=False)

    required = {'ID1', 'ID2', IBIS_KINSHIP_COL, IBIS_DEGREE_COL}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(
            f'IBIS pairs table {pairs_path} is missing required column(s): {sorted(missing)}',
        )

    no_degree = table[IBIS_DEGREE_COL] == IBIS_NA
    unknown = set(table.loc[~no_degree, IBIS_DEGREE_COL]) - IBIS_DEGREES
    if unknown:
        raise ValueError(
            f'IBIS pairs table {pairs_path} contains unrecognised {IBIS_DEGREE_COL} value(s): {sorted(unknown)}. '
            f'Expected one of {sorted(IBIS_DEGREES)} or {IBIS_NA}.',
        )

    n_no_degree = int(no_degree.sum())
    if n_no_degree:
        print(
            f'Dropping {n_no_degree} pair(s) from {pairs_path} with no {IBIS_DEGREE_COL} '
            '(below the 3rd-degree bound or without an IBIS run of that length)',
            file=sys.stderr,
        )
    related = table.loc[~no_degree]

    try:
        kinship = pd.Series(pd.to_numeric(related[IBIS_KINSHIP_COL], errors='raise'))
    except ValueError as e:
        raise ValueError(
            f'IBIS pairs table {pairs_path} has a related pair with a missing or non-numeric {IBIS_KINSHIP_COL}: {e}',
        ) from e

    pairs = pd.DataFrame(
        {
            'ID1': related['ID1'],
            'ID2': related['ID2'],
            'KINSHIP': kinship.round(4),
            'DEGREE': related[IBIS_DEGREE_COL],
        },
    )
    return format_relatives(pairs, IBIS_RELATIVES_COL)


def process_bafregress(paths: list[str]) -> pd.DataFrame:
    """Read and concatenate bafregress output files.

    Args:
        paths: List of paths to bafregress output files.

    Returns:
        Combined DataFrame with sample_id renamed to IID, or empty DataFrame.
    """
    dfs: list[pd.DataFrame] = []
    for path in paths:
        try:
            df = pd.read_csv(path, sep=r'\s+', engine='python')
            dfs.append(df)
        except (OSError, ValueError) as e:
            print(f'Warning: Could not read bafregress file {path}. Error: {e}')

    if not dfs:
        return pd.DataFrame()

    return pd.concat(dfs, ignore_index=True).rename(columns={'sample_id': 'IID'})


def process_relatedness_excluded(path: str, report_iids: pd.Series) -> pd.DataFrame:
    """Flag every report IID excluded from relatedness estimation with its exclusion reason.

    Args:
        path: Path to the excluded-samples TSV emitted by ``relatedness_exclusions.py``
            (columns ``IID``, ``REASON``, ``BAF_REGRESS``, ``F_MISS``; header-only
            when nothing was excluded).
        report_iids: The full set of IIDs in the QC report.

    Returns:
        DataFrame with ``IID`` and ``RELATEDNESS_EXCLUDED`` for every IID in
        ``report_iids``: the ``REASON`` string (``contamination``,
        ``missingness``, or ``contamination;missingness``) for excluded
        samples, ``NA`` otherwise.

    Raises:
        ValueError: If an excluded IID is not present in ``report_iids``.
    """
    excluded_df = read_qc_file(path)

    missing = sorted(set(excluded_df['IID']) - set(report_iids))
    if missing:
        raise ValueError(f'Relatedness-excluded sample(s) not found in QC report: {missing}')

    reason_by_iid = excluded_df.set_index('IID')['REASON']
    return pd.DataFrame(
        {
            'IID': report_iids,
            'RELATEDNESS_EXCLUDED': report_iids.map(reason_by_iid),
        },
    )


def main() -> None:
    """Parse arguments and run the QC merge pipeline."""
    parser = argparse.ArgumentParser(
        description=(
            'Merge PLINK2 QC files, relatedness from KING --ibdseg, KING-robust and IBIS (plus the '
            'contamination/missingness-based exclusion list), and bafregress into a CSV.'
        ),
    )
    parser.add_argument('--sexcheck', required=True, help='Path to .sexcheck file')
    parser.add_argument('--het', required=True, help='Path to .het file')
    parser.add_argument('--smiss', required=True, help='Path to .smiss file')
    parser.add_argument('--seg', required=True, help='Path to KING --ibdseg autosomal .seg file')
    parser.add_argument('--kin0', required=True, help='Path to the plink2 KING-robust .kin0 file (plain or .gz)')
    parser.add_argument('--ibis-pairs', required=True, help='Path to the IBIS relatedness pairs table (.tsv.gz)')
    parser.add_argument(
        '--relatedness-excluded',
        required=True,
        help='Path to the TSV of samples excluded from relatedness estimation for contamination or missingness',
    )
    parser.add_argument('--output', required=True, help='Output CSV path')
    parser.add_argument(
        '--bafregress',
        nargs='*',
        default=[],
        help='Paths to bafregress output files',
    )
    args = parser.parse_args()

    # Read PLINK2 QC files
    try:
        smiss = read_qc_file(args.smiss)
        het = read_qc_file(args.het)
        sexcheck = read_qc_file(args.sexcheck)
    except (OSError, ValueError) as e:
        print(f'Error reading QC files: {e}')
        sys.exit(1)

    # Identify merge keys
    merge_cols = [col for col in ['FID', 'IID'] if col in smiss.columns]

    # Merge QC dataframes
    merged = smiss.merge(het, on=merge_cols, suffixes=('', '_het'))
    merged = merged.merge(sexcheck, on=merge_cols, suffixes=('', '_sex'))

    # Process the three relatedness estimators. Errors from process_seg, process_king_robust and
    # process_ibis_pairs propagate so the Batch job fails fast on a malformed input rather than
    # emitting a silently-incomplete QC report.
    merged = merged.merge(process_seg(args.seg), on='IID', how='left')
    merged = merged.merge(process_king_robust(args.kin0), on='IID', how='left')
    merged = merged.merge(process_ibis_pairs(args.ibis_pairs), on='IID', how='left')

    # Flag samples excluded from relatedness estimation for contamination or missingness. No
    # estimator saw them, so they must not carry any relatedness — a non-null relatives value
    # here would mean the exclusion list and that estimator's output disagree about which
    # samples went in.
    relatedness_excluded = process_relatedness_excluded(args.relatedness_excluded, merged['IID'])
    merged = merged.merge(relatedness_excluded, on='IID', how='left')

    excluded_rows = merged[merged['RELATEDNESS_EXCLUDED'].notna()]
    for relatives_col in RELATIVES_COLS:
        has_relatives = excluded_rows[relatives_col].notna()
        if has_relatives.any():
            bad_iids = sorted(excluded_rows.loc[has_relatives, 'IID'])
            raise ValueError(
                f'Relatedness-excluded sample(s) unexpectedly have relatives in {relatives_col}: {bad_iids}'
            )

    # Place RELATEDNESS_EXCLUDED immediately after the last relatives column.
    cols = [c for c in merged.columns if c != 'RELATEDNESS_EXCLUDED']
    cols.insert(cols.index(RELATIVES_COLS[-1]) + 1, 'RELATEDNESS_EXCLUDED')
    merged = merged[cols]

    # Process bafregress files
    if args.bafregress:
        baf_df = process_bafregress(args.bafregress)
        if not baf_df.empty:
            merged = merged.merge(baf_df, on='IID', how='left', suffixes=('', '_baf'))

    merged.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
