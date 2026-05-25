"""
Streaming aggregator for per-cohort per-SNP cluster statistics.

Consumes the output of:
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID[\\t%GT,%THETA,%R]\\n'

Emits a bgzip-compressed TSV with one row per (variant, sex_subset).
Autosomes: one row per variant (sex_subset='all').
chrX / chrY / chrM(T): three rows per variant (all / female / male).
"""

from __future__ import annotations

import argparse
import io
import subprocess
import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

CLUSTERS: tuple[str, ...] = ('AA', 'AB', 'BB')
SUBSETS_ALL_ONLY: tuple[str, ...] = ('all',)
SUBSETS_STRATIFIED: tuple[str, ...] = ('all', 'female', 'male')
SEX_STRAT_CHROMS: frozenset[str] = frozenset(
    {'X', 'Y', 'M', 'MT', 'chrX', 'chrY', 'chrM', 'chrMT'},
)

GT_TO_CLUSTER: dict[str, str] = {
    '0/0': 'AA',
    '0|0': 'AA',
    '0/1': 'AB',
    '1/0': 'AB',
    '0|1': 'AB',
    '1|0': 'AB',
    '1/1': 'BB',
    '1|1': 'BB',
}

HEADER_COLS: tuple[str, ...] = (
    'chrom',
    'pos',
    'ref',
    'alt',
    'variant_id',
    'sex_subset',
    'n_AA',
    'n_AB',
    'n_BB',
    'n_nocall',
    'mean_THETA_AA',
    'var_THETA_AA',
    'mean_R_AA',
    'var_R_AA',
    'mean_THETA_AB',
    'var_THETA_AB',
    'mean_R_AB',
    'var_R_AB',
    'mean_THETA_BB',
    'var_THETA_BB',
    'mean_R_BB',
    'var_R_BB',
)

MIN_VARIANT_FIELDS: int = 5
MIN_OBS_FOR_VARIANCE: int = 2
NA_TOKEN: str = 'NA'  # noqa: S105
FLOAT_FMT: str = '.6g'


def load_sample_order(path: Path) -> list[str]:
    """
    Load BCF sample IDs in column order.

    Args:
        path: TSV/text file with one sample ID per line.

    Returns:
        Sample IDs in declaration order; blank lines are skipped.
    """
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def load_sex_map(path: Path) -> dict[str, str]:
    """
    Load the sample-id to sex mapping.

    Args:
        path: Two-column TSV `sample_id<TAB>sex_code` where 1=male, 2=female.
              Other codes are treated as unknown (not stratified).

    Returns:
        Mapping from sample ID to 'M' or 'F'. Samples with unknown sex
        are omitted from the mapping.
    """
    mapping: dict[str, str] = {}
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        sample_id = parts[0]
        sex_code = parts[1] if len(parts) > 1 else '0'
        if sex_code == '1':
            mapping[sample_id] = 'M'
        elif sex_code == '2':
            mapping[sample_id] = 'F'
    return mapping


def build_subset_indices(
    sample_order: Sequence[str],
    sex_map: dict[str, str],
) -> dict[str, list[int]]:
    """
    Compute the column-index list for each sex subset.

    Args:
        sample_order: Sample IDs in BCF column order.
        sex_map: Mapping from sample ID to 'M'/'F' (unknowns absent).

    Returns:
        Mapping subset -> sorted list of column indices into `sample_order`.
        Subsets: 'all', 'female', 'male'.
    """
    return {
        'all': list(range(len(sample_order))),
        'female': [i for i, s in enumerate(sample_order) if sex_map.get(s) == 'F'],
        'male': [i for i, s in enumerate(sample_order) if sex_map.get(s) == 'M'],
    }


def make_cell() -> dict[str, float]:
    """Return a fresh per-(cluster, subset) accumulator."""
    return {
        'n_geno': 0.0,
        'n_theta': 0.0,
        'sum_theta': 0.0,
        'sumsq_theta': 0.0,
        'n_r': 0.0,
        'sum_r': 0.0,
        'sumsq_r': 0.0,
    }


def update_cell(cell: dict[str, float], theta_s: str, r_s: str) -> None:
    """
    Fold one sample's THETA and R into the accumulator.

    Args:
        cell: Mutable accumulator returned by `make_cell`.
        theta_s: Raw THETA token; '.' or unparseable values are skipped.
        r_s: Raw R token; '.' or unparseable values are skipped.

    Returns:
        None; mutates `cell` in place.
    """
    cell['n_geno'] += 1
    try:
        t = float(theta_s)
    except ValueError:
        pass
    else:
        cell['n_theta'] += 1
        cell['sum_theta'] += t
        cell['sumsq_theta'] += t * t
    try:
        r = float(r_s)
    except ValueError:
        pass
    else:
        cell['n_r'] += 1
        cell['sum_r'] += r
        cell['sumsq_r'] += r * r


def moments(
    cell: dict[str, float],
    key_n: str,
    key_sum: str,
    key_sumsq: str,
) -> tuple[str, str]:
    """
    Derive (mean, sample variance) string tokens from a moment triple.

    Returns:
        Tuple of (mean_str, var_str). When n==0 both are 'NA';
        when n==1 mean is reported and variance is 'NA'. Tiny negative
        variances from float roundoff are clipped to 0.
    """
    n = cell[key_n]
    if n == 0:
        return NA_TOKEN, NA_TOKEN
    mean = cell[key_sum] / n
    if n < MIN_OBS_FOR_VARIANCE:
        return format(mean, FLOAT_FMT), NA_TOKEN
    var = (cell[key_sumsq] - cell[key_sum] * mean) / (n - 1)
    if var < 0:
        var = 0.0
    return format(mean, FLOAT_FMT), format(var, FLOAT_FMT)


def emit_row(
    writer: io.TextIOBase,
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    variant_id: str,
    subset: str,
    cells: dict[tuple[str, str], dict[str, float]],
    nocall_n: int,
) -> None:
    """Write one TSV row for a single (variant, subset) combination."""
    aa = cells[('AA', subset)]
    ab = cells[('AB', subset)]
    bb = cells[('BB', subset)]

    aa_theta_m, aa_theta_v = moments(aa, 'n_theta', 'sum_theta', 'sumsq_theta')
    aa_r_m, aa_r_v = moments(aa, 'n_r', 'sum_r', 'sumsq_r')
    ab_theta_m, ab_theta_v = moments(ab, 'n_theta', 'sum_theta', 'sumsq_theta')
    ab_r_m, ab_r_v = moments(ab, 'n_r', 'sum_r', 'sumsq_r')
    bb_theta_m, bb_theta_v = moments(bb, 'n_theta', 'sum_theta', 'sumsq_theta')
    bb_r_m, bb_r_v = moments(bb, 'n_r', 'sum_r', 'sumsq_r')

    writer.write(
        '\t'.join(
            (
                chrom,
                pos,
                ref,
                alt,
                variant_id,
                subset,
                str(int(aa['n_geno'])),
                str(int(ab['n_geno'])),
                str(int(bb['n_geno'])),
                str(nocall_n),
                aa_theta_m,
                aa_theta_v,
                aa_r_m,
                aa_r_v,
                ab_theta_m,
                ab_theta_v,
                ab_r_m,
                ab_r_v,
                bb_theta_m,
                bb_theta_v,
                bb_r_m,
                bb_r_v,
            ),
        ),
    )
    writer.write('\n')


def process_variant(
    line: str,
    sample_indices: dict[str, list[int]],
    writer: io.TextIOBase,
) -> None:
    """
    Aggregate one bcftools-query line and emit one or three output rows.

    Args:
        line: Raw query output. Tab-separated; first five fields are
              chrom/pos/ref/alt/id, remainder is one `GT,THETA,R` field
              per sample in BCF column order.
        sample_indices: Output of `build_subset_indices`.
        writer: Text-mode handle to the TSV stream.

    Returns:
        None. Lines with fewer than five tab-separated fields are skipped.
    """
    fields = line.rstrip('\n').split('\t')
    if len(fields) < MIN_VARIANT_FIELDS:
        return
    chrom, pos, ref, alt, variant_id = fields[0:5]
    per_sample = fields[5:]

    subsets_to_emit = SUBSETS_STRATIFIED if chrom in SEX_STRAT_CHROMS else SUBSETS_ALL_ONLY

    cells: dict[tuple[str, str], dict[str, float]] = {
        (cluster, subset): make_cell() for cluster in CLUSTERS for subset in subsets_to_emit
    }
    nocall_count: dict[str, int] = dict.fromkeys(subsets_to_emit, 0)

    # Parse once; reuse across subsets.
    parsed: list[tuple[str | None, str, str]] = []
    for s in per_sample:
        try:
            gt, theta_s, r_s = s.split(',', 2)
        except ValueError:
            parsed.append((None, '.', '.'))
            continue
        parsed.append((GT_TO_CLUSTER.get(gt), theta_s, r_s))

    for subset in subsets_to_emit:
        for col_idx in sample_indices[subset]:
            cluster, theta_s, r_s = parsed[col_idx]
            if cluster is None:
                nocall_count[subset] += 1
            else:
                update_cell(cells[(cluster, subset)], theta_s, r_s)

    for subset in subsets_to_emit:
        emit_row(writer, chrom, pos, ref, alt, variant_id, subset, cells, nocall_count[subset])


def run(
    input_lines: Iterable[str],
    sample_indices: dict[str, list[int]],
    writer: io.TextIOBase,
) -> None:
    """Write the header and stream every input line through `process_variant`."""
    writer.write('\t'.join(HEADER_COLS) + '\n')
    for line in input_lines:
        process_variant(line, sample_indices, writer)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point. Returns a process exit code."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--sample-order',
        required=True,
        type=Path,
        help='Text file with one sample ID per line, in BCF column order.',
    )
    parser.add_argument(
        '--sex-tsv',
        required=True,
        type=Path,
        help='Two-column TSV: sample_id<TAB>sex_code (1=M, 2=F).',
    )
    parser.add_argument(
        '--output',
        required=True,
        type=Path,
        help='Output path for the bgzip-compressed TSV.',
    )
    args = parser.parse_args(argv)

    sample_order = load_sample_order(args.sample_order)
    sex_map = load_sex_map(args.sex_tsv)
    sample_indices = build_subset_indices(sample_order, sex_map)

    with args.output.open('wb') as out_fh:
        proc = subprocess.Popen(
            ['bgzip', '-c'],  # noqa: S607
            stdin=subprocess.PIPE,
            stdout=out_fh,
        )
        assert proc.stdin is not None
        writer = io.TextIOWrapper(proc.stdin, encoding='utf-8', write_through=True)
        try:
            run(sys.stdin, sample_indices, writer)
        finally:
            writer.flush()
            writer.close()
            proc.wait()
        if proc.returncode != 0:
            raise RuntimeError(f'bgzip exited with status {proc.returncode}')

    return 0


if __name__ == '__main__':
    sys.exit(main())
