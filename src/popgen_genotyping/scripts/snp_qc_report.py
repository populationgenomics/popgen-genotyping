"""Per-SNP QC filter for the genotyping pipeline.

Joins the Illumina EGT INFO TSV (extracted from the references-repo
sample-less BCF) with the merged-set ``plink2 --missing`` per-variant output,
applies four threshold filters (``GenTrain_Score``, ``Cluster_Sep``,
``F_MISS``, and an optional ``{A,T}``/``{C,G}`` strand-ambiguity drop), and
emits an exclusion ``.snplist`` plus an audit TSV and per-filter summary.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd

EGT_TSV_COLUMNS: list[str] = [
    'CHROM',
    'POS',
    'REF',
    'ALT',
    'ID',
    'GenTrain_Score',
    'Cluster_Sep',
]

AUDIT_COLUMNS: list[str] = [
    'CHROM',
    'POS',
    'REF',
    'ALT',
    'ID',
    'GenTrain_Score',
    'Cluster_Sep',
    'F_MISS',
    'strand_ambiguous',
    'fail_gentrain',
    'fail_cluster_sep',
    'fail_fmiss',
    'fail_strand',
    'fail',
]


def load_egt_info(path: Path) -> pd.DataFrame:
    """Read a header-less EGT INFO TSV produced by ``bcftools query``.

    Args:
        path: TSV with columns CHROM, POS, REF, ALT, ID, GenTrain_Score, Cluster_Sep.

    Returns:
        DataFrame with numeric coercion applied to the score columns
        (``.`` and empty string become NaN).
    """
    df: pd.DataFrame = pd.read_csv(
        path,
        sep='\t',
        header=None,
        names=EGT_TSV_COLUMNS,
        na_values=['.', ''],
        dtype={'CHROM': str, 'ID': str, 'REF': str, 'ALT': str},
    )
    df['GenTrain_Score'] = pd.to_numeric(df['GenTrain_Score'], errors='coerce')
    df['Cluster_Sep'] = pd.to_numeric(df['Cluster_Sep'], errors='coerce')
    return df


def load_vmiss(path: Path) -> pd.DataFrame:
    """Read a ``plink2 --missing`` per-variant TSV.

    Args:
        path: Path to a ``.vmiss`` file (header line starts with ``#CHROM``).

    Returns:
        DataFrame with the ``#`` stripped from the header and a single
        ``ID``/``F_MISS`` projection.
    """
    df: pd.DataFrame = pd.read_csv(path, sep='\t')
    df.columns = df.columns.str.lstrip('#')
    return df[['ID', 'F_MISS']].copy()


def add_strand_ambiguous_flag(df: pd.DataFrame) -> pd.DataFrame:
    """Annotate ``df`` with a boolean ``strand_ambiguous`` column.

    A SNP is strand-ambiguous iff its ``{REF, ALT}`` is ``{A, T}`` or ``{C, G}``.
    Indels and missing alleles are never flagged.

    Args:
        df: Frame with ``REF`` and ``ALT`` string columns.

    Returns:
        Copy of ``df`` with a ``strand_ambiguous`` boolean column appended.
    """
    ref: pd.Series = df['REF'].astype('string').str.upper()
    alt: pd.Series = df['ALT'].astype('string').str.upper()
    snv: pd.Series = (ref.str.len() == 1) & (alt.str.len() == 1)
    at: pd.Series = ((ref == 'A') & (alt == 'T')) | ((ref == 'T') & (alt == 'A'))
    cg: pd.Series = ((ref == 'C') & (alt == 'G')) | ((ref == 'G') & (alt == 'C'))
    out: pd.DataFrame = df.copy()
    out['strand_ambiguous'] = (snv & (at | cg)).fillna(False).astype(bool)
    return out


def apply_filters(
    df: pd.DataFrame,
    *,
    gentrain_min: float,
    cluster_sep_min: float,
    fmiss_max: float,
    exclude_strand_ambiguous: bool,
) -> pd.DataFrame:
    """Append per-filter and aggregate fail columns to ``df``.

    NaN inputs (missing EGT scores or missing F_MISS) are treated as failures
    so the variant lands on the exclusion list rather than silently passing.

    Args:
        df: Joined EGT + vmiss frame, already annotated with ``strand_ambiguous``.
        gentrain_min: Inclusive lower bound for ``GenTrain_Score``.
        cluster_sep_min: Inclusive lower bound for ``Cluster_Sep``.
        fmiss_max: Inclusive upper bound for merged-set ``F_MISS``.
        exclude_strand_ambiguous: If True, ``strand_ambiguous`` SNPs set
            ``fail_strand``; otherwise the check is a no-op (always False).

    Returns:
        Copy of ``df`` with five boolean columns: ``fail_gentrain``,
        ``fail_cluster_sep``, ``fail_fmiss``, ``fail_strand``, ``fail``.
        ``True`` indicates the variant failed that filter (or any filter, for
        the aggregate ``fail`` column).
    """
    out: pd.DataFrame = df.copy()
    out['fail_gentrain'] = out['GenTrain_Score'].fillna(-1.0) < gentrain_min
    out['fail_cluster_sep'] = out['Cluster_Sep'].fillna(-1.0) < cluster_sep_min
    out['fail_fmiss'] = out['F_MISS'].fillna(1.0) > fmiss_max
    out['fail_strand'] = out['strand_ambiguous'] if exclude_strand_ambiguous else False
    out['fail'] = out['fail_gentrain'] | out['fail_cluster_sep'] | out['fail_fmiss'] | out['fail_strand']
    return out


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    """Build a per-filter drop-count summary.

    The output emits two complementary groups of per-filter rows:

    * ``first_fail_<filter>`` — count of variants for which ``<filter>`` was
      the first failing check in the priority cascade
      ``gentrain → cluster_sep → fmiss → strand``. Disjoint across filters:
      each excluded variant contributes to exactly one row, so the four
      ``first_fail_*`` rows sum to ``total_excluded``.
    * ``fail_<filter>`` — count of variants that fail ``<filter>`` regardless
      of any other filter. Overlapping across filters; useful when tuning a
      single threshold in isolation.

    Args:
        df: Annotated frame from :func:`apply_filters`.

    Returns:
        Two-column DataFrame (``metric``, ``value``) suitable for direct TSV write.
    """
    first_fail_gentrain: pd.Series = df['fail_gentrain']
    first_fail_cluster_sep: pd.Series = ~df['fail_gentrain'] & df['fail_cluster_sep']
    first_fail_fmiss: pd.Series = ~df['fail_gentrain'] & ~df['fail_cluster_sep'] & df['fail_fmiss']
    first_fail_strand: pd.Series = (
        ~df['fail_gentrain'] & ~df['fail_cluster_sep'] & ~df['fail_fmiss'] & df['fail_strand']
    )
    rows: list[tuple[str, int]] = [
        ('total_variants', len(df)),
        ('first_fail_gentrain', int(first_fail_gentrain.sum())),
        ('first_fail_cluster_sep', int(first_fail_cluster_sep.sum())),
        ('first_fail_fmiss', int(first_fail_fmiss.sum())),
        ('first_fail_strand', int(first_fail_strand.sum())),
        ('fail_gentrain', int(df['fail_gentrain'].sum())),
        ('fail_cluster_sep', int(df['fail_cluster_sep'].sum())),
        ('fail_fmiss', int(df['fail_fmiss'].sum())),
        ('strand_ambiguous', int(df['strand_ambiguous'].sum())),
        ('fail_strand_filter', int(df['fail_strand'].sum())),
        ('total_excluded', int(df['fail'].sum())),
        ('total_retained', int((~df['fail']).sum())),
    ]
    return pd.DataFrame(rows, columns=['metric', 'value'])


def write_outputs(
    df: pd.DataFrame,
    *,
    audit_path: Path,
    exclusion_path: Path,
    summary_path: Path,
) -> None:
    """Write the audit TSV (gzip), exclusion ``.snplist``, and summary TSV.

    Args:
        df: Annotated frame from :func:`apply_filters`.
        audit_path: Destination for the bgzippable per-variant audit TSV (``.gz``).
        exclusion_path: Destination for the newline-delimited failed-ID list.
        summary_path: Destination for the per-filter summary TSV.
    """
    with gzip.open(audit_path, 'wt') as audit_handle:
        df[AUDIT_COLUMNS].to_csv(audit_handle, sep='\t', index=False)

    excluded_ids: pd.Series = df.loc[df['fail'], 'ID']
    exclusion_path.write_text('\n'.join(excluded_ids.astype(str)) + ('\n' if len(excluded_ids) else ''))

    summarise(df).to_csv(summary_path, sep='\t', index=False)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse the CLI arguments.

    Args:
        argv: Optional override of ``sys.argv[1:]`` for testing.

    Returns:
        Parsed namespace.
    """
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--egt-info-tsv', type=Path, required=True)
    p.add_argument('--merged-vmiss', type=Path, required=True)
    p.add_argument('--gentrain-min', type=float, required=True)
    p.add_argument('--cluster-sep-min', type=float, required=True)
    p.add_argument('--fmiss-max', type=float, required=True)
    p.add_argument(
        '--exclude-strand-ambiguous',
        action='store_true',
        help='Also exclude {A,T}/{C,G} SNPs that cannot be strand-resolved.',
    )
    p.add_argument('--output-audit-tsv', type=Path, required=True)
    p.add_argument('--output-exclusion-list', type=Path, required=True)
    p.add_argument('--output-summary-tsv', type=Path, required=True)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Load the inputs, apply filters, and write the outputs.

    Args:
        argv: Optional override of ``sys.argv[1:]`` for testing.

    Returns:
        Process exit code (0 on success).
    """
    args: argparse.Namespace = parse_args(argv)
    df_egt: pd.DataFrame = load_egt_info(args.egt_info_tsv)
    df_vmiss: pd.DataFrame = load_vmiss(args.merged_vmiss)
    df: pd.DataFrame = df_egt.merge(df_vmiss, on='ID', how='left')
    df = add_strand_ambiguous_flag(df)
    df = apply_filters(
        df,
        gentrain_min=args.gentrain_min,
        cluster_sep_min=args.cluster_sep_min,
        fmiss_max=args.fmiss_max,
        exclude_strand_ambiguous=args.exclude_strand_ambiguous,
    )
    write_outputs(
        df,
        audit_path=args.output_audit_tsv,
        exclusion_path=args.output_exclusion_list,
        summary_path=args.output_summary_tsv,
    )
    return 0


if __name__ == '__main__':
    sys.exit(main())
