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
    'pass_gentrain',
    'pass_cluster_sep',
    'pass_fmiss',
    'pass_strand',
    'pass',
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
    """Append per-filter and aggregate pass columns to ``df``.

    NaN inputs (missing EGT scores or missing F_MISS) are treated as failures
    so the variant lands on the exclusion list rather than silently passing.

    Args:
        df: Joined EGT + vmiss frame, already annotated with ``strand_ambiguous``.
        gentrain_min: Inclusive lower bound for ``GenTrain_Score``.
        cluster_sep_min: Inclusive lower bound for ``Cluster_Sep``.
        fmiss_max: Inclusive upper bound for merged-set ``F_MISS``.
        exclude_strand_ambiguous: If True, ``strand_ambiguous`` SNPs fail
            ``pass_strand``; otherwise the check is a no-op.

    Returns:
        Copy of ``df`` with five boolean columns: ``pass_gentrain``,
        ``pass_cluster_sep``, ``pass_fmiss``, ``pass_strand``, ``pass``.
    """
    out: pd.DataFrame = df.copy()
    out['pass_gentrain'] = out['GenTrain_Score'].fillna(-1.0) >= gentrain_min
    out['pass_cluster_sep'] = out['Cluster_Sep'].fillna(-1.0) >= cluster_sep_min
    out['pass_fmiss'] = out['F_MISS'].fillna(1.0) <= fmiss_max
    out['pass_strand'] = ~out['strand_ambiguous'] if exclude_strand_ambiguous else True
    out['pass'] = out['pass_gentrain'] & out['pass_cluster_sep'] & out['pass_fmiss'] & out['pass_strand']
    return out


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    """Build a per-filter drop-count summary.

    Args:
        df: Annotated frame from :func:`apply_filters`.

    Returns:
        Two-column DataFrame (``metric``, ``value``) suitable for direct TSV write.
    """
    rows: list[tuple[str, int]] = [
        ('total_variants', len(df)),
        ('fail_gentrain', int((~df['pass_gentrain']).sum())),
        ('fail_cluster_sep', int((~df['pass_cluster_sep']).sum())),
        ('fail_fmiss', int((~df['pass_fmiss']).sum())),
        ('strand_ambiguous', int(df['strand_ambiguous'].sum())),
        ('fail_strand_filter', int((~df['pass_strand']).sum())),
        ('total_excluded', int((~df['pass']).sum())),
        ('total_retained', int(df['pass'].sum())),
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

    excluded_ids: pd.Series = df.loc[~df['pass'], 'ID']
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
