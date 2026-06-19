"""Per-SNP QC filter for the genotyping pipeline.

Joins the Illumina EGT INFO TSV (extracted from the references-repo
sample-less BCF) with the merged-set ``plink2 --missing`` per-variant output
and a ``plink2 --hwe`` pass snplist, applies the per-filter thresholds
(``GenTrain_Score``, ``Cluster_Sep``, ``F_MISS``, Hardy-Weinberg), and emits
an inclusion ``.snplist`` of passing variants (for ``plink2 --extract``) plus
an audit TSV and per-filter summary.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)


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
    'pass_gentrain',
    'pass_cluster_sep',
    'pass_fmiss',
    'pass_hwe',
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
    expected_id = (
        df['CHROM'].astype(str)
        + ':'
        + df['POS'].astype(str)
        + ':'
        + df['REF'].astype(str)
        + ':'
        + df['ALT'].astype(str)
    )

    bad_id = df['ID'].ne(expected_id)

    if bad_id.any():
        logging.info(f'Reformatted {bad_id.sum()} EGT IDs that did not match CHROM:POS:REF:ALT format.')
        df.loc[bad_id, 'ID'] = expected_id[bad_id]
    df['GenTrain_Score'] = pd.to_numeric(df['GenTrain_Score'], errors='coerce')
    df['Cluster_Sep'] = pd.to_numeric(df['Cluster_Sep'], errors='coerce')
    return df


def load_vmiss(path: Path) -> pd.DataFrame:
    """Read a ``plink2 --missing`` per-variant TSV.

    Args:
        path: Path to a ``.vmiss`` file (header line starts with ``#CHROM``).

    Returns:
        DataFrame with just the ``ID`` and ``F_MISS`` columns. ``ID`` carries no
        ``#`` prefix in a ``.vmiss`` header (only ``#CHROM`` does), so the two
        columns are read directly.
    """
    df: pd.DataFrame = pd.read_csv(path, sep='\t', usecols=['ID', 'F_MISS'])
    return df


def load_hwe_pass_ids(path: Path) -> set[str]:
    """Read a ``plink2 --write-snplist`` output as a set of passing variant IDs.

    Args:
        path: Newline-delimited ``.snplist`` (the post-filter variant set after
            ``plink2 --hwe``).

    Returns:
        Set of variant IDs that passed the HWE filter. Blank lines are skipped.
    """
    text: str = path.read_text()
    return {line.strip() for line in text.splitlines() if line.strip()}


def apply_filters(
    df: pd.DataFrame,
    *,
    gentrain_min: float,
    cluster_sep_min: float,
    fmiss_max: float,
    hwe_pass_ids: set[str],
) -> pd.DataFrame:
    """Append per-filter and aggregate pass columns to ``df``.

    A variant passes a filter when it meets the threshold; the aggregate
    ``pass`` column is the conjunction of all four. NaN inputs (missing EGT
    scores or missing F_MISS) make the threshold comparison ``False``, so a
    variant with missing data does not pass and is kept off the inclusion list.
    HWE pass is membership of the ``plink2 --hwe`` pass list.

    Args:
        df: Joined EGT + vmiss frame.
        gentrain_min: Inclusive lower bound for ``GenTrain_Score``.
        cluster_sep_min: Inclusive lower bound for ``Cluster_Sep``.
        fmiss_max: Inclusive upper bound for merged-set ``F_MISS``.
        hwe_pass_ids: Variant IDs that passed ``plink2 --hwe``.

    Returns:
        Copy of ``df`` with five boolean columns: ``pass_gentrain``,
        ``pass_cluster_sep``, ``pass_fmiss``, ``pass_hwe``, ``pass``.
        ``True`` indicates the variant passed that filter (or every filter, for
        the aggregate ``pass`` column).
    """
    out: pd.DataFrame = df.copy()
    out['pass_gentrain'] = out['GenTrain_Score'] >= gentrain_min
    out['pass_cluster_sep'] = out['Cluster_Sep'] >= cluster_sep_min
    out['pass_fmiss'] = out['F_MISS'] <= fmiss_max
    out['pass_hwe'] = out['ID'].astype(str).isin(hwe_pass_ids)
    out['pass'] = out['pass_gentrain'] & out['pass_cluster_sep'] & out['pass_fmiss'] & out['pass_hwe']
    return out


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    """Build a per-filter drop-count summary.

    The output emits two complementary groups of per-filter rows:

    * ``first_fail_<filter>`` — count of variants for which ``<filter>`` was
      the first failing check in the priority cascade
      ``gentrain → cluster_sep → fmiss → hwe``. Disjoint across filters: each
      excluded variant contributes to exactly one row, so the four
      ``first_fail_*`` rows sum to ``total_excluded``.
    * ``fail_<filter>`` — count of variants that fail ``<filter>`` regardless
      of any other filter. Overlapping across filters; useful when tuning a
      single threshold in isolation.

    Args:
        df: Annotated frame from :func:`apply_filters`.

    Returns:
        Two-column DataFrame (``metric``, ``value``) suitable for direct TSV write.
    """
    fail_gentrain: pd.Series = ~df['pass_gentrain']
    fail_cluster_sep: pd.Series = ~df['pass_cluster_sep']
    fail_fmiss: pd.Series = ~df['pass_fmiss']
    fail_hwe: pd.Series = ~df['pass_hwe']
    first_fail_gentrain: pd.Series = fail_gentrain
    first_fail_cluster_sep: pd.Series = df['pass_gentrain'] & fail_cluster_sep
    first_fail_fmiss: pd.Series = df['pass_gentrain'] & df['pass_cluster_sep'] & fail_fmiss
    first_fail_hwe: pd.Series = df['pass_gentrain'] & df['pass_cluster_sep'] & df['pass_fmiss'] & fail_hwe
    rows: list[tuple[str, int]] = [
        ('total_variants', len(df)),
        ('first_fail_gentrain', int(first_fail_gentrain.sum())),
        ('first_fail_cluster_sep', int(first_fail_cluster_sep.sum())),
        ('first_fail_fmiss', int(first_fail_fmiss.sum())),
        ('first_fail_hwe', int(first_fail_hwe.sum())),
        ('fail_gentrain', int(fail_gentrain.sum())),
        ('fail_cluster_sep', int(fail_cluster_sep.sum())),
        ('fail_fmiss', int(fail_fmiss.sum())),
        ('fail_hwe', int(fail_hwe.sum())),
        ('total_excluded', int((~df['pass']).sum())),
        ('total_retained', int(df['pass'].sum())),
    ]
    return pd.DataFrame(rows, columns=['metric', 'value'])


def write_outputs(
    df: pd.DataFrame,
    *,
    audit_path: Path,
    inclusion_path: Path,
    summary_path: Path,
) -> None:
    """Write the audit TSV (gzip), inclusion ``.snplist``, and summary TSV.

    Args:
        df: Annotated frame from :func:`apply_filters`.
        audit_path: Destination for the bgzippable per-variant audit TSV (``.gz``).
        inclusion_path: Destination for the newline-delimited passing-ID list
            (the variants to keep, for ``plink2 --extract``).
        summary_path: Destination for the per-filter summary TSV.
    """
    with gzip.open(audit_path, 'wt') as audit_handle:
        df[AUDIT_COLUMNS].to_csv(audit_handle, sep='\t', index=False)

    included_ids: pd.Series = df.loc[df['pass'], 'ID']
    inclusion_path.write_text('\n'.join(included_ids.astype(str)) + ('\n' if len(included_ids) else ''))

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
    p.add_argument('--hwe-pass-snplist', type=Path, required=True)
    p.add_argument('--gentrain-min', type=float, required=True)
    p.add_argument('--cluster-sep-min', type=float, required=True)
    p.add_argument('--fmiss-max', type=float, required=True)
    p.add_argument('--output-audit-tsv', type=Path, required=True)
    p.add_argument('--output-inclusion-list', type=Path, required=True)
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
    hwe_pass_ids: set[str] = load_hwe_pass_ids(args.hwe_pass_snplist)
    df: pd.DataFrame = df_egt.merge(df_vmiss, on='ID', how='left')
    df = apply_filters(
        df,
        gentrain_min=args.gentrain_min,
        cluster_sep_min=args.cluster_sep_min,
        fmiss_max=args.fmiss_max,
        hwe_pass_ids=hwe_pass_ids,
    )
    write_outputs(
        df,
        audit_path=args.output_audit_tsv,
        inclusion_path=args.output_inclusion_list,
        summary_path=args.output_summary_tsv,
    )
    return 0


if __name__ == '__main__':
    sys.exit(main())
