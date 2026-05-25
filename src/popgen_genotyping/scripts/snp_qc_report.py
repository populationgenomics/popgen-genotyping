"""
Multi-cohort per-SNP QC synthesis.

Joins four sources on variant ID and emits a per-SNP exclusion list plus a
wide audit TSV:

- Illumina EGT INFO (`GenTrain_Score`, `Cluster_Sep`, training cluster counts)
  from the references-repo sample-less BCF, pre-extracted as a TSV.
- Merged-set PLINK2 site stats from `Plink2Qc`: variant missingness (`.vmiss`),
  allele frequency (`.afreq`), HWE (`.hardy`).
- Per-cohort PLINK2 site stats from `CohortPlinkQc`, reduced to call-rate
  range, allele-frequency χ² heterogeneity, and HWE-flip count across batches.
- Per-cohort observed THETA cluster moments from `CohortClusterStats`, pooled
  to per-cluster (AA/AB/BB) batch ICC. Informational — does not drive exclusion.

Final pass/fail is the boolean OR of every `fail_*` column produced; the
`flag_*` columns surface informational issues without driving exclusion.

Outputs:
- ``--output-audit-tsv``: bgzip-style ``.tsv.gz`` (one row per variant).
- ``--output-exclusion-list``: one variant ID per line for the failed set.
- ``--output-summary-tsv``: per-filter drop counts and pairwise overlaps.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

if TYPE_CHECKING:
    from collections.abc import Iterable


VARIANT_KEY: str = 'ID'
CLUSTERS: tuple[str, ...] = ('AA', 'AB', 'BB')

EGT_INFO_COLS: tuple[str, ...] = (
    'CHROM',
    'POS',
    'REF',
    'ALT',
    VARIANT_KEY,
    'GenTrain_Score',
    'Cluster_Sep',
    'N_AA_egt',
    'N_AB_egt',
    'N_BB_egt',
    'meanTHETA_AA_egt',
    'meanTHETA_AB_egt',
    'meanTHETA_BB_egt',
    'devTHETA_AA_egt',
    'devTHETA_AB_egt',
    'devTHETA_BB_egt',
)

# Boolean OR of these columns drives the final `exclude` flag.
FAIL_COLS: tuple[str, ...] = (
    'fail_gentrain',
    'fail_cluster_sep',
    'fail_vmiss',
    'fail_hwe',
    'fail_cross_batch_call_rate',
    'fail_cross_batch_af',
    'fail_egt_centroid_delta',
    'fail_cluster_icc',
    'fail_cluster_spread',
    'fail_cluster_separation',
)

FLAG_COLS: tuple[str, ...] = (
    'flag_empty_egt_cluster',
    'flag_monomorphic',
    'flag_hwe_flip',
)


def _strip_plink_hash(df: pd.DataFrame) -> pd.DataFrame:
    """Strip the leading ``#`` from a PLINK2 header (e.g. ``#CHROM`` → ``CHROM``)."""
    df.columns = df.columns.str.lstrip('#')
    return df


def load_egt_info(tsv_path: str | Path) -> pd.DataFrame:
    """
    Load the EGT INFO TSV emitted by an upstream ``bcftools query`` step.

    Args:
        tsv_path: Path to a header-less tab-separated file with columns matching
            ``EGT_INFO_COLS`` in declaration order.

    Returns:
        DataFrame keyed implicitly on ``VARIANT_KEY`` (no index set). Numeric
        columns are coerced; missing values appear as NaN.
    """
    df = pd.read_csv(
        tsv_path,
        sep='\t',
        header=None,
        names=list(EGT_INFO_COLS),
        na_values=['.', ''],
        dtype={
            'CHROM': 'string',
            'POS': 'Int64',
            'REF': 'string',
            'ALT': 'string',
            VARIANT_KEY: 'string',
        },
    )
    float_cols = (
        'GenTrain_Score',
        'Cluster_Sep',
        'meanTHETA_AA_egt',
        'meanTHETA_AB_egt',
        'meanTHETA_BB_egt',
        'devTHETA_AA_egt',
        'devTHETA_AB_egt',
        'devTHETA_BB_egt',
    )
    for col in float_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    for col in ('N_AA_egt', 'N_AB_egt', 'N_BB_egt'):
        df[col] = pd.to_numeric(df[col], errors='coerce').astype('Int64')
    return df


def load_merged_qc(vmiss_path: str | Path, afreq_path: str | Path, hardy_path: str | Path) -> pd.DataFrame:
    """
    Load and join the merged-set PLINK2 site-QC outputs.

    Args:
        vmiss_path: Path to ``.vmiss`` (variant missingness).
        afreq_path: Path to ``.afreq`` (allele frequency).
        hardy_path: Path to ``.hardy`` (Hardy-Weinberg autosomal).

    Returns:
        DataFrame with columns ``ID``, ``F_MISS``, ``ALT_FREQS``, ``OBS_CT``,
        and ``HWE_P``. HWE is autosomal-only — plink2's separate ``.hardy.x``
        / ``.hardy.y`` files are not consumed.
    """
    vmiss = _strip_plink_hash(pd.read_csv(vmiss_path, sep=r'\s+', engine='python'))
    afreq = _strip_plink_hash(pd.read_csv(afreq_path, sep=r'\s+', engine='python'))
    hardy = _strip_plink_hash(pd.read_csv(hardy_path, sep=r'\s+', engine='python'))

    vmiss = vmiss[[VARIANT_KEY, 'F_MISS']]
    afreq = afreq[[VARIANT_KEY, 'ALT_FREQS', 'OBS_CT']]
    hardy = hardy[[VARIANT_KEY, 'P']].rename(columns={'P': 'HWE_P'})

    return vmiss.merge(afreq, on=VARIANT_KEY, how='outer').merge(hardy, on=VARIANT_KEY, how='outer')


def load_cohort_qc(
    cohort_id: str,
    vmiss_path: str | Path,
    afreq_path: str | Path,
    hardy_path: str | Path,
) -> pd.DataFrame:
    """
    Load per-cohort PLINK2 site-QC outputs and tag rows with the cohort ID.

    Args:
        cohort_id: Stable cohort identifier (used as a row tag, not in the column
            namespace).
        vmiss_path: Path to ``.vmiss``.
        afreq_path: Path to ``.afreq``.
        hardy_path: Path to ``.hardy``.

    Returns:
        Long-format DataFrame with one row per (cohort, variant); columns are
        ``cohort_id``, ``ID``, ``F_MISS``, ``ALT_FREQS``, ``OBS_CT``, ``HWE_P``.
    """
    df = load_merged_qc(vmiss_path, afreq_path, hardy_path)
    df.insert(0, 'cohort_id', cohort_id)
    return df


def load_cluster_stats(cohort_id: str, stats_tsv_gz: str | Path) -> pd.DataFrame:
    """
    Load the per-cohort observed cluster-stats TSV produced by `CohortClusterStats`.

    Only ``sex_subset == 'all'`` rows are returned: autosomal variants emit a
    single ``all`` row, and the synthesis does not pool sex-stratified rows.

    Args:
        cohort_id: Stable cohort identifier.
        stats_tsv_gz: Path to the bgzipped TSV emitted by ``aggregate_cluster_stats.py``.

    Returns:
        Long-format DataFrame keyed by (cohort_id, ID); one row per variant.
        Carries ``n_<C>``, ``mean_THETA_<C>``, ``var_THETA_<C>`` for each
        ``C in CLUSTERS``.
    """
    df = pd.read_csv(stats_tsv_gz, sep='\t', compression='gzip', na_values=['NA'])
    df = df.rename(columns={'variant_id': VARIANT_KEY})
    df = df[df['sex_subset'] == 'all']

    keep_cols = [VARIANT_KEY]
    for cluster in CLUSTERS:
        keep_cols.extend([f'n_{cluster}', f'mean_THETA_{cluster}', f'var_THETA_{cluster}'])
    df = df[keep_cols].copy()
    df.insert(0, 'cohort_id', cohort_id)
    return df


def compute_cross_batch_metrics(cohort_qc: pd.DataFrame, hwe_p_min: float) -> pd.DataFrame:
    """
    Reduce per-cohort site QC to one row per variant of cross-batch metrics.

    Args:
        cohort_qc: Long-format frame from concatenating multiple
            `load_cohort_qc` results.
        hwe_p_min: HWE p-value cutoff used only to count per-batch flips
            (the actual HWE exclusion fires on the merged HWE, not this count).

    Returns:
        DataFrame with one row per ``ID`` and columns
        ``F_MISS_range``, ``AF_chi2_P``, ``HWE_flip_count``, ``n_batches``.
        ``AF_chi2_P`` is ``NaN`` for variants with <2 informative batches or
        all-identical allele counts (degenerate contingency table).
    """
    if cohort_qc.empty:
        return pd.DataFrame(columns=[VARIANT_KEY, 'F_MISS_range', 'AF_chi2_P', 'HWE_flip_count', 'n_batches'])

    grouped = cohort_qc.groupby(VARIANT_KEY, sort=False)

    f_miss_range = grouped['F_MISS'].agg(lambda s: s.max() - s.min() if s.notna().any() else np.nan)
    n_batches = grouped['cohort_id'].nunique()
    hwe_flip_count = grouped['HWE_P'].agg(lambda s: int((s < hwe_p_min).sum()))

    af_chi2_p: dict[str, float] = {}
    for variant_id, sub in grouped:
        af_chi2_p[variant_id] = _af_chi2_pvalue(sub['ALT_FREQS'].to_numpy(), sub['OBS_CT'].to_numpy())

    out = pd.DataFrame(
        {
            VARIANT_KEY: f_miss_range.index,
            'F_MISS_range': f_miss_range.to_numpy(),
            'AF_chi2_P': [af_chi2_p[v] for v in f_miss_range.index],
            'HWE_flip_count': hwe_flip_count.to_numpy(),
            'n_batches': n_batches.to_numpy(),
        },
    )
    return out.reset_index(drop=True)


def _af_chi2_pvalue(alt_freqs: np.ndarray, obs_ct: np.ndarray) -> float:
    """
    Compute χ² p-value for allele-count homogeneity across batches.

    Args:
        alt_freqs: Per-batch alternate-allele frequencies (proportion).
        obs_ct: Per-batch observed sample count.

    Returns:
        p-value from `scipy.stats.chi2_contingency`, or NaN if the table is
        degenerate (one batch, zero coverage, or zero column / row totals).
    """
    mask = np.isfinite(alt_freqs) & np.isfinite(obs_ct) & (obs_ct > 0)
    if mask.sum() < 2:  # noqa: PLR2004
        return float('nan')

    alt_ct = np.rint(alt_freqs[mask] * obs_ct[mask] * 2).astype(int)
    ref_ct = (2 * obs_ct[mask]).astype(int) - alt_ct

    # Both rows / both columns must be non-zero somewhere for chi2 to be defined.
    if alt_ct.sum() == 0 or ref_ct.sum() == 0:
        return float('nan')

    table = np.vstack([alt_ct, ref_ct])
    try:
        _, p, _, _ = chi2_contingency(table)
    except ValueError:
        return float('nan')
    return float(p)


def compute_cluster_icc(cluster_stats: pd.DataFrame) -> pd.DataFrame:
    """
    Pool per-cohort THETA cluster moments to per-variant ICC across batches.

    Within-between variance decomposition:
      pooled_mean = Σ n_i x̄_i / N
      s²_within   = Σ (n_i - 1) s²_i / (N - B)
      s²_between  = Σ n_i (x̄_i - pooled_mean)² / (N - 1)
      ICC         = s²_between / (s²_within + s²_between)

    Args:
        cluster_stats: Long-format frame concatenating `load_cluster_stats`
            outputs across cohorts.

    Returns:
        Wide DataFrame with one row per variant and columns
        ``pooled_mean_THETA_<C>``, ``ICC_THETA_<C>``, and
        ``pooled_var_within_THETA_<C>`` for each ``C in CLUSTERS``. NaN where
        insufficient observations exist.
    """
    pooled_keys = ('pooled_mean', 'ICC', 'pooled_var_within')
    if cluster_stats.empty:
        cols = [VARIANT_KEY] + [f'{prefix}_THETA_{c}' for prefix in pooled_keys for c in CLUSTERS]
        return pd.DataFrame(columns=cols)

    rows: list[dict[str, float | str]] = []
    for variant_id, sub in cluster_stats.groupby(VARIANT_KEY, sort=False):
        row: dict[str, float | str] = {VARIANT_KEY: variant_id}
        for cluster in CLUSTERS:
            mean, icc, var_within = _pool_cluster_moments(
                sub[f'n_{cluster}'].to_numpy(),
                sub[f'mean_THETA_{cluster}'].to_numpy(),
                sub[f'var_THETA_{cluster}'].to_numpy(),
            )
            row[f'pooled_mean_THETA_{cluster}'] = mean
            row[f'ICC_THETA_{cluster}'] = icc
            row[f'pooled_var_within_THETA_{cluster}'] = var_within
        rows.append(row)
    return pd.DataFrame(rows)


def _pool_cluster_moments(n: np.ndarray, mean: np.ndarray, var: np.ndarray) -> tuple[float, float, float]:
    """
    Within-between variance decomposition for a single (variant, cluster).

    Args:
        n: Per-batch sample counts in this cluster.
        mean: Per-batch THETA means (NaN allowed; pairs with n==0).
        var: Per-batch THETA sample variances.

    Returns:
        Tuple ``(pooled_mean, ICC, s2_within)``. Returns
        ``(NaN, NaN, NaN)`` if N == 0 or fewer than two batches contribute.
        ICC may still be NaN even when pooled_mean is defined if N - B <= 0
        (one observation per contributing batch).
    """
    valid = (n > 0) & np.isfinite(mean)
    if valid.sum() < 2:  # noqa: PLR2004
        return float('nan'), float('nan'), float('nan')
    n_v = n[valid].astype(float)
    m_v = mean[valid].astype(float)
    v_v = np.nan_to_num(var[valid].astype(float), nan=0.0)

    total = n_v.sum()
    if total < 2:  # noqa: PLR2004
        return float('nan'), float('nan'), float('nan')

    pooled_mean = float((n_v * m_v).sum() / total)
    n_batches = int(valid.sum())

    within_df = total - n_batches
    s2_within = float(((n_v - 1) * v_v).sum() / within_df) if within_df > 0 else float('nan')

    s2_between = float((n_v * (m_v - pooled_mean) ** 2).sum() / (total - 1))

    denom = (0.0 if np.isnan(s2_within) else s2_within) + s2_between
    icc = float(s2_between / denom) if denom > 0 else float('nan')
    return pooled_mean, icc, s2_within


def apply_filters(
    df: pd.DataFrame,
    *,
    gentrain_min: float,
    cluster_sep_min: float,
    vmiss_max: float,
    hwe_p_min: float,
    call_rate_range_max: float,
    af_chi2_p_min: float,
    egt_centroid_delta_max: float,
    cluster_icc_max: float,
    cluster_spread_max: float,
    cluster_separation_min: float,
) -> pd.DataFrame:
    """
    Add ``fail_*`` / ``flag_*`` columns to the per-variant audit frame.

    Args:
        df: Frame already wide-joined across EGT INFO / merged QC / cross-batch.
        gentrain_min: Lower bound for `GenTrain_Score`.
        cluster_sep_min: Lower bound for `Cluster_Sep`.
        vmiss_max: Upper bound for merged ``F_MISS``.
        hwe_p_min: Lower bound for merged ``HWE_P`` (autosomal).
        call_rate_range_max: Upper bound for max-min per-batch ``F_MISS`` spread.
        af_chi2_p_min: Lower bound for cross-batch allele-frequency χ² p-value.
        egt_centroid_delta_max: Upper bound for the largest |observed -
            EGT| pooled-THETA mean delta across AA/AB/BB clusters.
        cluster_icc_max: Upper bound for the largest cross-batch THETA ICC
            across AA/AB/BB clusters.
        cluster_spread_max: Upper bound for the largest pooled within-cluster
            THETA variance across AA/AB/BB clusters.
        cluster_separation_min: Lower bound for the smallest gap between
            adjacent pooled cluster means (AB-AA, BB-AB).

    Returns:
        Same frame with the columns named in ``FAIL_COLS`` and ``FLAG_COLS``
        added, plus the underlying ``egt_centroid_delta_max``,
        ``cluster_icc_max``, ``cluster_spread_max`` and ``cluster_separation_min``
        metric columns, ``exclude`` (bool), and ``exclusion_reasons`` (str).
        NaN metrics do not trigger their fail column.
    """
    out = df.copy()

    out['fail_gentrain'] = (out['GenTrain_Score'] < gentrain_min).fillna(False)
    out['fail_cluster_sep'] = (out['Cluster_Sep'] < cluster_sep_min).fillna(False)
    out['flag_empty_egt_cluster'] = ((out['N_AA_egt'] == 0) | (out['N_AB_egt'] == 0) | (out['N_BB_egt'] == 0)).fillna(
        False,
    )

    out['fail_vmiss'] = (out['F_MISS'] > vmiss_max).fillna(False)
    out['fail_hwe'] = (out['HWE_P'] < hwe_p_min).fillna(False)
    out['flag_monomorphic'] = ((out['ALT_FREQS'] == 0) | (out['ALT_FREQS'] == 1)).fillna(False)

    out['fail_cross_batch_call_rate'] = (out['F_MISS_range'] > call_rate_range_max).fillna(False)
    out['fail_cross_batch_af'] = (out['AF_chi2_P'] < af_chi2_p_min).fillna(False)
    out['flag_hwe_flip'] = (out.get('HWE_flip_count', pd.Series(0, index=out.index)) > 0).fillna(False)

    centroid_deltas = pd.concat(
        [(out[f'pooled_mean_THETA_{c}'] - out[f'meanTHETA_{c}_egt']).abs() for c in CLUSTERS],
        axis=1,
    )
    out['egt_centroid_delta_max'] = centroid_deltas.max(axis=1)
    out['fail_egt_centroid_delta'] = (out['egt_centroid_delta_max'] > egt_centroid_delta_max).fillna(False)

    out['cluster_icc_max'] = out[[f'ICC_THETA_{c}' for c in CLUSTERS]].max(axis=1)
    out['fail_cluster_icc'] = (out['cluster_icc_max'] > cluster_icc_max).fillna(False)

    out['cluster_spread_max'] = out[[f'pooled_var_within_THETA_{c}' for c in CLUSTERS]].max(axis=1)
    out['fail_cluster_spread'] = (out['cluster_spread_max'] > cluster_spread_max).fillna(False)

    gap_aa_ab = out['pooled_mean_THETA_AB'] - out['pooled_mean_THETA_AA']
    gap_ab_bb = out['pooled_mean_THETA_BB'] - out['pooled_mean_THETA_AB']
    out['cluster_separation_min'] = pd.concat([gap_aa_ab, gap_ab_bb], axis=1).min(axis=1)
    out['fail_cluster_separation'] = (out['cluster_separation_min'] < cluster_separation_min).fillna(False)

    fail_matrix = out[list(FAIL_COLS)].fillna(False).astype(bool)
    out['exclude'] = fail_matrix.any(axis=1)
    out['exclusion_reasons'] = fail_matrix.apply(_join_reasons, axis=1)
    return out


def _join_reasons(row: pd.Series) -> str:
    """Comma-join the names of all truthy fail columns for one variant."""
    return ','.join(name for name, val in row.items() if bool(val))


def build_audit_frame(
    egt: pd.DataFrame,
    merged: pd.DataFrame,
    cross_batch: pd.DataFrame,
    cluster_icc: pd.DataFrame,
) -> pd.DataFrame:
    """
    Wide-join the four per-variant frames into a single audit table.

    Outer-joins on ``ID``. EGT is the natural anchor (one row per Illumina
    probe in the manifest) but we don't drop variants present in PLINK QC
    but missing from EGT — that's a signal worth surfacing rather than
    silently dropping.
    """
    return (
        egt.merge(merged, on=VARIANT_KEY, how='outer')
        .merge(cross_batch, on=VARIANT_KEY, how='left')
        .merge(cluster_icc, on=VARIANT_KEY, how='left')
    )


def write_outputs(
    audit_df: pd.DataFrame,
    audit_path: Path,
    exclusion_path: Path,
    summary_path: Path,
) -> None:
    """Persist the audit TSV, exclusion snplist, and summary stats TSV."""
    audit_df.to_csv(audit_path, sep='\t', index=False, compression='gzip')

    excluded = audit_df.loc[audit_df['exclude'], VARIANT_KEY].dropna().astype(str)
    with exclusion_path.open('w') as fh:
        for variant_id in excluded:
            fh.write(f'{variant_id}\n')

    summary = _summarise(audit_df)
    summary.to_csv(summary_path, sep='\t', index=False)


def _summarise(audit_df: pd.DataFrame) -> pd.DataFrame:
    """
    Per-filter drop counts.

    Returns a long-format frame with columns ``filter`` and ``n_variants``.
    Includes the headline ``exclude_any`` row plus a row per ``fail_*`` and
    ``flag_*`` column.
    """
    rows: list[dict[str, object]] = [{'filter': 'total', 'n_variants': len(audit_df)}]
    rows.append({'filter': 'exclude_any', 'n_variants': int(audit_df['exclude'].sum())})
    for col in (*FAIL_COLS, *FLAG_COLS):
        if col in audit_df.columns:
            rows.append({'filter': col, 'n_variants': int(audit_df[col].sum())})
    return pd.DataFrame(rows)


def _parse_triples(items: Iterable[str]) -> list[tuple[str, str]]:
    """Split each ``cohort_id:path`` argument into a tuple."""
    out: list[tuple[str, str]] = []
    for raw in items:
        cohort_id, _, path = raw.partition(':')
        if not cohort_id or not path:
            raise ValueError(f'expected cohort_id:path, got {raw!r}')
        out.append((cohort_id, path))
    return out


def main(argv: list[str] | None = None) -> int:
    """CLI entry point. Returns a process exit code."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--egt-info-tsv', required=True, type=Path)
    parser.add_argument('--merged-vmiss', required=True, type=Path)
    parser.add_argument('--merged-afreq', required=True, type=Path)
    parser.add_argument('--merged-hardy', required=True, type=Path)
    parser.add_argument(
        '--cohort-vmiss',
        action='append',
        default=[],
        help='Per-cohort vmiss as cohort_id:path. Repeat per cohort.',
    )
    parser.add_argument('--cohort-afreq', action='append', default=[])
    parser.add_argument('--cohort-hardy', action='append', default=[])
    parser.add_argument(
        '--cohort-cluster-stats',
        action='append',
        default=[],
        help='Per-cohort cluster_stats.tsv.gz as cohort_id:path. Repeat per cohort.',
    )
    parser.add_argument('--gentrain-min', type=float, required=True)
    parser.add_argument('--cluster-sep-min', type=float, required=True)
    parser.add_argument('--vmiss-max', type=float, required=True)
    parser.add_argument('--hwe-p-min', type=float, required=True)
    parser.add_argument('--call-rate-range-max', type=float, required=True)
    parser.add_argument('--af-chi2-p-min', type=float, required=True)
    parser.add_argument('--egt-centroid-delta-max', type=float, required=True)
    parser.add_argument('--cluster-icc-max', type=float, required=True)
    parser.add_argument('--cluster-spread-max', type=float, required=True)
    parser.add_argument('--cluster-separation-min', type=float, required=True)
    parser.add_argument('--output-audit-tsv', required=True, type=Path)
    parser.add_argument('--output-exclusion-list', required=True, type=Path)
    parser.add_argument('--output-summary-tsv', required=True, type=Path)
    args = parser.parse_args(argv)

    vmiss_triples = dict(_parse_triples(args.cohort_vmiss))
    afreq_triples = dict(_parse_triples(args.cohort_afreq))
    hardy_triples = dict(_parse_triples(args.cohort_hardy))
    cluster_stats_triples = dict(_parse_triples(args.cohort_cluster_stats))

    cohort_ids = sorted(set(vmiss_triples) | set(afreq_triples) | set(hardy_triples))
    if not (set(vmiss_triples) == set(afreq_triples) == set(hardy_triples)):
        print(
            f'warning: per-cohort plink QC inputs disagree on cohort set: '
            f'vmiss={sorted(vmiss_triples)} afreq={sorted(afreq_triples)} hardy={sorted(hardy_triples)}',
            file=sys.stderr,
        )

    cohort_qc_frames: list[pd.DataFrame] = []
    for cohort_id in cohort_ids:
        if cohort_id in vmiss_triples and cohort_id in afreq_triples and cohort_id in hardy_triples:
            cohort_qc_frames.append(
                load_cohort_qc(
                    cohort_id,
                    vmiss_triples[cohort_id],
                    afreq_triples[cohort_id],
                    hardy_triples[cohort_id],
                ),
            )
    cohort_qc = (
        pd.concat(cohort_qc_frames, ignore_index=True)
        if cohort_qc_frames
        else pd.DataFrame(columns=['cohort_id', VARIANT_KEY, 'F_MISS', 'ALT_FREQS', 'OBS_CT', 'HWE_P'])
    )

    cluster_stats_frames: list[pd.DataFrame] = [
        load_cluster_stats(cohort_id, path) for cohort_id, path in cluster_stats_triples.items()
    ]
    cluster_stats = pd.concat(cluster_stats_frames, ignore_index=True) if cluster_stats_frames else pd.DataFrame()

    egt = load_egt_info(args.egt_info_tsv)
    merged = load_merged_qc(args.merged_vmiss, args.merged_afreq, args.merged_hardy)
    cross_batch = compute_cross_batch_metrics(cohort_qc, hwe_p_min=args.hwe_p_min)
    cluster_icc = compute_cluster_icc(cluster_stats)

    audit = build_audit_frame(egt, merged, cross_batch, cluster_icc)
    audit = apply_filters(
        audit,
        gentrain_min=args.gentrain_min,
        cluster_sep_min=args.cluster_sep_min,
        vmiss_max=args.vmiss_max,
        hwe_p_min=args.hwe_p_min,
        call_rate_range_max=args.call_rate_range_max,
        af_chi2_p_min=args.af_chi2_p_min,
        egt_centroid_delta_max=args.egt_centroid_delta_max,
        cluster_icc_max=args.cluster_icc_max,
        cluster_spread_max=args.cluster_spread_max,
        cluster_separation_min=args.cluster_separation_min,
    )
    write_outputs(audit, args.output_audit_tsv, args.output_exclusion_list, args.output_summary_tsv)
    return 0


if __name__ == '__main__':
    sys.exit(main())
