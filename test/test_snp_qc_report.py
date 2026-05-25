"""Tests for scripts/snp_qc_report.py — pure DataFrame logic, no file IO."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from popgen_genotyping.scripts.snp_qc_report import (
    FAIL_COLS,
    apply_filters,
    build_audit_frame,
    compute_cluster_icc,
    compute_cross_batch_metrics,
    _af_chi2_pvalue,
    _pool_cluster_moments,
)

DEFAULT_THRESHOLDS: dict[str, float] = {
    'gentrain_min': 0.7,
    'cluster_sep_min': 0.4,
    'vmiss_max': 0.02,
    'hwe_p_min': 1e-6,
    'call_rate_range_max': 0.05,
    'af_chi2_p_min': 1e-6,
    'egt_centroid_delta_max': 0.1,
    'cluster_icc_max': 0.5,
    'cluster_spread_max': 0.01,
    'cluster_separation_min': 0.15,
}


# -- _af_chi2_pvalue ---------------------------------------------------------


class TestAfChi2:
    """The contingency-test helper covers most cross-batch AF edge cases."""

    def test_single_batch_returns_nan(self) -> None:
        """One batch is not a 2x≥2 table."""
        p = _af_chi2_pvalue(np.array([0.5]), np.array([100]))
        assert math.isnan(p)

    def test_identical_frequencies_pvalue_near_one(self) -> None:
        """Same AF across all batches -> p ≈ 1."""
        p = _af_chi2_pvalue(np.array([0.5, 0.5, 0.5]), np.array([100, 100, 100]))
        assert p > 0.99

    def test_divergent_frequencies_small_pvalue(self) -> None:
        """0.0 vs 1.0 with large N -> p ≈ 0."""
        p = _af_chi2_pvalue(np.array([0.0, 1.0]), np.array([500, 500]))
        assert p < 1e-100

    def test_zero_obs_count_skipped(self) -> None:
        """Batches with OBS_CT == 0 drop out of the test."""
        p = _af_chi2_pvalue(np.array([0.5, 0.5, 0.5]), np.array([100, 0, 100]))
        assert p > 0.99

    def test_monomorphic_returns_nan(self) -> None:
        """All-zero ALT counts -> degenerate table -> NaN."""
        p = _af_chi2_pvalue(np.array([0.0, 0.0]), np.array([100, 100]))
        assert math.isnan(p)


# -- _pool_cluster_moments ---------------------------------------------------


class TestPoolClusterMoments:
    """Within-between variance decomposition for a single (variant, cluster)."""

    def test_single_batch_returns_all_nan(self) -> None:
        """Need ≥2 batches for ICC to be defined."""
        m, icc, var_within = _pool_cluster_moments(np.array([10.0]), np.array([0.5]), np.array([0.01]))
        assert math.isnan(m) and math.isnan(icc) and math.isnan(var_within)

    def test_two_batches_with_identical_means_zero_icc(self) -> None:
        """No between-batch variation -> ICC = 0; within-cluster variance preserved."""
        m, icc, var_within = _pool_cluster_moments(
            np.array([10.0, 10.0]),
            np.array([0.5, 0.5]),
            np.array([0.01, 0.01]),
        )
        assert m == 0.5
        assert icc == 0.0
        assert var_within == pytest.approx(0.01)

    def test_two_batches_with_divergent_means_high_icc(self) -> None:
        """Strong between-batch shift swamps tiny within-batch variance."""
        m, icc, _ = _pool_cluster_moments(
            np.array([10.0, 10.0]),
            np.array([0.1, 0.9]),
            np.array([0.001, 0.001]),
        )
        assert m == 0.5
        assert icc > 0.95

    def test_zero_n_batch_dropped(self) -> None:
        """A batch with n==0 must not pollute the pool."""
        m, icc, _ = _pool_cluster_moments(
            np.array([10.0, 0.0, 10.0]),
            np.array([0.5, 0.5, 0.5]),
            np.array([0.01, np.nan, 0.01]),
        )
        assert m == 0.5
        assert icc == 0.0

    def test_pooled_mean_is_n_weighted(self) -> None:
        """Bigger batch dominates the pooled mean."""
        m, _, _ = _pool_cluster_moments(
            np.array([90.0, 10.0]),
            np.array([0.5, 1.5]),
            np.array([0.0, 0.0]),
        )
        # weighted average: (90*0.5 + 10*1.5) / 100 = 0.6
        assert m == pytest.approx(0.6)


# -- compute_cross_batch_metrics --------------------------------------------


class TestCrossBatchMetrics:
    """Reduction from per-(cohort, variant) rows to per-variant cross-batch."""

    def _frame(self, rows: list[dict[str, float | str]]) -> pd.DataFrame:
        return pd.DataFrame(rows)

    def test_empty_input_returns_empty_frame_with_columns(self) -> None:
        """Edge case: no cohorts."""
        empty = pd.DataFrame(columns=['cohort_id', 'ID', 'F_MISS', 'ALT_FREQS', 'OBS_CT', 'HWE_P'])
        out = compute_cross_batch_metrics(empty, hwe_p_min=1e-6)
        assert list(out.columns) == ['ID', 'F_MISS_range', 'AF_chi2_P', 'HWE_flip_count', 'n_batches']
        assert len(out) == 0

    def test_call_rate_range(self) -> None:
        """F_MISS_range = max - min across batches."""
        cohort_qc = self._frame(
            [
                {'cohort_id': 'A', 'ID': 'rs1', 'F_MISS': 0.01, 'ALT_FREQS': 0.5, 'OBS_CT': 100, 'HWE_P': 0.9},
                {'cohort_id': 'B', 'ID': 'rs1', 'F_MISS': 0.10, 'ALT_FREQS': 0.5, 'OBS_CT': 100, 'HWE_P': 0.9},
            ]
        )
        out = compute_cross_batch_metrics(cohort_qc, hwe_p_min=1e-6)
        assert out.loc[0, 'F_MISS_range'] == pytest.approx(0.09)
        assert out.loc[0, 'n_batches'] == 2

    def test_hwe_flip_count(self) -> None:
        """Per-batch HWE pass/fail under the threshold are counted."""
        cohort_qc = self._frame(
            [
                {'cohort_id': 'A', 'ID': 'rs1', 'F_MISS': 0.0, 'ALT_FREQS': 0.5, 'OBS_CT': 100, 'HWE_P': 1e-9},
                {'cohort_id': 'B', 'ID': 'rs1', 'F_MISS': 0.0, 'ALT_FREQS': 0.5, 'OBS_CT': 100, 'HWE_P': 1e-9},
                {'cohort_id': 'C', 'ID': 'rs1', 'F_MISS': 0.0, 'ALT_FREQS': 0.5, 'OBS_CT': 100, 'HWE_P': 0.9},
            ]
        )
        out = compute_cross_batch_metrics(cohort_qc, hwe_p_min=1e-6)
        assert out.loc[0, 'HWE_flip_count'] == 2

    def test_af_heterogeneity_picked_up(self) -> None:
        """Divergent AF across batches -> small p-value."""
        cohort_qc = self._frame(
            [
                {'cohort_id': 'A', 'ID': 'rs1', 'F_MISS': 0.0, 'ALT_FREQS': 0.05, 'OBS_CT': 500, 'HWE_P': 0.9},
                {'cohort_id': 'B', 'ID': 'rs1', 'F_MISS': 0.0, 'ALT_FREQS': 0.95, 'OBS_CT': 500, 'HWE_P': 0.9},
            ]
        )
        out = compute_cross_batch_metrics(cohort_qc, hwe_p_min=1e-6)
        assert out.loc[0, 'AF_chi2_P'] < 1e-100


# -- compute_cluster_icc -----------------------------------------------------


class TestClusterIcc:
    """Per-variant pooling across cohorts for THETA cluster moments."""

    def test_per_cluster_columns_present(self) -> None:
        """Output has pooled_mean and ICC per cluster."""
        stats = pd.DataFrame(
            [
                {
                    'cohort_id': 'A',
                    'ID': 'rs1',
                    'n_AA': 10,
                    'mean_THETA_AA': 0.1,
                    'var_THETA_AA': 0.001,
                    'n_AB': 10,
                    'mean_THETA_AB': 0.5,
                    'var_THETA_AB': 0.001,
                    'n_BB': 10,
                    'mean_THETA_BB': 0.9,
                    'var_THETA_BB': 0.001,
                },
                {
                    'cohort_id': 'B',
                    'ID': 'rs1',
                    'n_AA': 10,
                    'mean_THETA_AA': 0.1,
                    'var_THETA_AA': 0.001,
                    'n_AB': 10,
                    'mean_THETA_AB': 0.5,
                    'var_THETA_AB': 0.001,
                    'n_BB': 10,
                    'mean_THETA_BB': 0.9,
                    'var_THETA_BB': 0.001,
                },
            ]
        )
        out = compute_cluster_icc(stats)
        for cluster in ('AA', 'AB', 'BB'):
            assert f'pooled_mean_THETA_{cluster}' in out.columns
            assert f'ICC_THETA_{cluster}' in out.columns
            assert f'pooled_var_within_THETA_{cluster}' in out.columns

    def test_single_cohort_returns_nan_icc(self) -> None:
        """One cohort -> can't compute ICC."""
        stats = pd.DataFrame(
            [
                {
                    'cohort_id': 'A',
                    'ID': 'rs1',
                    'n_AA': 10,
                    'mean_THETA_AA': 0.1,
                    'var_THETA_AA': 0.001,
                    'n_AB': 10,
                    'mean_THETA_AB': 0.5,
                    'var_THETA_AB': 0.001,
                    'n_BB': 10,
                    'mean_THETA_BB': 0.9,
                    'var_THETA_BB': 0.001,
                },
            ]
        )
        out = compute_cluster_icc(stats)
        for cluster in ('AA', 'AB', 'BB'):
            assert math.isnan(out.loc[0, f'ICC_THETA_{cluster}'])

    def test_empty_input(self) -> None:
        """Empty cluster_stats -> empty frame, columns still defined."""
        out = compute_cluster_icc(pd.DataFrame())
        assert len(out) == 0
        assert 'ICC_THETA_AA' in out.columns


# -- apply_filters -----------------------------------------------------------


def _audit_row(**overrides: float | str) -> dict[str, float | str]:
    """One in-bounds variant row; overrides drive specific failures."""
    base: dict[str, float | str] = {
        'ID': 'rs1',
        'CHROM': '1',
        'POS': 1000,
        'REF': 'A',
        'ALT': 'G',
        'GenTrain_Score': 0.95,
        'Cluster_Sep': 0.8,
        'N_AA_egt': 10,
        'N_AB_egt': 10,
        'N_BB_egt': 10,
        'meanTHETA_AA_egt': 0.1,
        'meanTHETA_AB_egt': 0.5,
        'meanTHETA_BB_egt': 0.9,
        'devTHETA_AA_egt': 0.03,
        'devTHETA_AB_egt': 0.03,
        'devTHETA_BB_egt': 0.03,
        'F_MISS': 0.001,
        'ALT_FREQS': 0.3,
        'OBS_CT': 500,
        'HWE_P': 0.5,
        'F_MISS_range': 0.01,
        'AF_chi2_P': 0.5,
        'HWE_flip_count': 0,
        'n_batches': 3,
        'pooled_mean_THETA_AA': 0.1,
        'pooled_mean_THETA_AB': 0.5,
        'pooled_mean_THETA_BB': 0.9,
        'ICC_THETA_AA': 0.05,
        'ICC_THETA_AB': 0.05,
        'ICC_THETA_BB': 0.05,
        'pooled_var_within_THETA_AA': 0.001,
        'pooled_var_within_THETA_AB': 0.001,
        'pooled_var_within_THETA_BB': 0.001,
    }
    base.update(overrides)
    return base


class TestApplyFilters:
    """Threshold logic; one variant per test, one failure mode at a time."""

    def _apply(self, **overrides: float | str) -> pd.Series:
        df = pd.DataFrame([_audit_row(**overrides)])
        return apply_filters(df, **DEFAULT_THRESHOLDS).iloc[0]

    def test_clean_variant_not_excluded(self) -> None:
        """A variant within all bounds passes."""
        row = self._apply()
        assert row['exclude'] is np.False_ or row['exclude'] is False
        assert row['exclusion_reasons'] == ''

    def test_low_gentrain_fails(self) -> None:
        row = self._apply(GenTrain_Score=0.5)
        assert row['fail_gentrain']
        assert row['exclude']
        assert 'fail_gentrain' in row['exclusion_reasons']

    def test_low_cluster_sep_fails(self) -> None:
        row = self._apply(Cluster_Sep=0.2)
        assert row['fail_cluster_sep']
        assert row['exclude']

    def test_empty_egt_cluster_flags_but_does_not_exclude(self) -> None:
        """Empty training cluster is informational, not a fail."""
        row = self._apply(N_AB_egt=0)
        assert row['flag_empty_egt_cluster']
        assert not row['exclude']

    def test_high_vmiss_fails(self) -> None:
        row = self._apply(F_MISS=0.05)
        assert row['fail_vmiss']

    def test_low_hwe_p_fails(self) -> None:
        row = self._apply(HWE_P=1e-9)
        assert row['fail_hwe']

    def test_monomorphic_flags_but_does_not_exclude(self) -> None:
        row = self._apply(ALT_FREQS=0.0)
        assert row['flag_monomorphic']
        assert not row['exclude']

    def test_call_rate_range_fails(self) -> None:
        row = self._apply(F_MISS_range=0.1)
        assert row['fail_cross_batch_call_rate']

    def test_af_chi2_fails(self) -> None:
        row = self._apply(AF_chi2_P=1e-9)
        assert row['fail_cross_batch_af']

    def test_hwe_flip_count_flags_but_does_not_exclude(self) -> None:
        row = self._apply(HWE_flip_count=2)
        assert row['flag_hwe_flip']
        # flag_hwe_flip is not in FAIL_COLS -> shouldn't drive exclude
        assert not row['exclude']

    def test_multiple_failures_concatenated(self) -> None:
        """Reasons string lists every contributing fail."""
        row = self._apply(GenTrain_Score=0.5, F_MISS=0.5)
        reasons = row['exclusion_reasons'].split(',')
        assert 'fail_gentrain' in reasons
        assert 'fail_vmiss' in reasons

    def test_nan_metric_does_not_fail(self) -> None:
        """A missing metric must not trigger its fail column."""
        row = self._apply(GenTrain_Score=np.nan)
        assert not row['fail_gentrain']

    def test_egt_centroid_delta_fails(self) -> None:
        """Observed AA cluster centroid far from EGT expected -> fail."""
        row = self._apply(pooled_mean_THETA_AA=0.5, meanTHETA_AA_egt=0.1)
        assert row['egt_centroid_delta_max'] == pytest.approx(0.4)
        assert row['fail_egt_centroid_delta']
        assert row['exclude']

    def test_egt_centroid_delta_passes_when_close(self) -> None:
        """Small deltas across all clusters do not trigger fail."""
        row = self._apply(pooled_mean_THETA_AA=0.12)  # delta 0.02, under 0.1 default
        assert not row['fail_egt_centroid_delta']

    def test_cluster_icc_fails(self) -> None:
        """Cross-batch ICC exceeds the cap on any cluster -> fail."""
        row = self._apply(ICC_THETA_AB=0.9)
        assert row['cluster_icc_max'] == pytest.approx(0.9)
        assert row['fail_cluster_icc']
        assert row['exclude']

    def test_cluster_spread_fails(self) -> None:
        """Smeared within-cluster THETA variance on any cluster -> fail."""
        row = self._apply(pooled_var_within_THETA_BB=0.05)
        assert row['cluster_spread_max'] == pytest.approx(0.05)
        assert row['fail_cluster_spread']
        assert row['exclude']

    def test_cluster_separation_fails(self) -> None:
        """Collapsed adjacent cluster gap -> fail."""
        # AB at 0.55, BB at 0.6 -> gap = 0.05 < 0.15 default
        row = self._apply(pooled_mean_THETA_AB=0.55, pooled_mean_THETA_BB=0.6)
        assert row['cluster_separation_min'] == pytest.approx(0.05)
        assert row['fail_cluster_separation']
        assert row['exclude']

    def test_cluster_separation_uses_min_of_two_gaps(self) -> None:
        """fail triggers on the smaller of (AB-AA, BB-AB)."""
        # AA=0.05 AB=0.5 BB=0.55 -> gaps 0.45 and 0.05, min=0.05
        row = self._apply(pooled_mean_THETA_AA=0.05, pooled_mean_THETA_AB=0.5, pooled_mean_THETA_BB=0.55)
        assert row['cluster_separation_min'] == pytest.approx(0.05)
        assert row['fail_cluster_separation']

    def test_intensity_metrics_nan_when_inputs_missing(self) -> None:
        """Missing pooled/EGT inputs -> NaN metric, no fail."""
        row = self._apply(
            pooled_mean_THETA_AA=np.nan,
            pooled_mean_THETA_AB=np.nan,
            pooled_mean_THETA_BB=np.nan,
            ICC_THETA_AA=np.nan,
            ICC_THETA_AB=np.nan,
            ICC_THETA_BB=np.nan,
            pooled_var_within_THETA_AA=np.nan,
            pooled_var_within_THETA_AB=np.nan,
            pooled_var_within_THETA_BB=np.nan,
        )
        assert not row['fail_egt_centroid_delta']
        assert not row['fail_cluster_icc']
        assert not row['fail_cluster_spread']
        assert not row['fail_cluster_separation']


# -- build_audit_frame ------------------------------------------------------


class TestBuildAuditFrame:
    """Wide-join behaviour across the four sources."""

    def test_left_join_preserves_egt_rows_with_missing_qc(self) -> None:
        """Variants only in EGT still appear in the audit."""
        egt = pd.DataFrame(
            [
                {
                    'ID': 'rs1',
                    'CHROM': '1',
                    'POS': 100,
                    'REF': 'A',
                    'ALT': 'G',
                    'GenTrain_Score': 0.9,
                    'Cluster_Sep': 0.6,
                    'N_AA_egt': 10,
                    'N_AB_egt': 10,
                    'N_BB_egt': 10,
                }
            ]
        )
        merged = pd.DataFrame(columns=['ID', 'F_MISS', 'ALT_FREQS', 'OBS_CT', 'HWE_P'])
        cross_batch = pd.DataFrame(columns=['ID', 'F_MISS_range', 'AF_chi2_P', 'HWE_flip_count', 'n_batches'])
        icc = pd.DataFrame(columns=['ID'])
        out = build_audit_frame(egt, merged, cross_batch, icc)
        assert 'rs1' in out['ID'].to_numpy()


# -- FAIL_COLS contract -----------------------------------------------------


def test_fail_cols_drive_exclusion() -> None:
    """Every column in FAIL_COLS must end up surfaced in exclusion_reasons."""
    row = _audit_row(
        GenTrain_Score=0.0,
        Cluster_Sep=0.0,
        F_MISS=1.0,
        HWE_P=0.0,
        F_MISS_range=1.0,
        AF_chi2_P=0.0,
        # Intensity-space fails (deltas / ICC / spread / collapsed separation).
        pooled_mean_THETA_AA=0.9,  # delta vs egt 0.1 -> 0.8
        ICC_THETA_AA=0.99,
        pooled_var_within_THETA_AA=1.0,
        pooled_mean_THETA_AB=0.91,  # gap AB - AA = 0.01 < 0.15
        pooled_mean_THETA_BB=0.92,  # gap BB - AB = 0.01 < 0.15
    )
    out = apply_filters(pd.DataFrame([row]), **DEFAULT_THRESHOLDS)
    reasons = set(out.loc[0, 'exclusion_reasons'].split(','))
    assert set(FAIL_COLS) <= reasons


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
