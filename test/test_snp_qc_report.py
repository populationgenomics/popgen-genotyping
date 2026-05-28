"""Unit tests for ``scripts/snp_qc_report.py``."""

from __future__ import annotations

import gzip
from typing import TYPE_CHECKING

import pandas as pd
import pytest

from popgen_genotyping.scripts.snp_qc_report import (
    add_strand_ambiguous_flag,
    apply_filters,
    load_egt_info,
    load_vmiss,
    main,
    summarise,
    write_outputs,
)

if TYPE_CHECKING:
    from pathlib import Path

CONSERVATIVE_THRESHOLDS: dict[str, float | bool] = {
    'gentrain_min': 0.7,
    'cluster_sep_min': 0.4,
    'fmiss_max': 0.02,
    'exclude_strand_ambiguous': True,
}


@pytest.fixture
def egt_tsv(tmp_path: Path) -> Path:
    """Write a header-less EGT INFO TSV mimicking `bcftools query` output."""
    rows: list[str] = [
        # ID, alleles, GenTrain, Cluster_Sep
        'chr1\t100\tA\tG\trs1\t0.95\t0.80',  # pass everything
        'chr1\t200\tA\tT\trs2\t0.90\t0.70',  # strand-ambiguous (A/T)
        'chr1\t300\tC\tG\trs3\t0.85\t0.60',  # strand-ambiguous (C/G)
        'chr2\t100\tC\tT\trs4\t0.65\t0.50',  # fail gentrain
        'chr2\t200\tG\tA\trs5\t0.80\t0.30',  # fail cluster_sep
        'chr2\t300\tT\tC\trs6\t.\t.',  # missing scores → fail
        'chr3\t100\tCAG\tC\trs7\t0.90\t0.55',  # indel, not strand-ambiguous
    ]
    p = tmp_path / 'egt.tsv'
    p.write_text('\n'.join(rows) + '\n')
    return p


@pytest.fixture
def vmiss_tsv(tmp_path: Path) -> Path:
    """Write a `plink2 --missing` per-variant TSV with the `#CHROM` header."""
    lines: list[str] = [
        '#CHROM\tID\tMISSING_CT\tOBS_CT\tF_MISS',
        'chr1\trs1\t0\t100\t0',
        'chr1\trs2\t1\t100\t0.01',
        'chr1\trs3\t1\t100\t0.01',
        'chr2\trs4\t1\t100\t0.01',
        'chr2\trs5\t1\t100\t0.01',
        'chr2\trs6\t1\t100\t0.01',
        'chr3\trs7\t10\t100\t0.10',  # fail fmiss
    ]
    p = tmp_path / 'merged.vmiss'
    p.write_text('\n'.join(lines) + '\n')
    return p


def test_load_egt_info_parses_dot_as_nan(egt_tsv: Path) -> None:
    df = load_egt_info(egt_tsv)
    assert df.shape == (7, 7)
    assert df.loc[df['ID'] == 'rs6', 'GenTrain_Score'].isna().all()
    assert df.loc[df['ID'] == 'rs6', 'Cluster_Sep'].isna().all()


def test_load_vmiss_strips_header_hash(vmiss_tsv: Path) -> None:
    df = load_vmiss(vmiss_tsv)
    assert list(df.columns) == ['ID', 'F_MISS']
    assert df.loc[df['ID'] == 'rs7', 'F_MISS'].iloc[0] == pytest.approx(0.10)


def test_strand_ambiguous_flag_classifies_alleles() -> None:
    df = pd.DataFrame(
        {
            'REF': ['A', 'A', 'T', 'C', 'C', 'G', 'A', 'CAG'],
            'ALT': ['G', 'T', 'A', 'T', 'G', 'C', 'C', 'C'],
        },
    )
    out = add_strand_ambiguous_flag(df)
    assert out['strand_ambiguous'].tolist() == [False, True, True, False, True, True, False, False]


def test_strand_ambiguous_handles_missing_alleles() -> None:
    df = pd.DataFrame({'REF': [None, 'A'], 'ALT': ['T', None]})
    out = add_strand_ambiguous_flag(df)
    assert out['strand_ambiguous'].tolist() == [False, False]


def test_apply_filters_marks_each_failure_class(egt_tsv: Path, vmiss_tsv: Path) -> None:
    df_egt = load_egt_info(egt_tsv)
    df_vmiss = load_vmiss(vmiss_tsv)
    df = add_strand_ambiguous_flag(df_egt.merge(df_vmiss, on='ID', how='left'))
    out = apply_filters(df, **CONSERVATIVE_THRESHOLDS)  # type: ignore[arg-type]
    by_id: dict[str, dict[str, bool]] = out.set_index('ID').to_dict(orient='index')  # type: ignore[assignment]
    assert by_id['rs1']['pass'] is True
    assert by_id['rs2']['pass'] is False  # A/T
    assert by_id['rs3']['pass'] is False  # C/G
    assert by_id['rs4']['pass'] is False  # gentrain
    assert by_id['rs5']['pass'] is False  # cluster_sep
    assert by_id['rs6']['pass'] is False  # NaN scores
    assert by_id['rs7']['pass'] is False  # fmiss


def test_apply_filters_strand_flag_off_keeps_ambiguous(egt_tsv: Path, vmiss_tsv: Path) -> None:
    df_egt = load_egt_info(egt_tsv)
    df_vmiss = load_vmiss(vmiss_tsv)
    df = add_strand_ambiguous_flag(df_egt.merge(df_vmiss, on='ID', how='left'))
    out = apply_filters(
        df,
        gentrain_min=0.7,
        cluster_sep_min=0.4,
        fmiss_max=0.02,
        exclude_strand_ambiguous=False,
    )
    by_id: dict[str, dict[str, bool]] = out.set_index('ID').to_dict(orient='index')  # type: ignore[assignment]
    assert by_id['rs2']['pass'] is True
    assert by_id['rs3']['pass'] is True
    # The bookkeeping column still tracks ambiguity even when not enforced.
    assert by_id['rs2']['strand_ambiguous'] is True


def test_apply_filters_nan_fmiss_fails() -> None:
    df = pd.DataFrame(
        {
            'REF': ['A'],
            'ALT': ['G'],
            'ID': ['x'],
            'GenTrain_Score': [0.9],
            'Cluster_Sep': [0.6],
            'F_MISS': [float('nan')],
        },
    )
    df = add_strand_ambiguous_flag(df)
    out = apply_filters(df, **CONSERVATIVE_THRESHOLDS)  # type: ignore[arg-type]
    assert out['pass_fmiss'].iloc[0] is False or not out['pass_fmiss'].iloc[0]
    assert not out['pass'].iloc[0]


def test_summarise_counts_match_filter_columns(egt_tsv: Path, vmiss_tsv: Path) -> None:
    df_egt = load_egt_info(egt_tsv)
    df_vmiss = load_vmiss(vmiss_tsv)
    df = add_strand_ambiguous_flag(df_egt.merge(df_vmiss, on='ID', how='left'))
    out = apply_filters(df, **CONSERVATIVE_THRESHOLDS)  # type: ignore[arg-type]
    s = summarise(out).set_index('metric')['value'].to_dict()
    assert s['total_variants'] == 7
    assert s['fail_gentrain'] == 2  # rs4 + rs6 (NaN)
    assert s['fail_cluster_sep'] == 2  # rs5 + rs6 (NaN)
    assert s['fail_fmiss'] == 1  # rs7
    assert s['strand_ambiguous'] == 2  # rs2 + rs3
    assert s['fail_strand_filter'] == 2
    assert s['total_excluded'] + s['total_retained'] == s['total_variants']
    # Priority cascade gentrain → cluster_sep → fmiss → strand: rs6 fails both
    # gentrain and cluster_sep but is attributed only to gentrain.
    assert s['first_fail_gentrain'] == 2  # rs4, rs6
    assert s['first_fail_cluster_sep'] == 1  # rs5 (rs6 absorbed by gentrain)
    assert s['first_fail_fmiss'] == 1  # rs7
    assert s['first_fail_strand'] == 2  # rs2, rs3
    cascade_keys = ('first_fail_gentrain', 'first_fail_cluster_sep', 'first_fail_fmiss', 'first_fail_strand')
    assert sum(s[k] for k in cascade_keys) == s['total_excluded']


def test_main_end_to_end(tmp_path: Path, egt_tsv: Path, vmiss_tsv: Path) -> None:
    audit = tmp_path / 'snp_qc.audit.tsv.gz'
    excl = tmp_path / 'snp_qc.exclude.snplist'
    summary = tmp_path / 'snp_qc.summary.tsv'
    rc = main(
        [
            '--egt-info-tsv',
            str(egt_tsv),
            '--merged-vmiss',
            str(vmiss_tsv),
            '--gentrain-min',
            '0.7',
            '--cluster-sep-min',
            '0.4',
            '--fmiss-max',
            '0.02',
            '--exclude-strand-ambiguous',
            '--output-audit-tsv',
            str(audit),
            '--output-exclusion-list',
            str(excl),
            '--output-summary-tsv',
            str(summary),
        ],
    )
    assert rc == 0
    with gzip.open(audit, 'rt') as f:
        audit_df = pd.read_csv(f, sep='\t')
    assert set(audit_df['ID']) == {'rs1', 'rs2', 'rs3', 'rs4', 'rs5', 'rs6', 'rs7'}
    excluded: list[str] = excl.read_text().strip().splitlines()
    assert set(excluded) == {'rs2', 'rs3', 'rs4', 'rs5', 'rs6', 'rs7'}
    assert audit_df.loc[audit_df['ID'] == 'rs1', 'pass'].iloc[0]


def test_write_outputs_empty_exclusion(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            'CHROM': ['chr1'],
            'POS': [100],
            'REF': ['A'],
            'ALT': ['G'],
            'ID': ['rs1'],
            'GenTrain_Score': [0.99],
            'Cluster_Sep': [0.9],
            'F_MISS': [0.0],
            'strand_ambiguous': [False],
            'pass_gentrain': [True],
            'pass_cluster_sep': [True],
            'pass_fmiss': [True],
            'pass_strand': [True],
            'pass': [True],
        },
    )
    audit = tmp_path / 'a.tsv.gz'
    excl = tmp_path / 'e.snplist'
    summary = tmp_path / 's.tsv'
    write_outputs(df, audit_path=audit, exclusion_path=excl, summary_path=summary)
    assert excl.read_text() == ''
