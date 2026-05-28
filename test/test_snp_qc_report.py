"""Unit tests for ``scripts/snp_qc_report.py``."""

from __future__ import annotations

import gzip
from typing import TYPE_CHECKING

import pandas as pd
import pytest

from popgen_genotyping.scripts.snp_qc_report import (
    apply_filters,
    load_egt_info,
    load_hwe_pass_ids,
    load_vmiss,
    main,
    summarise,
    write_outputs,
)

if TYPE_CHECKING:
    from pathlib import Path

CONSERVATIVE_THRESHOLDS: dict[str, float] = {
    'gentrain_min': 0.7,
    'cluster_sep_min': 0.4,
    'fmiss_max': 0.02,
}


@pytest.fixture
def egt_tsv(tmp_path: Path) -> Path:
    """Write a header-less EGT INFO TSV mimicking `bcftools query` output."""
    rows: list[str] = [
        # ID, alleles, GenTrain, Cluster_Sep
        'chr1\t100\tA\tG\trs1\t0.95\t0.80',  # pass everything
        'chr1\t200\tA\tT\trs2\t0.90\t0.70',  # pass (strand-ambiguity not filtered here)
        'chr1\t300\tC\tG\trs3\t0.85\t0.60',  # pass (strand-ambiguity not filtered here)
        'chr2\t100\tC\tT\trs4\t0.65\t0.50',  # fail gentrain
        'chr2\t200\tG\tA\trs5\t0.80\t0.30',  # fail cluster_sep
        'chr2\t300\tT\tC\trs6\t.\t.',  # missing scores → fail
        'chr3\t100\tCAG\tC\trs7\t0.90\t0.55',  # indel
        'chr3\t200\tA\tG\trs8\t0.95\t0.70',  # passes all other filters, fails HWE
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
        'chr3\trs8\t0\t100\t0',
    ]
    p = tmp_path / 'merged.vmiss'
    p.write_text('\n'.join(lines) + '\n')
    return p


@pytest.fixture
def hwe_pass_snplist(tmp_path: Path) -> Path:
    """Write a ``plink2 --write-snplist`` output omitting rs8 to mark it as HWE-fail."""
    ids: list[str] = ['rs1', 'rs2', 'rs3', 'rs4', 'rs5', 'rs6', 'rs7']
    p = tmp_path / 'hwe_pass.snplist'
    p.write_text('\n'.join(ids) + '\n')
    return p


def test_load_egt_info_parses_dot_as_nan(egt_tsv: Path) -> None:
    df = load_egt_info(egt_tsv)
    assert df.shape == (8, 7)
    assert df.loc[df['ID'] == 'rs6', 'GenTrain_Score'].isna().all()
    assert df.loc[df['ID'] == 'rs6', 'Cluster_Sep'].isna().all()


def test_load_vmiss_strips_header_hash(vmiss_tsv: Path) -> None:
    df = load_vmiss(vmiss_tsv)
    assert list(df.columns) == ['ID', 'F_MISS']
    assert df.loc[df['ID'] == 'rs7', 'F_MISS'].iloc[0] == pytest.approx(0.10)


def test_load_hwe_pass_ids_skips_blanks(tmp_path: Path) -> None:
    p = tmp_path / 'hwe.snplist'
    p.write_text('rs1\n\nrs2\n   \nrs3\n')
    assert load_hwe_pass_ids(p) == {'rs1', 'rs2', 'rs3'}


def test_apply_filters_marks_each_failure_class(
    egt_tsv: Path,
    vmiss_tsv: Path,
    hwe_pass_snplist: Path,
) -> None:
    df_egt = load_egt_info(egt_tsv)
    df_vmiss = load_vmiss(vmiss_tsv)
    df = df_egt.merge(df_vmiss, on='ID', how='left')
    out = apply_filters(df, hwe_pass_ids=load_hwe_pass_ids(hwe_pass_snplist), **CONSERVATIVE_THRESHOLDS)
    by_id: dict[str, dict[str, bool]] = out.set_index('ID').to_dict(orient='index')  # type: ignore[assignment]
    assert by_id['rs1']['fail'] is False
    assert by_id['rs2']['fail'] is False  # strand-ambiguity is not a per-SNP QC filter
    assert by_id['rs3']['fail'] is False  # strand-ambiguity is not a per-SNP QC filter
    assert by_id['rs4']['fail'] is True  # gentrain
    assert by_id['rs5']['fail'] is True  # cluster_sep
    assert by_id['rs6']['fail'] is True  # NaN scores
    assert by_id['rs7']['fail'] is True  # fmiss
    assert by_id['rs8']['fail'] is True  # absent from hwe pass list
    assert by_id['rs8']['fail_hwe'] is True
    assert by_id['rs1']['fail_hwe'] is False


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
    out = apply_filters(df, hwe_pass_ids={'x'}, **CONSERVATIVE_THRESHOLDS)
    assert bool(out['fail_fmiss'].iloc[0])
    assert bool(out['fail'].iloc[0])
    assert not bool(out['fail_hwe'].iloc[0])


def test_summarise_counts_match_filter_columns(
    egt_tsv: Path,
    vmiss_tsv: Path,
    hwe_pass_snplist: Path,
) -> None:
    df_egt = load_egt_info(egt_tsv)
    df_vmiss = load_vmiss(vmiss_tsv)
    df = df_egt.merge(df_vmiss, on='ID', how='left')
    out = apply_filters(df, hwe_pass_ids=load_hwe_pass_ids(hwe_pass_snplist), **CONSERVATIVE_THRESHOLDS)
    s = summarise(out).set_index('metric')['value'].to_dict()
    assert s['total_variants'] == 8
    assert s['fail_gentrain'] == 2  # rs4 + rs6 (NaN)
    assert s['fail_cluster_sep'] == 2  # rs5 + rs6 (NaN)
    assert s['fail_fmiss'] == 1  # rs7
    assert s['fail_hwe'] == 1  # rs8
    assert s['total_excluded'] + s['total_retained'] == s['total_variants']
    # Priority cascade gentrain → cluster_sep → fmiss → hwe: rs6 fails gentrain
    # and cluster_sep but is attributed only to gentrain.
    assert s['first_fail_gentrain'] == 2  # rs4, rs6
    assert s['first_fail_cluster_sep'] == 1  # rs5 (rs6 absorbed by gentrain)
    assert s['first_fail_fmiss'] == 1  # rs7
    assert s['first_fail_hwe'] == 1  # rs8
    cascade_keys = ('first_fail_gentrain', 'first_fail_cluster_sep', 'first_fail_fmiss', 'first_fail_hwe')
    assert sum(s[k] for k in cascade_keys) == s['total_excluded']


def test_main_end_to_end(
    tmp_path: Path,
    egt_tsv: Path,
    vmiss_tsv: Path,
    hwe_pass_snplist: Path,
) -> None:
    audit = tmp_path / 'snp_qc.audit.tsv.gz'
    excl = tmp_path / 'snp_qc.exclude.snplist'
    summary = tmp_path / 'snp_qc.summary.tsv'
    rc = main(
        [
            '--egt-info-tsv',
            str(egt_tsv),
            '--merged-vmiss',
            str(vmiss_tsv),
            '--hwe-pass-snplist',
            str(hwe_pass_snplist),
            '--gentrain-min',
            '0.7',
            '--cluster-sep-min',
            '0.4',
            '--fmiss-max',
            '0.02',
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
    assert set(audit_df['ID']) == {'rs1', 'rs2', 'rs3', 'rs4', 'rs5', 'rs6', 'rs7', 'rs8'}
    excluded: list[str] = excl.read_text().strip().splitlines()
    assert set(excluded) == {'rs4', 'rs5', 'rs6', 'rs7', 'rs8'}
    assert not audit_df.loc[audit_df['ID'] == 'rs1', 'fail'].iloc[0]


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
            'fail_gentrain': [False],
            'fail_cluster_sep': [False],
            'fail_fmiss': [False],
            'fail_hwe': [False],
            'fail': [False],
        },
    )
    audit = tmp_path / 'a.tsv.gz'
    excl = tmp_path / 'e.snplist'
    summary = tmp_path / 's.tsv'
    write_outputs(df, audit_path=audit, exclusion_path=excl, summary_path=summary)
    assert excl.read_text() == ''
