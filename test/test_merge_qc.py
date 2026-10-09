"""Tests for scripts/merge_qc.py."""

import csv
import gzip
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from popgen_genotyping.scripts.merge_qc import (
    IBIS_DEGREE_COL,
    IBIS_DEGREES,
    IBIS_KINSHIP_COL,
    IBIS_RELATIVES_COL,
    INFTYPE_BEYOND_DEGREE,
    KING_IBDSEG_INFTYPES,
    KING_IBDSEG_RELATIVES_COL,
    KING_ROBUST_DEGREE_THRESHOLDS,
    KING_ROBUST_RELATIVES_COL,
    RELATIVES_COLS,
    format_relatives,
    kinship_to_degree,
    main,
    process_bafregress,
    process_ibis_pairs,
    process_king_robust,
    process_relatedness_excluded,
    process_seg,
    read_qc_file,
)


# -- Fixtures ------------------------------------------------------------------


@pytest.fixture
def tmp_dir(tmp_path: Path) -> Path:
    """Return a temporary directory for test files."""
    return tmp_path


def _write(path: Path, text: str) -> str:
    """Write text to a file and return its string path."""
    path.write_text(text)
    return str(path)


def _write_gz(path: Path, text: str) -> str:
    """Write gzipped text to a file and return its string path."""
    with gzip.open(path, 'wt') as f:
        f.write(text)
    return str(path)


_SEG_HEADER = 'FID1\tID1\tFID2\tID2\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType\n'
_KIN0_HEADER = '#IID1\tIID2\tNSNP\tHETHET\tIBS0\tKINSHIP\n'
_IBIS_PAIRS_HEADER = (
    'ID1\tID2\tKING_INFTYPE\tKING_KINSHIP\tIBIS_KINSHIP_GE5CM\tIBIS_DEGREE_GE5CM\t'
    'IBIS_KINSHIP_GE10CM\tIBIS_DEGREE_GE10CM\tIBIS_N_SEG\n'
)


# -- read_qc_file -------------------------------------------------------------


class TestReadQcFile:
    """Tests for read_qc_file."""

    def test_basic_parsing(self, tmp_dir: Path) -> None:
        """Parse a simple whitespace-delimited file."""
        path = _write(
            tmp_dir / 'basic.txt',
            '#FID\tIID\tVAL\nFAM1\tS1\t0.5\nFAM2\tS2\t0.8\n',
        )
        df = read_qc_file(path)
        assert len(df) == 2
        assert list(df.columns) == ['FID', 'IID', 'VAL']
        assert df.iloc[0]['IID'] == 'S1'

    def test_strips_hash_from_headers(self, tmp_dir: Path) -> None:
        """Verify that leading '#' is stripped from column headers."""
        path = _write(tmp_dir / 'hash.txt', '#IID\t#METRIC\nS1\t1.0\n')
        df = read_qc_file(path)
        assert list(df.columns) == ['IID', 'METRIC']


# -- Column names and degree labels ---------------------------------------------


class TestRelativesColumns:
    """The three relatives columns and the labels their DEGREE field can take."""

    def test_one_column_per_estimator(self) -> None:
        """Three columns, one per estimator, in report order."""
        assert RELATIVES_COLS == ['KING_IBDSEG_RELATIVES', 'KING_ROBUST_RELATIVES', 'IBIS_RELATIVES_GE5CM']

    def test_known_inftypes(self) -> None:
        """KING --ibdseg InfType values kept as the DEGREE label."""
        assert {'Dup/MZ', 'PO', 'FS', '2nd', '3rd'} == KING_IBDSEG_INFTYPES

    def test_unknown_inftype(self) -> None:
        """Unknown InfType values (e.g. 4th, UN) are not kept."""
        assert '4th' not in KING_IBDSEG_INFTYPES
        assert 'UN' not in KING_IBDSEG_INFTYPES

    def test_beyond_degree_is_only_4th(self) -> None:
        """Only 4th is tolerated as beyond the requested degree; UN still raises."""
        assert {'4th'} == INFTYPE_BEYOND_DEGREE

    def test_ibis_degrees_match_king_robust_labels(self) -> None:
        """IBIS and KING-robust bin kinship into the same degree labels."""
        assert {label for _, label in KING_ROBUST_DEGREE_THRESHOLDS} == IBIS_DEGREES


# -- format_relatives -----------------------------------------------------------


class TestFormatRelatives:
    """Tests for format_relatives."""

    def test_pair_is_listed_under_both_samples(self) -> None:
        """Each member of a pair lists the other as REL_ID:KINSHIP:DEGREE."""
        pairs = pd.DataFrame({'ID1': ['S1'], 'ID2': ['S2'], 'KINSHIP': [0.25], 'DEGREE': ['1st']})
        result = format_relatives(pairs, 'COL').set_index('IID')['COL']
        assert result['S1'] == 'S2:0.25:1st'
        assert result['S2'] == 'S1:0.25:1st'

    def test_relatives_are_joined_with_semicolons(self) -> None:
        """A sample with several relatives gets them joined with ';', in pair order."""
        pairs = pd.DataFrame(
            {'ID1': ['S1', 'S3'], 'ID2': ['S2', 'S1'], 'KINSHIP': [0.25, 0.1], 'DEGREE': ['1st', '2nd']},
        )
        result = format_relatives(pairs, 'COL').set_index('IID')['COL']
        assert result['S1'] == 'S2:0.25:1st;S3:0.1:2nd'

    def test_no_pairs_gives_empty_frame_with_columns(self) -> None:
        """No pairs gives an empty frame with IID and the named column."""
        pairs = pd.DataFrame(columns=['ID1', 'ID2', 'KINSHIP', 'DEGREE'])
        result = format_relatives(pairs, 'COL')
        assert result.empty
        assert list(result.columns) == ['IID', 'COL']


# -- process_seg ---------------------------------------------------------------


class TestProcessSeg:
    """Tests for process_seg."""

    def test_basic_pair(self, tmp_dir: Path) -> None:
        """Process a .seg file with one full-sib pair."""
        # Full sibs: IBD1Seg ≈ 0.5, IBD2Seg ≈ 0.25 → kinship = 0.5/4 + 0.25/2 = 0.25
        path = _write(
            tmp_dir / 'test.seg',
            _SEG_HEADER + 'F\tS1\tF\tS2\t0.5000\t0.2500\t0.7500\tFS\n',
        )
        result = process_seg(path)
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ['IID', KING_IBDSEG_RELATIVES_COL]

        s1_row = result[result['IID'] == 'S1'].iloc[0]
        assert s1_row[KING_IBDSEG_RELATIVES_COL] == 'S2:0.25:FS'

        s2_row = result[result['IID'] == 'S2'].iloc[0]
        assert s2_row[KING_IBDSEG_RELATIVES_COL] == 'S1:0.25:FS'

    def test_parent_offspring(self, tmp_dir: Path) -> None:
        """PO pairs carry InfType=PO as their DEGREE."""
        # PO: IBD1Seg ≈ 1.0, IBD2Seg ≈ 0 → kinship = 0.25
        path = _write(
            tmp_dir / 'po.seg',
            _SEG_HEADER + 'F\tParent\tF\tChild\t1.0000\t0.0000\t0.5000\tPO\n',
        )
        result = process_seg(path)
        parent_row = result[result['IID'] == 'Parent'].iloc[0]
        assert parent_row[KING_IBDSEG_RELATIVES_COL] == 'Child:0.25:PO'

    def test_unknown_inftype_raises(self, tmp_dir: Path) -> None:
        """Unrecognised InfType values raise ValueError."""
        path = _write(
            tmp_dir / 'unknown.seg',
            _SEG_HEADER + 'F\tS1\tF\tS2\t0.0100\t0.0000\t0.0050\tUN\n',
        )
        with pytest.raises(ValueError, match='unrecognised InfType'):
            process_seg(path)

    def test_4th_degree_dropped(self, tmp_dir: Path, capsys: pytest.CaptureFixture[str]) -> None:
        """4th-degree pairs are dropped with a stderr note; other pairs are kept."""
        path = _write(
            tmp_dir / 'fourth.seg',
            _SEG_HEADER + 'F\tS1\tF\tS2\t0.5000\t0.2500\t0.7500\tFS\n' + 'F\tS3\tF\tS4\t0.1743\t0.0012\t0.0884\t4th\n',
        )
        result = process_seg(path)
        assert set(result['IID']) == {'S1', 'S2'}
        assert 'Dropping 1 pair(s)' in capsys.readouterr().err

    def test_only_4th_degree_returns_empty(self, tmp_dir: Path) -> None:
        """A .seg holding only 4th-degree pairs yields the empty relatedness frame."""
        path = _write(
            tmp_dir / 'only_fourth.seg',
            _SEG_HEADER + 'F\tS3\tF\tS4\t0.1743\t0.0012\t0.0884\t4th\n',
        )
        result = process_seg(path)
        assert result.empty
        assert list(result.columns) == ['IID', KING_IBDSEG_RELATIVES_COL]

    def test_missing_columns_raises(self, tmp_dir: Path) -> None:
        """Missing required columns raise ValueError."""
        path = _write(
            tmp_dir / 'bad.seg',
            'COL_A\tCOL_B\nX\tY\n',
        )
        with pytest.raises(ValueError, match='missing required column'):
            process_seg(path)

    def test_empty_seg(self, tmp_dir: Path) -> None:
        """Return empty DataFrame for a header-only .seg file."""
        path = _write(tmp_dir / 'empty.seg', _SEG_HEADER)
        result = process_seg(path)
        assert result.empty

    def test_multiple_degrees(self, tmp_dir: Path) -> None:
        """A sample with relationships at different degrees lists them all in the one column."""
        path = _write(
            tmp_dir / 'multi.seg',
            _SEG_HEADER + 'F\tS1\tF\tS2\t0.5000\t0.2500\t0.7500\tFS\n' + 'F\tS1\tF\tS3\t0.1800\t0.0000\t0.0900\t3rd\n',
        )
        result = process_seg(path)
        s1_row = result[result['IID'] == 'S1'].iloc[0]
        assert s1_row[KING_IBDSEG_RELATIVES_COL] == 'S2:0.25:FS;S3:0.045:3rd'


# -- kinship_to_degree ---------------------------------------------------------


class TestKinshipToDegree:
    """Tests for kinship_to_degree (Manichaikul et al. 2010 bounds)."""

    def test_interior_values(self) -> None:
        """A kinship inside each range gets that range's degree."""
        kinship = pd.Series([0.45, 0.25, 0.125, 0.06])
        assert kinship_to_degree(kinship).tolist() == ['Dup/MZ', '1st', '2nd', '3rd']

    def test_dup_mz_bound_is_exclusive(self) -> None:
        """Dup/MZ is strictly above 0.354; 0.354 itself is 1st degree."""
        kinship = pd.Series([0.3541, 0.354])
        assert kinship_to_degree(kinship).tolist() == ['Dup/MZ', '1st']

    def test_lower_bounds_are_inclusive(self) -> None:
        """Each lower degree includes its lower bound: 0.177 is 1st, 0.0884 is 2nd, 0.0442 is 3rd."""
        kinship = pd.Series([0.177, 0.1769, 0.0884, 0.0883, 0.0442])
        assert kinship_to_degree(kinship).tolist() == ['1st', '2nd', '2nd', '3rd', '3rd']

    def test_below_third_degree_has_no_degree(self) -> None:
        """Kinship below 0.0442, including negative values, is beyond 3rd degree."""
        kinship = pd.Series([0.0441, 0.0, -0.1])
        assert kinship_to_degree(kinship).isna().all()


# -- process_king_robust --------------------------------------------------------


class TestProcessKingRobust:
    """Tests for process_king_robust."""

    def test_basic_pair(self, tmp_dir: Path) -> None:
        """A .kin0 pair is listed under both samples with kinship (4dp) and its binned degree."""
        path = _write(tmp_dir / 'test.kin0', _KIN0_HEADER + 'S1\tS2\t800\t0.16875\t0\t0.213908\n')
        result = process_king_robust(path)
        assert list(result.columns) == ['IID', KING_ROBUST_RELATIVES_COL]
        values = result.set_index('IID')[KING_ROBUST_RELATIVES_COL]
        assert values['S1'] == 'S2:0.2139:1st'
        assert values['S2'] == 'S1:0.2139:1st'

    def test_reads_gzipped_table(self, tmp_dir: Path) -> None:
        """The stage writes the table gzipped; compression is inferred from the .gz name."""
        path = _write_gz(tmp_dir / 'test.kin0.gz', _KIN0_HEADER + 'S1\tS2\t800\t0.40375\t0\t0.5\n')
        values = process_king_robust(path).set_index('IID')[KING_ROBUST_RELATIVES_COL]
        assert values['S1'] == 'S2:0.5:Dup/MZ'

    def test_each_degree_is_binned(self, tmp_dir: Path) -> None:
        """Pairs across all four degrees get the matching labels."""
        path = _write(
            tmp_dir / 'test.kin0',
            _KIN0_HEADER
            + 'A\tB\t800\t0.4\t0\t0.49\n'
            + 'A\tC\t800\t0.2\t0\t0.25\n'
            + 'A\tD\t800\t0.1\t0\t0.125\n'
            + 'A\tE\t800\t0.05\t0\t0.06\n',
        )
        values = process_king_robust(path).set_index('IID')[KING_ROBUST_RELATIVES_COL]
        assert values['A'] == 'B:0.49:Dup/MZ;C:0.25:1st;D:0.125:2nd;E:0.06:3rd'

    def test_pair_below_third_degree_dropped(self, tmp_dir: Path, capsys: pytest.CaptureFixture[str]) -> None:
        """A pair below the 3rd-degree bound (a lower table floor) is dropped with a stderr note."""
        path = _write(
            tmp_dir / 'test.kin0',
            _KIN0_HEADER + 'S1\tS2\t800\t0.16875\t0\t0.213908\n' + 'S3\tS4\t800\t0.1\t0.07\t0.0301\n',
        )
        result = process_king_robust(path)
        assert set(result['IID']) == {'S1', 'S2'}
        assert 'Dropping 1 pair(s)' in capsys.readouterr().err

    def test_header_only_table_returns_empty(self, tmp_dir: Path) -> None:
        """A header-only .kin0 (no pair reached the kinship floor) yields the empty frame."""
        path = _write(tmp_dir / 'empty.kin0', _KIN0_HEADER)
        result = process_king_robust(path)
        assert result.empty
        assert list(result.columns) == ['IID', KING_ROBUST_RELATIVES_COL]

    def test_header_only_gzipped_table_returns_empty(self, tmp_dir: Path) -> None:
        """A header-only gzipped .kin0 yields the empty frame."""
        path = _write_gz(tmp_dir / 'empty.kin0.gz', _KIN0_HEADER)
        assert process_king_robust(path).empty

    def test_missing_columns_raises(self, tmp_dir: Path) -> None:
        """A .kin0 without IID1/IID2/KINSHIP raises ValueError, even with no rows."""
        path = _write(tmp_dir / 'bad.kin0', '#FID1\tIID1\tKINSHIP\nF\tS1\t0.5\n')
        with pytest.raises(ValueError, match=r'missing required column.*IID2'):
            process_king_robust(path)
        header_only = _write(tmp_dir / 'bad_empty.kin0', '#IID1\tIID2\tNSNP\n')
        with pytest.raises(ValueError, match=r'missing required column.*KINSHIP'):
            process_king_robust(header_only)

    def test_non_numeric_kinship_raises(self, tmp_dir: Path) -> None:
        """A non-numeric or missing KINSHIP raises ValueError rather than being binned."""
        path = _write(tmp_dir / 'bad.kin0', _KIN0_HEADER + 'S1\tS2\t800\t0.1\t0\tnan\n')
        with pytest.raises(ValueError, match='non-numeric KINSHIP'):
            process_king_robust(path)
        text = _write(tmp_dir / 'text.kin0', _KIN0_HEADER + 'S1\tS2\t800\t0.1\t0\tabc\n')
        with pytest.raises(ValueError, match='non-numeric KINSHIP'):
            process_king_robust(text)


# -- process_ibis_pairs ---------------------------------------------------------


class TestProcessIbisPairs:
    """Tests for process_ibis_pairs."""

    def test_basic_pair(self, tmp_dir: Path) -> None:
        """A related pair is listed under both samples with the 5 cM kinship (4dp) and degree."""
        path = _write_gz(
            tmp_dir / 'pairs.tsv.gz',
            _IBIS_PAIRS_HEADER + 'S1\tS2\tFS\t0.25\t0.246913\t1st\t0.24\t1st\t12\n',
        )
        result = process_ibis_pairs(path)
        assert list(result.columns) == ['IID', IBIS_RELATIVES_COL]
        values = result.set_index('IID')[IBIS_RELATIVES_COL]
        assert values['S1'] == 'S2:0.2469:1st'
        assert values['S2'] == 'S1:0.2469:1st'

    def test_reads_the_5cm_columns_only(self, tmp_dir: Path) -> None:
        """The 10 cM columns are ignored; a pair related at 5 cM but not 10 cM is still listed."""
        path = _write_gz(
            tmp_dir / 'pairs.tsv.gz',
            _IBIS_PAIRS_HEADER + 'S1\tS2\tNA\tNA\t0.0625\t3rd\t0.0\tNA\t3\n',
        )
        values = process_ibis_pairs(path).set_index('IID')[IBIS_RELATIVES_COL]
        assert values['S1'] == 'S2:0.0625:3rd'

    def test_every_degree_label_is_kept(self, tmp_dir: Path) -> None:
        """Dup/MZ, 1st, 2nd and 3rd are all kept."""
        path = _write_gz(
            tmp_dir / 'pairs.tsv.gz',
            _IBIS_PAIRS_HEADER
            + 'A\tB\tNA\tNA\t0.45\tDup/MZ\tNA\tNA\t1\n'
            + 'A\tC\tNA\tNA\t0.25\t1st\tNA\tNA\t1\n'
            + 'A\tD\tNA\tNA\t0.125\t2nd\tNA\tNA\t1\n'
            + 'A\tE\tNA\tNA\t0.06\t3rd\tNA\tNA\t1\n',
        )
        values = process_ibis_pairs(path).set_index('IID')[IBIS_RELATIVES_COL]
        assert values['A'] == 'B:0.45:Dup/MZ;C:0.25:1st;D:0.125:2nd;E:0.06:3rd'

    def test_pairs_without_a_degree_are_dropped(self, tmp_dir: Path, capsys: pytest.CaptureFixture[str]) -> None:
        """NA degree (below the 3rd-degree bound, or a KING pair with no IBIS run) is dropped with a note."""
        path = _write_gz(
            tmp_dir / 'pairs.tsv.gz',
            _IBIS_PAIRS_HEADER
            + 'S1\tS2\tFS\t0.25\t0.246913\t1st\t0.24\t1st\t12\n'
            + 'S3\tS4\t3rd\t0.05\t0.02\tNA\t0.0\tNA\t2\n'
            + 'S5\tS6\t3rd\t0.05\tNA\tNA\tNA\tNA\tNA\n',
        )
        result = process_ibis_pairs(path)
        assert set(result['IID']) == {'S1', 'S2'}
        assert 'Dropping 2 pair(s)' in capsys.readouterr().err

    def test_unknown_degree_raises(self, tmp_dir: Path) -> None:
        """A degree label that is neither a known degree nor NA raises ValueError."""
        path = _write_gz(
            tmp_dir / 'pairs.tsv.gz',
            _IBIS_PAIRS_HEADER + 'S1\tS2\tNA\tNA\t0.02\t4th\tNA\tNA\t1\n',
        )
        with pytest.raises(ValueError, match='unrecognised IBIS_DEGREE_GE5CM'):
            process_ibis_pairs(path)

    def test_lowercase_na_is_not_a_missing_degree(self, tmp_dir: Path) -> None:
        """Only the exact string NA marks a missing degree; any other spelling is an unknown label."""
        path = _write_gz(
            tmp_dir / 'pairs.tsv.gz',
            _IBIS_PAIRS_HEADER + 'S1\tS2\tNA\tNA\t0.02\tna\tNA\tNA\t1\n',
        )
        with pytest.raises(ValueError, match='unrecognised IBIS_DEGREE_GE5CM'):
            process_ibis_pairs(path)

    def test_related_pair_without_kinship_raises(self, tmp_dir: Path) -> None:
        """A degree with an NA or non-numeric kinship is a malformed table."""
        path = _write_gz(
            tmp_dir / 'pairs.tsv.gz',
            _IBIS_PAIRS_HEADER + 'S1\tS2\tNA\tNA\tNA\t1st\tNA\tNA\t1\n',
        )
        with pytest.raises(ValueError, match='non-numeric IBIS_KINSHIP_GE5CM'):
            process_ibis_pairs(path)

    def test_missing_columns_raises(self, tmp_dir: Path) -> None:
        """A table without the 5 cM kinship or degree column raises ValueError."""
        path = _write_gz(
            tmp_dir / 'pairs.tsv.gz',
            'ID1\tID2\tIBIS_KINSHIP_GE10CM\tIBIS_DEGREE_GE10CM\nS1\tS2\t0.25\t1st\n',
        )
        with pytest.raises(ValueError, match=rf'missing required column.*{IBIS_DEGREE_COL}.*{IBIS_KINSHIP_COL}'):
            process_ibis_pairs(path)

    def test_header_only_table_returns_empty(self, tmp_dir: Path) -> None:
        """A header-only table (no related pairs at all) yields the empty frame."""
        path = _write_gz(tmp_dir / 'pairs.tsv.gz', _IBIS_PAIRS_HEADER)
        result = process_ibis_pairs(path)
        assert result.empty
        assert list(result.columns) == ['IID', IBIS_RELATIVES_COL]


# -- process_bafregress --------------------------------------------------------


class TestProcessBafregress:
    """Tests for process_bafregress."""

    def test_basic_bafregress(self, tmp_dir: Path) -> None:
        """Read and concatenate bafregress files, renaming sample_id to IID."""
        path1 = _write(
            tmp_dir / 'baf1.txt',
            'sample_id\tLRR_mean\nS1\t0.01\n',
        )
        path2 = _write(
            tmp_dir / 'baf2.txt',
            'sample_id\tLRR_mean\nS2\t0.02\n',
        )
        result = process_bafregress([path1, path2])
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2
        assert 'IID' in result.columns
        assert list(result['IID']) == ['S1', 'S2']

    def test_missing_file(self, tmp_dir: Path) -> None:
        """Missing files produce a warning but don't raise."""
        result = process_bafregress([str(tmp_dir / 'nonexistent.txt')])
        assert result.empty

    def test_empty_list(self) -> None:
        """Empty input list returns empty DataFrame."""
        result = process_bafregress([])
        assert result.empty


# -- process_relatedness_excluded -----------------------------------------------


class TestProcessRelatednessExcluded:
    """Tests for process_relatedness_excluded."""

    def test_flags_only_excluded_iids_with_reason(self, tmp_dir: Path) -> None:
        """Only IIDs present in the exclusion TSV are flagged, with their REASON string."""
        path = _write(
            tmp_dir / 'excluded.tsv',
            'IID\tREASON\tBAF_REGRESS\tF_MISS\nS2\tcontamination\t0.1125\t0.01\n',
        )
        report_iids = pd.Series(['S1', 'S2', 'S3'])
        result = process_relatedness_excluded(path, report_iids)
        assert list(result.columns) == ['IID', 'RELATEDNESS_EXCLUDED']
        values = result.set_index('IID')['RELATEDNESS_EXCLUDED']
        assert pd.isna(values['S1'])
        assert values['S2'] == 'contamination'
        assert pd.isna(values['S3'])

    def test_combined_reason_string(self, tmp_dir: Path) -> None:
        """A sample excluded for both criteria carries the combined REASON string."""
        path = _write(
            tmp_dir / 'excluded.tsv',
            'IID\tREASON\tBAF_REGRESS\tF_MISS\nS1\tcontamination;missingness\t0.1125\t0.05\n',
        )
        report_iids = pd.Series(['S1'])
        result = process_relatedness_excluded(path, report_iids)
        assert result.set_index('IID')['RELATEDNESS_EXCLUDED']['S1'] == 'contamination;missingness'

    def test_header_only_file_flags_nothing(self, tmp_dir: Path) -> None:
        """A header-only exclusion TSV flags no report IID."""
        path = _write(tmp_dir / 'excluded.tsv', 'IID\tREASON\tBAF_REGRESS\tF_MISS\n')
        report_iids = pd.Series(['S1', 'S2'])
        result = process_relatedness_excluded(path, report_iids)
        assert result['RELATEDNESS_EXCLUDED'].isna().all()

    def test_excluded_iid_not_in_report_raises(self, tmp_dir: Path) -> None:
        """An excluded IID absent from the report raises ValueError."""
        path = _write(tmp_dir / 'excluded.tsv', 'IID\tREASON\tBAF_REGRESS\tF_MISS\nS99\tcontamination\t0.5\t0.01\n')
        report_iids = pd.Series(['S1', 'S2'])
        with pytest.raises(ValueError, match='not found in QC report'):
            process_relatedness_excluded(path, report_iids)


# -- end-to-end ----------------------------------------------------------------


_NO_RELATIVES_SEG = _SEG_HEADER
_NO_RELATIVES_KIN0 = _KIN0_HEADER
_NO_RELATIVES_IBIS = _IBIS_PAIRS_HEADER
_NOBODY_EXCLUDED = 'IID\tREASON\tBAF_REGRESS\tF_MISS\n'


def _run_main(
    tmp_dir: Path,
    *,
    seg: str = _NO_RELATIVES_SEG,
    kin0: str = _NO_RELATIVES_KIN0,
    ibis_pairs: str = _NO_RELATIVES_IBIS,
    excluded: str = _NOBODY_EXCLUDED,
    bafregress: str | None = None,
    samples: tuple[str, ...] = ('S1', 'S2'),
) -> tuple[list[dict[str, str]], list[str]]:
    """Write the QC inputs for the given samples, run main(), and return the CSV rows and header.

    Args:
        tmp_dir: Directory for the input and output files.
        seg: Text of the KING --ibdseg .seg file.
        kin0: Text of the KING-robust .kin0 file (written gzipped, as the stage does).
        ibis_pairs: Text of the IBIS pairs table (written gzipped, as the stage does).
        excluded: Text of the relatedness-excluded TSV.
        bafregress: Text of a bafregress file, or None to omit --bafregress.
        samples: Sample IDs in the PLINK2 QC files.

    Returns:
        The output CSV rows as dicts, and its column names.
    """
    families = {sample: f'FAM{i + 1}' for i, sample in enumerate(samples)}
    smiss_path = _write(
        tmp_dir / 'test.smiss',
        '#FID\tIID\tMISS_PHENO_CT\tMISSING_CT\tOBS_CT\tF_MISS\n'
        + ''.join(f'{fid}\t{iid}\t0\t10\t1000\t0.01\n' for iid, fid in families.items()),
    )
    het_path = _write(
        tmp_dir / 'test.het',
        '#FID\tIID\tO_HOM\tE_HOM\tN_NM\tF\n'
        + ''.join(f'{fid}\t{iid}\t500\t490\t1000\t0.02\n' for iid, fid in families.items()),
    )
    sexcheck_path = _write(
        tmp_dir / 'test.sexcheck',
        '#FID\tIID\tPEDSEX\tSNPSEX\tSTATUS\tF\n'
        + ''.join(f'{fid}\t{iid}\t1\t1\tOK\t0.99\n' for iid, fid in families.items()),
    )
    seg_path = _write(tmp_dir / 'test.seg', seg)
    kin0_path = _write_gz(tmp_dir / 'test.kin0.gz', kin0)
    ibis_path = _write_gz(tmp_dir / 'pairs.tsv.gz', ibis_pairs)
    excluded_path = _write(tmp_dir / 'excluded.tsv', excluded)
    output_path = str(tmp_dir / 'output.csv')

    sys_argv = [
        'merge_qc.py',
        '--smiss',
        smiss_path,
        '--het',
        het_path,
        '--sexcheck',
        sexcheck_path,
        '--seg',
        seg_path,
        '--kin0',
        kin0_path,
        '--ibis-pairs',
        ibis_path,
        '--relatedness-excluded',
        excluded_path,
        '--output',
        output_path,
    ]
    if bafregress is not None:
        sys_argv += ['--bafregress', _write(tmp_dir / 'baf.txt', bafregress)]

    original_argv = sys.argv
    try:
        sys.argv = sys_argv
        main()
    finally:
        sys.argv = original_argv

    with open(output_path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    return rows, list(fieldnames)


class TestEndToEnd:
    """End-to-end test calling main() with temp files."""

    def test_full_merge(self, tmp_dir: Path) -> None:
        """Run a full merge with all QC file types and verify CSV output."""
        rows, fieldnames = _run_main(
            tmp_dir,
            seg=_SEG_HEADER + 'FAM1\tS1\tFAM2\tS2\t0.5000\t0.2500\t0.7500\tFS\n',
            kin0=_KIN0_HEADER + 'S1\tS2\t800\t0.16875\t0\t0.213908\n',
            ibis_pairs=_IBIS_PAIRS_HEADER + 'S1\tS2\tFS\t0.25\t0.246913\t1st\t0.24\t1st\t12\n',
            bafregress='sample_id\tLRR_mean\nS1\t0.005\nS2\t0.006\n',
        )

        assert len(rows) == 2
        assert rows[0]['IID'] == 'S1'
        assert rows[1]['IID'] == 'S2'

        # One relatives column per estimator, replacing the old per-degree columns.
        for col in RELATIVES_COLS:
            assert col in rows[0]
        assert not any(col.startswith('RELATED_') for col in fieldnames)

        # Each estimator lists S2 as S1's relative (and vice versa) with its own kinship and degree.
        assert rows[0][KING_IBDSEG_RELATIVES_COL] == 'S2:0.25:FS'
        assert rows[1][KING_IBDSEG_RELATIVES_COL] == 'S1:0.25:FS'
        assert rows[0][KING_ROBUST_RELATIVES_COL] == 'S2:0.2139:1st'
        assert rows[1][KING_ROBUST_RELATIVES_COL] == 'S1:0.2139:1st'
        assert rows[0][IBIS_RELATIVES_COL] == 'S2:0.2469:1st'
        assert rows[1][IBIS_RELATIVES_COL] == 'S1:0.2469:1st'

        # Bafregress data should be merged
        assert 'LRR_mean' in rows[0]
        assert rows[0]['LRR_mean'] == '0.005'

        # Nobody was excluded from relatedness estimation in this run.
        assert rows[0]['RELATEDNESS_EXCLUDED'] == ''
        assert rows[1]['RELATEDNESS_EXCLUDED'] == ''

        # RELATEDNESS_EXCLUDED sits immediately after the three relatives columns.
        assert fieldnames.index(KING_ROBUST_RELATIVES_COL) == fieldnames.index(KING_IBDSEG_RELATIVES_COL) + 1
        assert fieldnames.index(IBIS_RELATIVES_COL) == fieldnames.index(KING_ROBUST_RELATIVES_COL) + 1
        assert fieldnames.index('RELATEDNESS_EXCLUDED') == fieldnames.index(IBIS_RELATIVES_COL) + 1

    def test_estimators_are_independent(self, tmp_dir: Path) -> None:
        """A pair one estimator calls does not appear in the other columns."""
        rows, _fieldnames = _run_main(
            tmp_dir,
            kin0=_KIN0_HEADER + 'S1\tS2\t800\t0.06\t0.01\t0.05\n',
        )
        assert rows[0][KING_ROBUST_RELATIVES_COL] == 'S2:0.05:3rd'
        assert rows[0][KING_IBDSEG_RELATIVES_COL] == ''
        assert rows[0][IBIS_RELATIVES_COL] == ''

    @pytest.mark.parametrize(
        ('estimator', 'inputs', 'relatives_col'),
        [
            (
                'KING --ibdseg',
                {'seg': _SEG_HEADER + 'FAM1\tS1\tFAM2\tS2\t0.5000\t0.2500\t0.7500\tFS\n'},
                KING_IBDSEG_RELATIVES_COL,
            ),
            (
                'KING-robust',
                {'kin0': _KIN0_HEADER + 'S1\tS2\t800\t0.16875\t0\t0.213908\n'},
                KING_ROBUST_RELATIVES_COL,
            ),
            (
                'IBIS',
                {'ibis_pairs': _IBIS_PAIRS_HEADER + 'S1\tS2\tNA\tNA\t0.246913\t1st\t0.24\t1st\t12\n'},
                IBIS_RELATIVES_COL,
            ),
        ],
    )
    def test_excluded_sample_with_relatives_raises(
        self, tmp_dir: Path, estimator: str, inputs: dict[str, Any], relatives_col: str
    ) -> None:
        """A sample flagged RELATEDNESS_EXCLUDED but still present in any estimator's output raises.

        This is the invariant that no estimator saw an excluded sample: the exclusion list
        and an estimator's output disagreeing means that estimator's --remove didn't actually
        drop the sample before it ran. Each of the three estimators is checked.
        """
        del estimator  # only labels the parametrised case
        with pytest.raises(ValueError, match=f'unexpectedly have relatives in {relatives_col}'):
            _run_main(
                tmp_dir,
                excluded='IID\tREASON\tBAF_REGRESS\tF_MISS\nS1\tcontamination\t0.1125\t0.01\n',
                **inputs,
            )

    def test_excluded_sample_flagged_with_reason(self, tmp_dir: Path) -> None:
        """An excluded sample with no relatives is flagged with its reason; others are not."""
        rows, _fieldnames = _run_main(
            tmp_dir,
            excluded='IID\tREASON\tBAF_REGRESS\tF_MISS\nS1\tcontamination;missingness\t0.1125\t0.05\n',
        )
        assert rows[0]['RELATEDNESS_EXCLUDED'] == 'contamination;missingness'
        assert rows[1]['RELATEDNESS_EXCLUDED'] == ''

    def test_excluded_iid_not_in_report_raises_from_main(self, tmp_dir: Path) -> None:
        """main() propagates the error when the exclusion list names an unknown IID."""
        with pytest.raises(ValueError, match='not found in QC report'):
            _run_main(
                tmp_dir,
                excluded='IID\tREASON\tBAF_REGRESS\tF_MISS\nS99\tcontamination\t0.5\t0.01\n',
                samples=('S1',),
            )

    def test_malformed_estimator_input_fails_the_report(self, tmp_dir: Path) -> None:
        """A malformed input of any estimator raises instead of emitting an incomplete report."""
        with pytest.raises(ValueError, match='unrecognised IBIS_DEGREE_GE5CM'):
            _run_main(tmp_dir, ibis_pairs=_IBIS_PAIRS_HEADER + 'S1\tS2\tNA\tNA\t0.02\t4th\tNA\tNA\t1\n')
