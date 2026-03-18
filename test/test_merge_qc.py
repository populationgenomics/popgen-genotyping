"""Tests for scripts/merge_qc.py."""

import csv
import sys
from pathlib import Path

import pandas as pd
import pytest

from popgen_genotyping.scripts.merge_qc import (
    DEGREE_COLS,
    get_degree,
    main,
    process_bafregress,
    process_kinship,
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


# -- get_degree ----------------------------------------------------------------


class TestGetDegree:
    """Tests for get_degree."""

    def test_mz(self) -> None:
        """Kinship >= 0.354 is MZ."""
        assert get_degree(0.5) == 'RELATED_MZ'
        assert get_degree(0.354) == 'RELATED_MZ'

    def test_first_degree(self) -> None:
        """Kinship >= 0.177 and < 0.354 is 1st degree."""
        assert get_degree(0.25) == 'RELATED_1ST'
        assert get_degree(0.177) == 'RELATED_1ST'

    def test_second_degree(self) -> None:
        """Kinship >= 0.0884 and < 0.177 is 2nd degree."""
        assert get_degree(0.125) == 'RELATED_2ND'
        assert get_degree(0.0884) == 'RELATED_2ND'

    def test_third_degree(self) -> None:
        """Kinship >= 0.0442 and < 0.0884 is 3rd degree."""
        assert get_degree(0.06) == 'RELATED_3RD'
        assert get_degree(0.0442) == 'RELATED_3RD'

    def test_below_threshold(self) -> None:
        """Kinship below 0.0442 returns None."""
        assert get_degree(0.01) is None


# -- process_kinship -----------------------------------------------------------


class TestProcessKinship:
    """Tests for process_kinship."""

    def test_basic_kinship(self, tmp_dir: Path) -> None:
        """Process a .kin0 file with one related pair."""
        path = _write(
            tmp_dir / 'test.kin0',
            '#IID1\tIID2\tKINSHIP\tIBS0\n'
            'S1\tS2\t0.25\t0.01\n',
        )
        result = process_kinship(path)
        assert isinstance(result, pd.DataFrame)
        assert set(DEGREE_COLS).issubset(result.columns)

        s1_row = result[result['IID'] == 'S1'].iloc[0]
        assert 'S2:0.25:0.01' in s1_row['RELATED_1ST']

        s2_row = result[result['IID'] == 'S2'].iloc[0]
        assert 'S1:0.25:0.01' in s2_row['RELATED_1ST']

    def test_filters_below_threshold(self, tmp_dir: Path) -> None:
        """Pairs below 3rd degree threshold are excluded."""
        path = _write(
            tmp_dir / 'low.kin0',
            '#IID1\tIID2\tKINSHIP\tIBS0\n'
            'S1\tS2\t0.01\t0.5\n',
        )
        result = process_kinship(path)
        assert result.empty

    def test_missing_columns(self, tmp_dir: Path) -> None:
        """Return empty DataFrame when required columns are missing."""
        path = _write(
            tmp_dir / 'bad.kin0',
            'COL_A\tCOL_B\nX\tY\n',
        )
        result = process_kinship(path)
        assert result.empty

    def test_empty_kin0(self, tmp_dir: Path) -> None:
        """Return empty DataFrame for a header-only .kin0 file."""
        path = _write(tmp_dir / 'empty.kin0', '#IID1\tIID2\tKINSHIP\tIBS0\n')
        result = process_kinship(path)
        assert result.empty

    def test_multiple_degrees(self, tmp_dir: Path) -> None:
        """A sample with relationships at different degrees gets separate columns."""
        path = _write(
            tmp_dir / 'multi.kin0',
            '#IID1\tIID2\tKINSHIP\tIBS0\n'
            'S1\tS2\t0.25\t0.01\n'
            'S1\tS3\t0.06\t0.1\n',
        )
        result = process_kinship(path)
        s1_row = result[result['IID'] == 'S1'].iloc[0]
        assert 'S2:0.25:0.01' in s1_row['RELATED_1ST']
        assert 'S3:0.06:0.1' in s1_row['RELATED_3RD']


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


# -- end-to-end ----------------------------------------------------------------


class TestEndToEnd:
    """End-to-end test calling main() with temp files."""

    def test_full_merge(self, tmp_dir: Path) -> None:
        """Run a full merge with all QC file types and verify CSV output."""
        smiss_path = _write(
            tmp_dir / 'test.smiss',
            '#FID\tIID\tMISS_PHENO_CT\tMISSING_CT\tOBS_CT\tF_MISS\n'
            'FAM1\tS1\t0\t10\t1000\t0.01\n'
            'FAM2\tS2\t0\t20\t1000\t0.02\n',
        )
        het_path = _write(
            tmp_dir / 'test.het',
            '#FID\tIID\tO_HOM\tE_HOM\tN_NM\tF\n'
            'FAM1\tS1\t500\t490\t1000\t0.02\n'
            'FAM2\tS2\t510\t490\t1000\t0.04\n',
        )
        sexcheck_path = _write(
            tmp_dir / 'test.sexcheck',
            '#FID\tIID\tPEDSEX\tSNPSEX\tSTATUS\tF\n'
            'FAM1\tS1\t1\t1\tOK\t0.99\n'
            'FAM2\tS2\t2\t2\tOK\t0.01\n',
        )
        kin0_path = _write(
            tmp_dir / 'test.kin0',
            '#IID1\tIID2\tKINSHIP\tIBS0\n'
            'S1\tS2\t0.25\t0.01\n',
        )
        baf_path = _write(
            tmp_dir / 'baf.txt',
            'sample_id\tLRR_mean\nS1\t0.005\nS2\t0.006\n',
        )
        output_path = str(tmp_dir / 'output.csv')

        # Call main with sys.argv override
        sys_argv = [
            'merge_qc.py',
            '--smiss', smiss_path,
            '--het', het_path,
            '--sexcheck', sexcheck_path,
            '--kin0', kin0_path,
            '--output', output_path,
            '--bafregress', baf_path,
        ]

        original_argv = sys.argv
        try:
            sys.argv = sys_argv
            main()
        finally:
            sys.argv = original_argv

        # Verify output
        with open(output_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == 2
        assert rows[0]['IID'] == 'S1'
        assert rows[1]['IID'] == 'S2'

        # Check kinship columns present
        for col in DEGREE_COLS:
            assert col in rows[0]

        # S1 should have a 1st degree relationship with S2
        assert 'S2' in rows[0]['RELATED_1ST']

        # Bafregress data should be merged
        assert 'LRR_mean' in rows[0]
        assert rows[0]['LRR_mean'] == '0.005'
