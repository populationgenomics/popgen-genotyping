"""Tests for scripts/merge_qc.py."""

import csv
from pathlib import Path

import pytest

from popgen_genotyping.scripts.merge_qc import (
    get_degree,
    left_merge,
    process_bafregress,
    process_kinship,
    read_whitespace_file,
    write_csv,
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


# -- read_whitespace_file ------------------------------------------------------


class TestReadWhitespaceFile:
    """Tests for read_whitespace_file."""

    def test_basic_parsing(self, tmp_dir: Path) -> None:
        """Parse a simple whitespace-delimited file."""
        path = _write(
            tmp_dir / 'basic.txt',
            '#FID\tIID\tVAL\nFAM1\tS1\t0.5\nFAM2\tS2\t0.8\n',
        )
        rows = read_whitespace_file(path)
        assert len(rows) == 2
        assert rows[0] == {'FID': 'FAM1', 'IID': 'S1', 'VAL': '0.5'}
        assert rows[1] == {'FID': 'FAM2', 'IID': 'S2', 'VAL': '0.8'}

    def test_strips_hash_from_headers(self, tmp_dir: Path) -> None:
        """Verify that leading '#' is stripped from column headers."""
        path = _write(tmp_dir / 'hash.txt', '#IID\t#METRIC\nS1\t1.0\n')
        rows = read_whitespace_file(path)
        assert list(rows[0].keys()) == ['IID', 'METRIC']

    def test_empty_file(self, tmp_dir: Path) -> None:
        """Return an empty list for an empty file."""
        path = _write(tmp_dir / 'empty.txt', '')
        assert read_whitespace_file(path) == []

    def test_header_only(self, tmp_dir: Path) -> None:
        """Return an empty list for a file with only a header."""
        path = _write(tmp_dir / 'header.txt', 'COL1\tCOL2\n')
        assert read_whitespace_file(path) == []


# -- left_merge ----------------------------------------------------------------


class TestLeftMerge:
    """Tests for left_merge."""

    def test_basic_merge(self) -> None:
        """Merge two simple datasets on IID."""
        left = [{'IID': 'S1', 'A': '1'}, {'IID': 'S2', 'A': '2'}]
        right = [{'IID': 'S1', 'B': 'x'}, {'IID': 'S2', 'B': 'y'}]
        merged = left_merge(left, right, key='IID')
        assert merged[0] == {'IID': 'S1', 'A': '1', 'B': 'x'}
        assert merged[1] == {'IID': 'S2', 'A': '2', 'B': 'y'}

    def test_left_preserves_unmatched(self) -> None:
        """Rows in left without a match in right are preserved."""
        left = [{'IID': 'S1', 'A': '1'}, {'IID': 'S3', 'A': '3'}]
        right = [{'IID': 'S1', 'B': 'x'}]
        merged = left_merge(left, right, key='IID')
        assert len(merged) == 2
        assert 'B' not in merged[1]

    def test_suffix_on_duplicate_columns(self) -> None:
        """Duplicate column names get the specified suffix."""
        left = [{'IID': 'S1', 'VAL': '1'}]
        right = [{'IID': 'S1', 'VAL': '2'}]
        merged = left_merge(left, right, key='IID', suffix='_r')
        assert merged[0] == {'IID': 'S1', 'VAL': '1', 'VAL_r': '2'}


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
            '#IID1\tIID2\tKINSHIP\tIBS0\nS1\tS2\t0.25\t0.01\n',
        )
        result = process_kinship(path)
        assert 'S1' in result
        assert 'S2' in result
        assert 'S2:0.2500:0.0100' in result['S1']['RELATED_1ST']
        assert 'S1:0.2500:0.0100' in result['S2']['RELATED_1ST']

    def test_filters_below_threshold(self, tmp_dir: Path) -> None:
        """Pairs below 3rd degree threshold are excluded."""
        path = _write(
            tmp_dir / 'low.kin0',
            '#IID1\tIID2\tKINSHIP\tIBS0\nS1\tS2\t0.01\t0.5\n',
        )
        result = process_kinship(path)
        assert result == {}

    def test_missing_columns(self, tmp_dir: Path) -> None:
        """Return empty dict when required columns are missing."""
        path = _write(
            tmp_dir / 'bad.kin0',
            'COL_A\tCOL_B\nX\tY\n',
        )
        result = process_kinship(path)
        assert result == {}

    def test_empty_kin0(self, tmp_dir: Path) -> None:
        """Return empty dict for an empty .kin0 file."""
        path = _write(tmp_dir / 'empty.kin0', '')
        result = process_kinship(path)
        assert result == {}

    def test_multiple_degrees(self, tmp_dir: Path) -> None:
        """A sample with relationships at different degrees gets separate columns."""
        path = _write(
            tmp_dir / 'multi.kin0',
            '#IID1\tIID2\tKINSHIP\tIBS0\nS1\tS2\t0.25\t0.01\nS1\tS3\t0.06\t0.1\n',
        )
        result = process_kinship(path)
        assert 'S2:0.2500:0.0100' in result['S1']['RELATED_1ST']
        assert 'S3:0.0600:0.1000' in result['S1']['RELATED_3RD']


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
        rows = process_bafregress([path1, path2])
        assert len(rows) == 2
        assert rows[0] == {'IID': 'S1', 'LRR_mean': '0.01'}
        assert rows[1] == {'IID': 'S2', 'LRR_mean': '0.02'}

    def test_missing_file(self, tmp_dir: Path) -> None:
        """Missing files produce a warning but don't raise."""
        rows = process_bafregress([str(tmp_dir / 'nonexistent.txt')])
        assert rows == []

    def test_empty_list(self) -> None:
        """Empty input list returns empty output."""
        assert process_bafregress([]) == []


# -- write_csv -----------------------------------------------------------------


class TestWriteCsv:
    """Tests for write_csv."""

    def test_basic_write(self, tmp_dir: Path) -> None:
        """Write rows to CSV and verify contents."""
        out_path = str(tmp_dir / 'out.csv')
        rows = [
            {'IID': 'S1', 'VAL': '1'},
            {'IID': 'S2', 'VAL': '2'},
        ]
        write_csv(rows, out_path)

        with open(out_path) as f:
            reader = csv.DictReader(f)
            result = list(reader)
        assert len(result) == 2
        assert result[0] == {'IID': 'S1', 'VAL': '1'}

    def test_empty_rows(self, tmp_dir: Path, capsys: pytest.CaptureFixture[str]) -> None:
        """Empty rows list produces a warning and no file."""
        out_path = str(tmp_dir / 'empty.csv')
        write_csv([], out_path)
        captured = capsys.readouterr()
        assert 'Warning' in captured.out

    def test_extra_columns_in_later_rows(self, tmp_dir: Path) -> None:
        """Columns appearing only in later rows are included in the output."""
        out_path = str(tmp_dir / 'extra.csv')
        rows = [
            {'IID': 'S1', 'A': '1'},
            {'IID': 'S2', 'A': '2', 'B': '3'},
        ]
        write_csv(rows, out_path)

        with open(out_path) as f:
            reader = csv.DictReader(f)
            result = list(reader)
        assert 'B' in result[0]  # Column exists even for first row
        assert result[1]['B'] == '3'
