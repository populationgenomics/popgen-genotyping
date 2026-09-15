"""Tests for scripts/relatedness_exclusions.py."""

import csv
import sys
from pathlib import Path

import pandas as pd
import pytest

from popgen_genotyping.scripts.relatedness_exclusions import (
    main,
    read_bafregress,
    read_smiss,
    select_excluded_samples,
)


@pytest.fixture
def tmp_dir(tmp_path: Path) -> Path:
    """Return a temporary directory for test files."""
    return tmp_path


def _write(path: Path, text: str) -> str:
    """Write text to a file and return its string path."""
    path.write_text(text)
    return str(path)


# -- read_bafregress ------------------------------------------------------------


class TestReadBafregress:
    """Tests for read_bafregress."""

    def test_concatenates_multiple_files(self, tmp_dir: Path) -> None:
        """Multiple BAFRegress files are concatenated into one DataFrame."""
        path1 = _write(
            tmp_dir / 'plate1.txt',
            'sample_id\tbaf_regress\tNhom\nCPG1\t0.001\t1000\nCPG2\t0.03\t1000\n',
        )
        path2 = _write(
            tmp_dir / 'plate2.txt',
            'sample_id\tbaf_regress\tNhom\nCPG3\t0.12\t1000\n',
        )
        result = read_bafregress([path1, path2])
        assert list(result['sample_id']) == ['CPG1', 'CPG2', 'CPG3']
        assert list(result['baf_regress']) == [0.001, 0.03, 0.12]

    def test_missing_sample_id_column_raises(self, tmp_dir: Path) -> None:
        """A file lacking sample_id raises ValueError."""
        path = _write(tmp_dir / 'bad.txt', 'iid\tbaf_regress\nCPG1\t0.01\n')
        with pytest.raises(ValueError, match='missing required column'):
            read_bafregress([path])

    def test_missing_baf_regress_column_raises(self, tmp_dir: Path) -> None:
        """A file lacking baf_regress raises ValueError."""
        path = _write(tmp_dir / 'bad.txt', 'sample_id\tNhom\nCPG1\t1000\n')
        with pytest.raises(ValueError, match='missing required column'):
            read_bafregress([path])


# -- read_smiss ------------------------------------------------------------------


class TestReadSmiss:
    """Tests for read_smiss."""

    def test_strips_hash_and_reads_columns(self, tmp_dir: Path) -> None:
        """Leading '#' is stripped, and IID/F_MISS are present."""
        path = _write(
            tmp_dir / 'test.smiss',
            '#FID\tIID\tMISSING_CT\tOBS_CT\tF_MISS\nFAM1\tCPG1\t10\t1000\t0.01\n',
        )
        result = read_smiss(path)
        assert list(result.columns) == ['FID', 'IID', 'MISSING_CT', 'OBS_CT', 'F_MISS']
        assert result.iloc[0]['IID'] == 'CPG1'

    def test_missing_iid_column_raises(self, tmp_dir: Path) -> None:
        """A file lacking IID raises ValueError."""
        path = _write(tmp_dir / 'bad.smiss', '#FID\tF_MISS\nFAM1\t0.01\n')
        with pytest.raises(ValueError, match='missing required column'):
            read_smiss(path)

    def test_missing_fmiss_column_raises(self, tmp_dir: Path) -> None:
        """A file lacking F_MISS raises ValueError."""
        path = _write(tmp_dir / 'bad.smiss', '#FID\tIID\nFAM1\tCPG1\n')
        with pytest.raises(ValueError, match='missing required column'):
            read_smiss(path)


# -- select_excluded_samples ------------------------------------------------------


class TestSelectExcludedSamples:
    """Tests for select_excluded_samples."""

    def _bafregress(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                'sample_id': ['CPG1', 'CPG2', 'CPG3', 'CPG4'],
                'baf_regress': [0.001, 0.1125, 0.001, 0.001],
            },
        )

    def _smiss(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                'IID': ['CPG1', 'CPG2', 'CPG3', 'CPG4'],
                'F_MISS': [0.001, 0.001, 0.09, 0.001],
            },
        )

    def test_contamination_only(self) -> None:
        """A sample above contamination_max but below fmiss_max is flagged contamination."""
        result = select_excluded_samples(self._bafregress(), self._smiss(), contamination_max=0.02, fmiss_max=0.5)
        assert list(result['IID']) == ['CPG2']
        assert result.iloc[0]['REASON'] == 'contamination'

    def test_missingness_only(self) -> None:
        """A sample above fmiss_max but below contamination_max is flagged missingness."""
        result = select_excluded_samples(self._bafregress(), self._smiss(), contamination_max=0.5, fmiss_max=0.02)
        assert list(result['IID']) == ['CPG3']
        assert result.iloc[0]['REASON'] == 'missingness'

    def test_both_reasons_combined(self) -> None:
        """A sample exceeding both thresholds carries the combined reason string."""
        bafregress = pd.DataFrame({'sample_id': ['CPG1'], 'baf_regress': [0.5]})
        smiss = pd.DataFrame({'IID': ['CPG1'], 'F_MISS': [0.5]})
        result = select_excluded_samples(bafregress, smiss, contamination_max=0.02, fmiss_max=0.02)
        assert result.iloc[0]['REASON'] == 'contamination;missingness'

    def test_threshold_is_strictly_greater_than(self) -> None:
        """A sample exactly at either threshold is not excluded."""
        bafregress = pd.DataFrame({'sample_id': ['CPG1'], 'baf_regress': [0.02]})
        smiss = pd.DataFrame({'IID': ['CPG1'], 'F_MISS': [0.02]})
        result = select_excluded_samples(bafregress, smiss, contamination_max=0.02, fmiss_max=0.02)
        assert result.empty

    def test_sorted_by_iid(self) -> None:
        """Results are sorted by IID regardless of input order."""
        bafregress = pd.DataFrame({'sample_id': ['CPG9', 'CPG1'], 'baf_regress': [0.9, 0.8]})
        smiss = pd.DataFrame({'IID': ['CPG9', 'CPG1'], 'F_MISS': [0.001, 0.001]})
        result = select_excluded_samples(bafregress, smiss, contamination_max=0.1, fmiss_max=0.5)
        assert list(result['IID']) == ['CPG1', 'CPG9']

    def test_missing_bafregress_coverage_raises(self) -> None:
        """A .smiss IID with no BAFRegress estimate raises ValueError."""
        bafregress = pd.DataFrame({'sample_id': ['CPG1'], 'baf_regress': [0.001]})
        smiss = pd.DataFrame({'IID': ['CPG1', 'CPG2'], 'F_MISS': [0.001, 0.001]})
        with pytest.raises(ValueError, match='No BAFRegress estimate found'):
            select_excluded_samples(bafregress, smiss, contamination_max=0.02, fmiss_max=0.02)

    def test_empty_result_has_expected_columns(self) -> None:
        """No samples above either threshold returns an empty frame with the right columns."""
        result = select_excluded_samples(self._bafregress(), self._smiss(), contamination_max=0.99, fmiss_max=0.99)
        assert result.empty
        assert list(result.columns) == ['IID', 'REASON', 'BAF_REGRESS', 'F_MISS']


# -- end-to-end (main) ------------------------------------------------------------


class TestEndToEnd:
    """End-to-end test calling main() with temp files."""

    def test_writes_tsv_and_remove_list(self, tmp_dir: Path) -> None:
        """main() writes both the audit TSV and the plink2 #IID remove list."""
        baf_path = _write(
            tmp_dir / 'baf.txt',
            'sample_id\tbaf_regress\tNhom\nCPG1\t0.001\t1000\nCPG2\t0.1125\t1000\n',
        )
        smiss_path = _write(
            tmp_dir / 'test.smiss',
            '#FID\tIID\tMISSING_CT\tOBS_CT\tF_MISS\nFAM1\tCPG1\t1\t1000\t0.001\nFAM2\tCPG2\t1\t1000\t0.001\n',
        )
        output_tsv = str(tmp_dir / 'excluded.tsv')
        output_remove_list = str(tmp_dir / 'remove.txt')

        sys_argv = [
            'relatedness_exclusions.py',
            '--bafregress',
            baf_path,
            '--smiss',
            smiss_path,
            '--contamination-max',
            '0.02',
            '--fmiss-max',
            '0.02',
            '--output-tsv',
            output_tsv,
            '--output-remove-list',
            output_remove_list,
        ]

        original_argv = sys.argv
        try:
            sys.argv = sys_argv
            main()
        finally:
            sys.argv = original_argv

        with open(output_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)
        assert len(rows) == 1
        assert rows[0]['IID'] == 'CPG2'
        assert rows[0]['REASON'] == 'contamination'

        remove_lines = Path(output_remove_list).read_text().splitlines()
        assert remove_lines[0] == '#IID'
        assert remove_lines[1:] == ['CPG2']

    def test_header_only_outputs_when_nothing_excluded(self, tmp_dir: Path) -> None:
        """With nothing excluded, both outputs are written header-only."""
        baf_path = _write(
            tmp_dir / 'baf.txt',
            'sample_id\tbaf_regress\tNhom\nCPG1\t0.001\t1000\n',
        )
        smiss_path = _write(
            tmp_dir / 'test.smiss',
            '#FID\tIID\tMISSING_CT\tOBS_CT\tF_MISS\nFAM1\tCPG1\t1\t1000\t0.001\n',
        )
        output_tsv = str(tmp_dir / 'excluded.tsv')
        output_remove_list = str(tmp_dir / 'remove.txt')

        sys_argv = [
            'relatedness_exclusions.py',
            '--bafregress',
            baf_path,
            '--smiss',
            smiss_path,
            '--contamination-max',
            '0.02',
            '--fmiss-max',
            '0.02',
            '--output-tsv',
            output_tsv,
            '--output-remove-list',
            output_remove_list,
        ]

        original_argv = sys.argv
        try:
            sys.argv = sys_argv
            main()
        finally:
            sys.argv = original_argv

        with open(output_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)
        assert rows == []

        assert Path(output_remove_list).read_text() == '#IID\n'
