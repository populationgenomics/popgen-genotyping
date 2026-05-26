"""Tests for scripts/merge_qc.py."""

import csv
import sys
from pathlib import Path

import pandas as pd
import pytest

from popgen_genotyping.scripts.merge_qc import (
    DEGREE_COLS,
    INFTYPE_TO_DEGREE,
    main,
    process_bafregress,
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


_SEG_HEADER = 'FID1\tID1\tFID2\tID2\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType\n'


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


# -- INFTYPE_TO_DEGREE ---------------------------------------------------------


class TestInftypeToDegree:
    """Tests for the InfType → RELATED_* mapping."""

    def test_known_inftypes(self) -> None:
        """KING --ibdseg InfType values bin into the expected degree columns."""
        assert INFTYPE_TO_DEGREE['Dup/MZ'] == 'RELATED_MZ'
        assert INFTYPE_TO_DEGREE['PO'] == 'RELATED_1ST'
        assert INFTYPE_TO_DEGREE['FS'] == 'RELATED_1ST'
        assert INFTYPE_TO_DEGREE['2nd'] == 'RELATED_2ND'
        assert INFTYPE_TO_DEGREE['3rd'] == 'RELATED_3RD'

    def test_unknown_inftype(self) -> None:
        """Unknown InfType values (e.g. 4th, UN) are not in the mapping."""
        assert '4th' not in INFTYPE_TO_DEGREE
        assert 'UN' not in INFTYPE_TO_DEGREE


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
        assert set(DEGREE_COLS).issubset(result.columns)

        s1_row = result[result['IID'] == 'S1'].iloc[0]
        assert 'S2:0.25:FS' in s1_row['RELATED_1ST']

        s2_row = result[result['IID'] == 'S2'].iloc[0]
        assert 'S1:0.25:FS' in s2_row['RELATED_1ST']

    def test_parent_offspring(self, tmp_dir: Path) -> None:
        """PO pairs land in RELATED_1ST and carry InfType=PO."""
        # PO: IBD1Seg ≈ 1.0, IBD2Seg ≈ 0 → kinship = 0.25
        path = _write(
            tmp_dir / 'po.seg',
            _SEG_HEADER + 'F\tParent\tF\tChild\t1.0000\t0.0000\t0.5000\tPO\n',
        )
        result = process_seg(path)
        parent_row = result[result['IID'] == 'Parent'].iloc[0]
        assert 'Child:0.25:PO' in parent_row['RELATED_1ST']

    def test_unknown_inftype_raises(self, tmp_dir: Path) -> None:
        """Unrecognised InfType values raise ValueError."""
        path = _write(
            tmp_dir / 'unknown.seg',
            _SEG_HEADER + 'F\tS1\tF\tS2\t0.0100\t0.0000\t0.0050\tUN\n',
        )
        with pytest.raises(ValueError, match='unrecognised InfType'):
            process_seg(path)

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
        """A sample with relationships at different degrees gets separate columns."""
        path = _write(
            tmp_dir / 'multi.seg',
            _SEG_HEADER + 'F\tS1\tF\tS2\t0.5000\t0.2500\t0.7500\tFS\n' + 'F\tS1\tF\tS3\t0.1800\t0.0000\t0.0900\t3rd\n',
        )
        result = process_seg(path)
        s1_row = result[result['IID'] == 'S1'].iloc[0]
        assert 'S2:0.25:FS' in s1_row['RELATED_1ST']
        assert 'S3:0.045:3rd' in s1_row['RELATED_3RD']


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
            '#FID\tIID\tO_HOM\tE_HOM\tN_NM\tF\nFAM1\tS1\t500\t490\t1000\t0.02\nFAM2\tS2\t510\t490\t1000\t0.04\n',
        )
        sexcheck_path = _write(
            tmp_dir / 'test.sexcheck',
            '#FID\tIID\tPEDSEX\tSNPSEX\tSTATUS\tF\nFAM1\tS1\t1\t1\tOK\t0.99\nFAM2\tS2\t2\t2\tOK\t0.01\n',
        )
        seg_path = _write(
            tmp_dir / 'test.seg',
            _SEG_HEADER + 'FAM1\tS1\tFAM2\tS2\t0.5000\t0.2500\t0.7500\tFS\n',
        )
        baf_path = _write(
            tmp_dir / 'baf.txt',
            'sample_id\tLRR_mean\nS1\t0.005\nS2\t0.006\n',
        )
        output_path = str(tmp_dir / 'output.csv')

        # Call main with sys.argv override
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
            '--output',
            output_path,
            '--bafregress',
            baf_path,
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

        # Check degree columns present
        for col in DEGREE_COLS:
            assert col in rows[0]

        # S1 should have a 1st degree FS relationship with S2
        assert 'S2:0.25:FS' in rows[0]['RELATED_1ST']

        # Bafregress data should be merged
        assert 'LRR_mean' in rows[0]
        assert rows[0]['LRR_mean'] == '0.005'
