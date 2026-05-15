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
    process_ibdseg,
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


# -- INFTYPE_TO_DEGREE ---------------------------------------------------------


class TestInfTypeMap:
    """Sanity checks on the InfType → degree-column map."""

    def test_known_inftypes(self) -> None:
        """All canonical KING InfType labels map to a degree column."""
        assert INFTYPE_TO_DEGREE['Dup/MZ'] == 'RELATED_MZ'
        assert INFTYPE_TO_DEGREE['PO'] == 'RELATED_1ST'
        assert INFTYPE_TO_DEGREE['FS'] == 'RELATED_1ST'
        assert INFTYPE_TO_DEGREE['2nd'] == 'RELATED_2ND'
        assert INFTYPE_TO_DEGREE['3rd'] == 'RELATED_3RD'

    def test_map_only_covers_known(self) -> None:
        """The map exposes exactly the inference types the report tracks."""
        assert set(INFTYPE_TO_DEGREE.values()) == set(DEGREE_COLS)


# -- process_ibdseg ------------------------------------------------------------


def _seg_header() -> str:
    """Header line for a KING --ibdseg .seg file."""
    return 'FID1\tID1\tFID2\tID2\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType'


class TestProcessIbdseg:
    """Tests for process_ibdseg."""

    def test_basic_ibdseg(self, tmp_dir: Path) -> None:
        """Process a .seg file with one MZ, one PO and one 2nd-degree pair."""
        path = _write(
            tmp_dir / 'test.seg',
            _seg_header() + '\n'
            # Dup/MZ: full identity
            'F1\tS1\tF2\tS2\t0.0\t1.0\t1.0\tDup/MZ\n'
            # PO: ~half IBD1, no IBD2
            'F3\tS3\tF4\tS4\t1.0\t0.0\t0.5\tPO\n'
            # 2nd degree
            'F5\tS5\tF6\tS6\t0.5\t0.0\t0.25\t2nd\n',
        )
        result = process_ibdseg(path)
        assert isinstance(result, pd.DataFrame)
        assert set(DEGREE_COLS).issubset(result.columns)

        # S1 should carry the MZ payload referencing S2.
        s1_row = result[result['IID'] == 'S1'].iloc[0]
        assert 'S2:0.0:0.0:1.0:1.0' in s1_row['RELATED_MZ']

        # S2 mirrors S1.
        s2_row = result[result['IID'] == 'S2'].iloc[0]
        assert 'S1:0.0:0.0:1.0:1.0' in s2_row['RELATED_MZ']

        # PO pair lands in 1st-degree column with IBD0 ≈ 0.
        s3_row = result[result['IID'] == 'S3'].iloc[0]
        assert 'S4:0.0:1.0:0.0:0.5' in s3_row['RELATED_1ST']

        # 2nd-degree pair lands in 2nd column with IBD0 = 0.5.
        s5_row = result[result['IID'] == 'S5'].iloc[0]
        assert 'S6:0.5:0.5:0.0:0.25' in s5_row['RELATED_2ND']

    def test_unknown_inftype_excluded(self, tmp_dir: Path) -> None:
        """Pairs with unmapped InfType (e.g. '4th', 'UN') are dropped."""
        path = _write(
            tmp_dir / 'unmapped.seg',
            _seg_header() + '\nF1\tS1\tF2\tS2\t0.1\t0.0\t0.05\t4th\nF3\tS3\tF4\tS4\t0.0\t0.0\t0.0\tUN\n',
        )
        result = process_ibdseg(path)
        assert result.empty

    def test_missing_columns(self, tmp_dir: Path) -> None:
        """Return empty DataFrame when required columns are missing."""
        path = _write(
            tmp_dir / 'bad.seg',
            'COL_A\tCOL_B\nX\tY\n',
        )
        result = process_ibdseg(path)
        assert result.empty

    def test_empty_seg(self, tmp_dir: Path) -> None:
        """Return empty DataFrame for a header-only .seg file."""
        path = _write(tmp_dir / 'empty.seg', _seg_header() + '\n')
        result = process_ibdseg(path)
        assert result.empty

    def test_multiple_degrees(self, tmp_dir: Path) -> None:
        """A sample with relationships at different degrees gets separate columns."""
        path = _write(
            tmp_dir / 'multi.seg',
            _seg_header() + '\nF1\tS1\tF2\tS2\t1.0\t0.0\t0.5\tPO\nF1\tS1\tF3\tS3\t0.2\t0.0\t0.1\t3rd\n',
        )
        result = process_ibdseg(path)
        s1_row = result[result['IID'] == 'S1'].iloc[0]
        assert 'S2:0.0:1.0:0.0:0.5' in s1_row['RELATED_1ST']
        assert 'S3:0.8:0.2:0.0:0.1' in s1_row['RELATED_3RD']

    def test_ibd0_clipped_to_zero(self, tmp_dir: Path) -> None:
        """IBD0 derivation clamps at 0 when IBD1+IBD2 > 1 due to noise."""
        path = _write(
            tmp_dir / 'clip.seg',
            _seg_header() + '\n'
            # IBD1+IBD2 = 1.01, would give IBD0 = -0.01 without clipping
            'F1\tS1\tF2\tS2\t0.51\t0.5\t0.755\tDup/MZ\n',
        )
        result = process_ibdseg(path)
        s1_row = result[result['IID'] == 'S1'].iloc[0]
        assert ':0.0:' in s1_row['RELATED_MZ']


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
            _seg_header() + '\nFAM1\tS1\tFAM2\tS2\t1.0\t0.0\t0.5\tPO\n',
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
            '--ibdseg',
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

        # Check relatedness columns present
        for col in DEGREE_COLS:
            assert col in rows[0]

        # S1 should have a 1st-degree (PO) relationship with S2, encoded
        # as REL_ID:IBD0:IBD1:IBD2:PropIBD.
        assert 'S2:0.0:1.0:0.0:0.5' in rows[0]['RELATED_1ST']

        # Bafregress data should be merged
        assert 'LRR_mean' in rows[0]
        assert rows[0]['LRR_mean'] == '0.005'
