"""
Tests for utils.py.
"""

from unittest.mock import patch
import pytest
from popgen_genotyping.utils import parse_psam


@pytest.fixture
def synthetic_psam(tmp_path):
    """
    Fixture to generate a synthetic PSAM file for testing.
    """
    psam_path = tmp_path / 'test.psam'
    content = '#IID\tFID\tSEX\nCPG001\tFAM1\t1\nCPG002\tFAM2\t2\n'
    psam_path.write_text(content)
    return psam_path


@patch('popgen_genotyping.utils.to_path')
def test_parse_psam(mock_to_path, synthetic_psam):
    # Mock to_path to return the local synthetic psam path
    mock_to_path.return_value = synthetic_psam

    ids = parse_psam(str(synthetic_psam))

    assert ids == ['CPG001', 'CPG002']


@patch('popgen_genotyping.utils.to_path')
def test_parse_psam_no_header(mock_to_path, tmp_path):
    psam_path = tmp_path / 'no_header.psam'
    content = 'CPG123\tOTHER\nCPG456\tOTHER\n'
    psam_path.write_text(content)

    mock_to_path.return_value = psam_path

    ids = parse_psam(str(psam_path))
    assert ids == ['CPG123', 'CPG456']
