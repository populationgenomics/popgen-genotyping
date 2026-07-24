"""
Tests for utils.py.
"""

from unittest.mock import MagicMock, patch
import pytest
from popgen_genotyping.utils import get_output_prefix, parse_psam


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


@patch('popgen_genotyping.utils.config_retrieve')
@patch('popgen_genotyping.utils.get_workflow')
def test_get_output_prefix_versioned(mock_get_workflow, mock_config_retrieve, tmp_path):
    """Default behaviour appends the workflow.version segment."""
    mock_get_workflow.return_value.name = 'wf'
    mock_config_retrieve.return_value = 3  # workflow.version
    dataset = MagicMock()
    dataset.prefix.return_value = tmp_path

    result = get_output_prefix(dataset=dataset, stage_name='MyStage')

    dataset.prefix.assert_called_once_with(category='default')
    assert result == tmp_path / 'wf' / 'MyStage' / '3'


@patch('popgen_genotyping.utils.config_retrieve')
@patch('popgen_genotyping.utils.get_workflow')
def test_get_output_prefix_unversioned(mock_get_workflow, mock_config_retrieve, tmp_path):
    """versioned=False drops the version segment and never reads workflow.version."""
    mock_get_workflow.return_value.name = 'wf'
    dataset = MagicMock()
    dataset.prefix.return_value = tmp_path

    result = get_output_prefix(dataset=dataset, stage_name='CohortBcfToPlink', versioned=False)

    assert result == tmp_path / 'wf' / 'CohortBcfToPlink'
    mock_config_retrieve.assert_not_called()


@patch('popgen_genotyping.utils.config_retrieve')
@patch('popgen_genotyping.utils.get_workflow')
def test_get_output_prefix_tmp(mock_get_workflow, mock_config_retrieve, tmp_path):
    """tmp=True resolves the prefix from the 'tmp' category."""
    mock_get_workflow.return_value.name = 'wf'
    mock_config_retrieve.return_value = 1
    dataset = MagicMock()
    dataset.prefix.return_value = tmp_path

    result = get_output_prefix(dataset=dataset, stage_name='MyStage', tmp=True)

    dataset.prefix.assert_called_once_with(category='tmp')
    assert result == tmp_path / 'wf' / 'MyStage' / '1'
