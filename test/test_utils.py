"""
Tests for utils.py.
"""

from pathlib import Path
from unittest.mock import patch
import pytest
from popgen_genotyping.utils import parse_psam, reconstruct_cohort_output_paths


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


@patch('popgen_genotyping.utils.get_output_prefix')
def test_reconstruct_cohort_output_paths_plink(mock_prefix):
    # Simulate the phase-1 CohortBcfToPlink tmp prefix.
    mock_prefix.return_value = Path('/cpg-ds-main-tmp/popgen_genotyping/CohortBcfToPlink/1')
    dataset = object()

    result = reconstruct_cohort_output_paths(
        dataset=dataset,
        stage_name='CohortBcfToPlink',
        cohort_ids=['COH1', 'COH2'],
        suffixes={'bed': '.bed', 'bim': '.bim', 'fam': '.fam'},
        tmp=True,
    )

    mock_prefix.assert_called_once_with(dataset=dataset, stage_name='CohortBcfToPlink', tmp=True)
    assert result == [
        {
            'bed': '/cpg-ds-main-tmp/popgen_genotyping/CohortBcfToPlink/1/COH1.bed',
            'bim': '/cpg-ds-main-tmp/popgen_genotyping/CohortBcfToPlink/1/COH1.bim',
            'fam': '/cpg-ds-main-tmp/popgen_genotyping/CohortBcfToPlink/1/COH1.fam',
        },
        {
            'bed': '/cpg-ds-main-tmp/popgen_genotyping/CohortBcfToPlink/1/COH2.bed',
            'bim': '/cpg-ds-main-tmp/popgen_genotyping/CohortBcfToPlink/1/COH2.bim',
            'fam': '/cpg-ds-main-tmp/popgen_genotyping/CohortBcfToPlink/1/COH2.fam',
        },
    ]


@patch('popgen_genotyping.utils.get_output_prefix')
def test_reconstruct_cohort_output_paths_bafregress(mock_prefix):
    # BafRegress writes to retained (non-tmp) storage.
    mock_prefix.return_value = Path('/cpg-ds-main/popgen_genotyping/BafRegress/1')
    dataset = object()

    result = reconstruct_cohort_output_paths(
        dataset=dataset,
        stage_name='BafRegress',
        cohort_ids=['COH1'],
        suffixes={'txt': '.BAFRegress.txt'},
        tmp=False,
    )

    mock_prefix.assert_called_once_with(dataset=dataset, stage_name='BafRegress', tmp=False)
    assert result == [{'txt': '/cpg-ds-main/popgen_genotyping/BafRegress/1/COH1.BAFRegress.txt'}]


@patch('popgen_genotyping.utils.get_output_prefix')
def test_reconstruct_cohort_output_paths_empty(mock_prefix):
    mock_prefix.return_value = Path('/cpg-ds-main-tmp/popgen_genotyping/CohortBcfToPlink/1')

    result = reconstruct_cohort_output_paths(
        dataset=object(),
        stage_name='CohortBcfToPlink',
        cohort_ids=[],
        suffixes={'bed': '.bed'},
        tmp=True,
    )

    assert result == []
