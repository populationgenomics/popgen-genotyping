"""
Tests for metamist_utils.py.
"""

from unittest.mock import patch, MagicMock
import pytest
from popgen_genotyping.metamist_utils import (
    query_genotyping_manifests,
    parse_genotyping_manifest,
    query_previous_aggregate,
    resolve_gtc_path
)


def generate_manifest(output_path, num_samples=10, prefix='CPGSYN'):
    """
    Generate a synthetic manifest CSV with all required columns.
    """
    import csv  # noqa: PLC0415
    headers = [
        'file_name', 'sample_sheet_id', 'sample_id', 'md5sum', 'bbv_barcode',
        'sentrix_barcode_a', 'sentrix_position_a', 'sample_plate', 'sample_well',
        'sample_plate_position', 'cpg_sample_id_internal', 'cpg_gcp_filepath',
        'cpg_sequencing_group_id', 'cpg_cohort_id'
    ]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        for i in range(1, num_samples + 1):
            sg_id = f'{prefix}{i:03d}'
            writer.writerow({
                'file_name': 'barcode_pos.gtc',
                'sample_sheet_id': '123',
                'sample_id': f'AGDB{i:05d}',
                'md5sum': 'md5',
                'bbv_barcode': '123',
                'sentrix_barcode_a': 'barcode',
                'sentrix_position_a': 'pos',
                'sample_plate': 'plate',
                'sample_well': 'well',
                'sample_plate_position': 'pos',
                'cpg_sample_id_internal': f'INT{i}',
                'cpg_gcp_filepath': f'gs://cpg-test-main/gtc/{sg_id}.gtc',
                'cpg_sequencing_group_id': sg_id,
                'cpg_cohort_id': 'COH1'
            })


@pytest.fixture
def synthetic_manifest(tmp_path):
    """
    Fixture to generate a synthetic manifest CSV for testing.
    """
    manifest_path = tmp_path / 'genotyping_array_manifest.csv'
    generate_manifest(manifest_path, num_samples=3, prefix='CPGSYN')
    return manifest_path


@patch('popgen_genotyping.metamist_utils.query')
@patch('popgen_genotyping.metamist_utils.config_retrieve')
def test_query_genotyping_manifests(mock_config, mock_query):
    # Mock config
    mock_config.return_value = 'ourdna'

    # Mock GraphQL response - id is now at the top level of the analysis object
    mock_query.return_value = {
        'project': {
            'analyses': [
                {
                    'id': '224611',
                    'type': 'manifest',
                    'outputs': {
                        'path': 'gs://cpg-ourdna-main/manifests/production_manifests/COH8495_production_manifest.csv',
                        'basename': 'COH8495_production_manifest.csv'
                    }
                },
                {
                    'id': '231239',
                    'type': 'manifest',
                    'outputs': {
                        'path': 'gs://cpg-ourdna-main/gtc_genotyping_array/manifests/genotyping_array_manifest_cohort_COH10152.csv',
                        'basename': 'genotyping_array_manifest_cohort_COH10152.csv'
                    }
                }
            ]
        }
    }

    # Execute
    manifests = query_genotyping_manifests('ourdna')

    # Verify filtering
    assert len(manifests) == 1
    assert manifests[0]['basename'] == 'genotyping_array_manifest_cohort_COH10152.csv'
    assert manifests[0]['id'] == '231239'


@patch('popgen_genotyping.metamist_utils.query_genotyping_manifests')
@patch('popgen_genotyping.metamist_utils.get_sequencing_group_cohort')
@patch('popgen_genotyping.metamist_utils.parse_genotyping_manifest')
def test_resolve_gtc_path(mock_parse, mock_get_cohort, mock_query):
    # Mock manifests with different IDs and basenames
    mock_query.return_value = [
        {'id': '10', 'basename': 'genotyping_array_COH001_old.csv', 'path': 'gs://path/v1.csv'},
        {'id': '20', 'basename': 'genotyping_array_COH001_latest.csv', 'path': 'gs://path/v2.csv'},
        {'id': '30', 'basename': 'genotyping_array_COH002.csv', 'path': 'gs://path/other.csv'},
    ]

    # Mock cohort
    mock_cohort = MagicMock()
    mock_cohort.id = 'COH001'
    mock_get_cohort.return_value = mock_cohort

    # Mock sequencing group
    mock_sg = MagicMock()
    mock_sg.id = 'CPG001'
    mock_sg.dataset.name = 'test_dataset'

    # Mock manifest parsing
    mock_parse.return_value = {'CPG001': 'gs://bucket/CPG001.gtc'}

    # Execute
    path = resolve_gtc_path(mock_sg)

    # Verify
    mock_parse.assert_called_with('gs://path/v2.csv')
    assert path == 'gs://bucket/CPG001.gtc'


@patch('popgen_genotyping.metamist_utils.query_genotyping_manifests')
@patch('popgen_genotyping.metamist_utils.get_sequencing_group_cohort')
def test_resolve_gtc_path_no_manifest(mock_get_cohort, mock_query):
    mock_query.return_value = [
        {'id': '10', 'basename': 'some_other_manifest.csv', 'path': 'gs://path/v1.csv'},
    ]
    mock_cohort = MagicMock()
    mock_cohort.id = 'COH001'
    mock_get_cohort.return_value = mock_cohort
    mock_sg = MagicMock()
    mock_sg.id = 'CPG001'
    mock_sg.dataset.name = 'test_dataset'

    with pytest.raises(ValueError, match='No manifest found for cohort COH001'):
        resolve_gtc_path(mock_sg)


@patch('popgen_genotyping.metamist_utils.to_path')
def test_parse_genotyping_manifest(mock_to_path, synthetic_manifest):
    # Mock to_path to return the path to the synthetic manifest
    mock_to_path.return_value = synthetic_manifest

    # Execute
    mapping = parse_genotyping_manifest(str(synthetic_manifest))

    # Verify
    assert mapping == {
        'CPGSYN001': 'gs://cpg-test-main/gtc/CPGSYN001.gtc',
        'CPGSYN002': 'gs://cpg-test-main/gtc/CPGSYN002.gtc',
        'CPGSYN003': 'gs://cpg-test-main/gtc/CPGSYN003.gtc'
    }


@patch('popgen_genotyping.metamist_utils.query')
@patch('popgen_genotyping.metamist_utils.config_retrieve')
def test_query_genotyping_manifests_no_results(mock_config, mock_query):
    mock_config.return_value = 'ourdna'
    mock_query.return_value = {'project': {'analyses': []}}

    manifests = query_genotyping_manifests('ourdna')
    assert manifests == []


@patch('popgen_genotyping.metamist_utils.query')
def test_query_previous_aggregate(mock_query):
    mock_query.return_value = {
        'analyses': [
            {
                'outputs': {'pgen': 'gs://path/merged.pgen'},
                'project': {
                    'sequencingGroups': [
                        {'id': 'CPG001'},
                        {'id': 'CPG002'}
                    ]
                }
            }
        ]
    }

    outputs, active_sgs = query_previous_aggregate(123)

    assert outputs == {'pgen': 'gs://path/merged.pgen'}
    assert active_sgs == ['CPG001', 'CPG002']
    mock_query.assert_called_once()
