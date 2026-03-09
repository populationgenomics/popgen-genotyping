"""
Tests for metamist_utils.py.
"""

from unittest.mock import patch
import pytest
from popgen_genotyping.metamist_utils import (
    query_genotyping_manifests,
    parse_genotyping_manifest,
    query_previous_aggregate
)
from popgen_genotyping.scripts.generate_synthetic_manifest import generate_manifest


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

    # Mock GraphQL response
    mock_query.return_value = {
        'project': {
            'analyses': [
                {
                    'type': 'manifest',
                    'outputs': {
                        'id': 224611,
                        'path': 'gs://cpg-ourdna-main/manifests/production_manifests/COH8495_production_manifest.csv',
                        'basename': 'COH8495_production_manifest.csv'
                    }
                },
                {
                    'type': 'manifest',
                    'outputs': {
                        'id': 231239,
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
