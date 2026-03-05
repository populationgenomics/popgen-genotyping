"""
Tests for metamist_utils module.
"""

from unittest.mock import patch, MagicMock
from popgen_genotyping.metamist_utils import query_reported_sex


def test_query_reported_sex_parsing():
    """
    Test that query_reported_sex correctly parses the nested Metamist response.
    """
    mock_response = {
        'project': {
            'sequencingGroups': [
                {
                    'id': 'CPG001',
                    'sample': {
                        'participant': {
                            'reportedSex': 'Female'
                        }
                    }
                },
                {
                    'id': 'CPG002',
                    'sample': {
                        'participant': {
                            'reportedSex': 'Male'
                        }
                    }
                },
                {
                    'id': 'CPG003',
                    'sample': None # Test handling of missing sample
                },
                {
                    'id': 'CPG004',
                    'sample': {
                        'participant': None # Test handling of missing participant
                    }
                }
            ]
        }
    }

    with patch('popgen_genotyping.metamist_utils.query', return_value=mock_response), \
         patch('popgen_genotyping.metamist_utils.config_retrieve', return_value='ourdna'):

        result = query_reported_sex()

        assert result == {
            'CPG001': 'Female',
            'CPG002': 'Male'
        }
        assert len(result) == 2


def test_query_reported_sex_empty_response():
    """
    Test that query_reported_sex returns an empty dict on empty project or sequencingGroups.
    """
    mock_responses = [
        {'project': None},
        {'project': {'sequencingGroups': []}}
    ]

    for mock_resp in mock_responses:
        with patch('popgen_genotyping.metamist_utils.query', return_value=mock_resp), \
             patch('popgen_genotyping.metamist_utils.config_retrieve', return_value='ourdna'):

            result = query_reported_sex()
            assert result == {}
