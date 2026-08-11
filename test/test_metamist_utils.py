"""
Tests for metamist_utils.py.
"""

from typing import Any
from unittest.mock import MagicMock, patch

import pytest
from metamist.model.new_cohort import NewCohort

from popgen_genotyping.metamist_utils import (
    create_custom_cohort,
    find_cohort_by_membership,
    format_merge_plan,
    parse_genotyping_manifest,
    query_cohorts_with_analyses,
    query_genotyping_manifests,
    resolve_bafregress_map,
    resolve_cohort_bed_map,
    resolve_gtc_path,
    resolve_merge_inputs,
    resolve_super_cohort_membership,
    wait_for_cohort_analyses,
)


def _analysis(analysis_type, path, ts='2026-01-01T00:00:00', analysis_id=1) -> dict[str, Any]:
    """Build a synthetic analysis record for cohort-first query mocks."""
    return {'id': analysis_id, 'type': analysis_type, 'timestampCompleted': ts, 'outputs': {'path': path}}


def _cohort(cohort_id, sg_ids, analyses=None) -> dict[str, Any]:
    """Build a synthetic cohort record as returned by query_cohorts_with_analyses."""
    return {
        'id': cohort_id,
        'name': cohort_id,
        'sequencingGroups': [{'id': s} for s in sg_ids],
        'analyses': analyses or [],
    }


def generate_manifest(output_path, num_samples=10, prefix='CPGSYN'):
    """
    Generate a synthetic manifest CSV with all required columns.
    """
    import csv  # noqa: PLC0415

    headers = [
        'file_name',
        'sample_sheet_id',
        'sample_id',
        'md5sum',
        'bbv_barcode',
        'sentrix_barcode_a',
        'sentrix_position_a',
        'sample_plate',
        'sample_well',
        'sample_plate_position',
        'cpg_sample_id_internal',
        'cpg_gcp_filepath',
        'cpg_sequencing_group_id',
        'cpg_cohort_id',
    ]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        for i in range(1, num_samples + 1):
            sg_id = f'{prefix}{i:03d}'
            writer.writerow(
                {
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
                    'cpg_cohort_id': 'COH1',
                }
            )


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
                        'basename': 'COH8495_production_manifest.csv',
                    },
                },
                {
                    'id': '231239',
                    'type': 'manifest',
                    'outputs': {
                        'path': 'gs://cpg-ourdna-main/gtc_genotyping_array/manifests/genotyping_array_manifest_cohort_COH10152.csv',
                        'basename': 'genotyping_array_manifest_cohort_COH10152.csv',
                    },
                },
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
        'CPGSYN003': 'gs://cpg-test-main/gtc/CPGSYN003.gtc',
    }


@patch('popgen_genotyping.metamist_utils.query')
@patch('popgen_genotyping.metamist_utils.config_retrieve')
def test_query_genotyping_manifests_no_results(mock_config, mock_query):
    mock_config.return_value = 'ourdna'
    mock_query.return_value = {'project': {'analyses': []}}

    manifests = query_genotyping_manifests('ourdna')
    assert manifests == []


# ---------------------------------------------------------------------------
# query_cohorts_with_analyses
# ---------------------------------------------------------------------------


@patch('popgen_genotyping.metamist_utils.query')
@patch('popgen_genotyping.metamist_utils.config_retrieve')
def test_query_cohorts_with_analyses(mock_config, mock_query):
    mock_config.return_value = 'full'  # access_level (not 'test', so no suffix)
    mock_query.return_value = {'project': {'cohorts': [_cohort('COH5', ['CPG1'])]}}

    cohorts = query_cohorts_with_analyses('ourdna')

    assert len(cohorts) == 1
    assert cohorts[0]['id'] == 'COH5'


@patch('popgen_genotyping.metamist_utils.query')
@patch('popgen_genotyping.metamist_utils.config_retrieve')
def test_query_cohorts_with_analyses_no_project(mock_config, mock_query):
    mock_config.return_value = 'full'
    mock_query.return_value = {'project': None}

    assert query_cohorts_with_analyses('ourdna') == []


# ---------------------------------------------------------------------------
# resolve_cohort_bed_map
# ---------------------------------------------------------------------------


def test_resolve_cohort_bed_map():
    cohorts = [
        _cohort('COH5', ['CPG1', 'CPG2'], [_analysis('array_cohort_bed', 'gs://b/COH5.bed')]),
        _cohort('COH6', ['CPG3'], [_analysis('array_cohort_bed', 'gs://b/COH6.bed')]),
        # Aggregate cohort shares SGs but carries no array_cohort_bed → ignored, no clash.
        _cohort('COH900', ['CPG1', 'CPG3'], [_analysis('array_aggregate_pgen', 'gs://b/agg.pgen')]),
    ]

    result = resolve_cohort_bed_map(cohorts=cohorts)

    assert set(result) == {'CPG1', 'CPG2', 'CPG3'}
    assert result['CPG1'] == {
        'bed': 'gs://b/COH5.bed',
        'bim': 'gs://b/COH5.bim',
        'fam': 'gs://b/COH5.fam',
        'cohort_id': 'COH5',
    }
    assert result['CPG3']['cohort_id'] == 'COH6'


def test_resolve_cohort_bed_map_latest_timestamp_wins():
    cohorts = [
        _cohort(
            'COH5',
            ['CPG1'],
            [
                _analysis('array_cohort_bed', 'gs://b/old.bed', ts='2026-01-01T00:00:00', analysis_id=1),
                _analysis('array_cohort_bed', 'gs://b/new.bed', ts='2026-06-01T00:00:00', analysis_id=2),
            ],
        ),
    ]

    result = resolve_cohort_bed_map(cohorts=cohorts)

    assert result['CPG1']['bed'] == 'gs://b/new.bed'


def test_resolve_cohort_bed_map_duplicate_sg_raises():
    cohorts = [
        _cohort('COH5', ['CPG1'], [_analysis('array_cohort_bed', 'gs://b/COH5.bed')]),
        _cohort('COH6', ['CPG1'], [_analysis('array_cohort_bed', 'gs://b/COH6.bed')]),
    ]

    with pytest.raises(ValueError, match='maps to multiple array_cohort_bed cohorts'):
        resolve_cohort_bed_map(cohorts=cohorts)


def test_resolve_cohort_bed_map_non_bed_path_raises():
    cohorts = [_cohort('COH5', ['CPG1'], [_analysis('array_cohort_bed', 'gs://b/COH5.pgen')])]

    with pytest.raises(ValueError, match=r'is not a \.bed path'):
        resolve_cohort_bed_map(cohorts=cohorts)


# ---------------------------------------------------------------------------
# resolve_bafregress_map
# ---------------------------------------------------------------------------


def test_resolve_bafregress_map():
    cohorts = [
        _cohort('COH5', ['CPG1', 'CPG2'], [_analysis('array_bafregress', 'gs://b/COH5.BAFRegress.txt')]),
        _cohort('COH6', ['CPG3'], [_analysis('array_bafregress', 'gs://b/COH6.BAFRegress.txt')]),
    ]

    result = resolve_bafregress_map(['CPG1', 'CPG3'], cohorts=cohorts)

    assert result == {'CPG1': 'gs://b/COH5.BAFRegress.txt', 'CPG3': 'gs://b/COH6.BAFRegress.txt'}


def test_resolve_bafregress_map_missing_raises():
    cohorts = [_cohort('COH5', ['CPG1'], [_analysis('array_bafregress', 'gs://b/COH5.BAFRegress.txt')])]

    with pytest.raises(ValueError, match='No BafRegress analysis found'):
        resolve_bafregress_map(['CPG1', 'CPG_missing'], cohorts=cohorts)


# ---------------------------------------------------------------------------
# resolve_merge_inputs
# ---------------------------------------------------------------------------


@patch('popgen_genotyping.metamist_utils.query_cohorts_with_analyses')
def test_resolve_merge_inputs_rolling(mock_cohorts):
    # COH_old was already aggregated into COH900 (its SGs are in AGG) → not re-merged
    # (aggregate-priority). COH5/COH6 are new plates, disjoint from AGG. CPG3 was withdrawn.
    mock_cohorts.return_value = [
        _cohort('COH900', ['CPG1', 'CPG2', 'CPG3'], [_analysis('array_aggregate_pgen', 'gs://b/agg.pgen')]),
        _cohort('COH_old', ['CPG1', 'CPG2'], [_analysis('array_cohort_bed', 'gs://b/COH_old.bed')]),
        _cohort('COH5', ['CPG4', 'CPG5'], [_analysis('array_cohort_bed', 'gs://b/COH5.bed')]),
        _cohort('COH6', ['CPG6'], [_analysis('array_cohort_bed', 'gs://b/COH6.bed')]),
    ]

    result = resolve_merge_inputs(['CPG1', 'CPG2', 'CPG4', 'CPG5', 'CPG6'], 'COH900')

    assert result['previous_aggregate_paths'] == {
        'pgen': 'gs://b/agg.pgen',
        'pvar': 'gs://b/agg.pvar',
        'psam': 'gs://b/agg.psam',
    }
    assert result['samples_to_remove'] == ['CPG3']  # in aggregate, dropped from super cohort
    assert result['super_cohort_size'] == 5

    plates = {p['cohort_id']: p for p in result['plate_merge_list']}
    assert set(plates) == {'COH5', 'COH6'}  # COH_old NOT pulled — already aggregated
    assert plates['COH5']['new_count'] == 2
    assert plates['COH6']['new_count'] == 1
    assert plates['COH5']['bim'] == 'gs://b/COH5.bim'


@patch('popgen_genotyping.metamist_utils.query_cohorts_with_analyses')
def test_resolve_merge_inputs_bootstrap(mock_cohorts):
    mock_cohorts.return_value = [
        _cohort('COH5', ['CPG1', 'CPG2'], [_analysis('array_cohort_bed', 'gs://b/COH5.bed')]),
    ]

    result = resolve_merge_inputs(['CPG1', 'CPG2'], None)

    assert result['previous_aggregate_paths'] is None
    assert result['samples_to_remove'] == []
    assert len(result['plate_merge_list']) == 1
    assert result['plate_merge_list'][0]['new_count'] == 2


@patch('popgen_genotyping.metamist_utils.query_cohorts_with_analyses')
def test_resolve_merge_inputs_missing_prev_cohort_raises(mock_cohorts):
    mock_cohorts.return_value = [_cohort('COH5', ['CPG1'], [_analysis('array_cohort_bed', 'gs://b/COH5.bed')])]

    with pytest.raises(ValueError, match='Previous aggregate cohort COH999 not found'):
        resolve_merge_inputs(['CPG1'], 'COH999')


@patch('popgen_genotyping.metamist_utils.query_cohorts_with_analyses')
def test_resolve_merge_inputs_prev_no_aggregate_analysis_raises(mock_cohorts):
    mock_cohorts.return_value = [_cohort('COH900', ['CPG1'], [])]

    with pytest.raises(ValueError, match='no array_aggregate_pgen analysis'):
        resolve_merge_inputs(['CPG1'], 'COH900')


@patch('popgen_genotyping.metamist_utils.query_cohorts_with_analyses')
def test_resolve_merge_inputs_missing_bed_for_new_sg_raises(mock_cohorts):
    mock_cohorts.return_value = [
        _cohort('COH900', ['CPG1'], [_analysis('array_aggregate_pgen', 'gs://b/agg.pgen')]),
    ]

    # CPG2 is new (in super cohort, not in aggregate) but has no registered bed anywhere.
    with pytest.raises(ValueError, match='No array_cohort_bed analysis found for new sequencing groups'):
        resolve_merge_inputs(['CPG1', 'CPG2'], 'COH900')


# ---------------------------------------------------------------------------
# format_merge_plan
# ---------------------------------------------------------------------------


def test_format_merge_plan():
    resolved = {
        'previous_aggregate_paths': {'pgen': 'x'},
        'samples_to_remove': ['CPG3'],
        'plate_merge_list': [
            {'cohort_id': 'COH6', 'new_count': 1},
            {'cohort_id': 'COH5', 'new_count': 2},
        ],
        'super_cohort_size': 5,
    }

    text = format_merge_plan(resolved, 'COH900')

    assert 'previous aggregate cohort COH900' in text
    assert 'super cohort:          5 SGs' in text
    assert 'carried forward:       2 SGs' in text
    assert 'withdrawn:             1 SGs' in text
    assert 'CPG3' in text
    assert 'COH5:  2 new SGs' in text
    assert 'COH6:  1 new SGs' in text
    assert 'expected merged total: 5 SGs' in text
    # Plates are listed sorted by cohort ID.
    assert text.index('COH5:') < text.index('COH6:')


def test_format_merge_plan_bootstrap():
    resolved = {
        'previous_aggregate_paths': None,
        'samples_to_remove': [],
        'plate_merge_list': [{'cohort_id': 'COH5', 'new_count': 2}],
        'super_cohort_size': 2,
    }

    text = format_merge_plan(resolved, None)

    assert 'bootstrap' in text
    assert 'carried forward:       0 SGs' in text
    assert 'expected merged total: 2 SGs' in text


def test_resolve_super_cohort_membership_rolling():
    cohorts = [_cohort('COH_AGG', ['CPG1', 'CPG2', 'CPG3'])]

    membership = resolve_super_cohort_membership(
        plate_sg_ids=['CPG4', 'CPG5'],
        previous_aggregate_cohort_id='COH_AGG',
        cohorts=cohorts,
    )

    assert membership == ['CPG1', 'CPG2', 'CPG3', 'CPG4', 'CPG5']


def test_resolve_super_cohort_membership_bootstrap():
    membership = resolve_super_cohort_membership(
        plate_sg_ids=['CPG2', 'CPG1'],
        previous_aggregate_cohort_id=None,
    )

    assert membership == ['CPG1', 'CPG2']


def test_resolve_super_cohort_membership_missing_aggregate_raises():
    with pytest.raises(ValueError, match='COH_MISSING'):
        resolve_super_cohort_membership(
            plate_sg_ids=['CPG1'],
            previous_aggregate_cohort_id='COH_MISSING',
            cohorts=[_cohort('COH_OTHER', ['CPG9'])],
        )


def test_find_cohort_by_membership():
    cohorts = [
        _cohort('COH1', ['CPG1', 'CPG2']),
        _cohort('COH2', ['CPG1', 'CPG2', 'CPG3']),
    ]

    match = find_cohort_by_membership(['CPG3', 'CPG2', 'CPG1'], cohorts=cohorts)

    assert match is not None
    assert match['id'] == 'COH2'


def test_find_cohort_by_membership_no_match():
    cohorts = [_cohort('COH1', ['CPG1', 'CPG2'])]

    assert find_cohort_by_membership(['CPG1'], cohorts=cohorts) is None


def test_find_cohort_by_membership_latest_wins():
    # COH99 vs COH100 crosses a digit-length boundary, where lexicographic
    # comparison would wrongly pick COH99.
    cohorts = [
        _cohort('COH100', ['CPG1', 'CPG2']),
        _cohort('COH99', ['CPG1', 'CPG2']),
        _cohort('COH9', ['CPG1', 'CPG2']),
    ]

    match = find_cohort_by_membership(['CPG1', 'CPG2'], cohorts=cohorts)

    assert match is not None
    assert match['id'] == 'COH100'


@patch('popgen_genotyping.metamist_utils.CohortApi')
@patch('popgen_genotyping.metamist_utils.metamist_project')
def test_create_custom_cohort(mock_project, mock_api):
    mock_project.return_value = 'ourdna'
    mock_api.return_value.create_cohort_from_criteria.return_value = NewCohort(
        cohort_id='COH42',
        sequencing_group_ids=['CPG1', 'CPG2'],
        dry_run=False,
    )

    cohort_id = create_custom_cohort(
        name='array-aggregate-2026-08-10',
        description='test cohort',
        sg_ids=['CPG2', 'CPG1'],
    )

    assert cohort_id == 'COH42'
    _, kwargs = mock_api.return_value.create_cohort_from_criteria.call_args
    assert kwargs['project'] == 'ourdna'
    body = kwargs['body_create_cohort_from_criteria']
    assert body.cohort_criteria.sg_ids_internal == ['CPG1', 'CPG2']
    # The server rejects an SG list combined with any other criterion (SET-839).
    assert getattr(body.cohort_criteria, 'projects', None) in (None, [])
    assert body.cohort_spec.name == 'array-aggregate-2026-08-10'


@patch('popgen_genotyping.metamist_utils.CohortApi')
@patch('popgen_genotyping.metamist_utils.metamist_project')
def test_create_custom_cohort_no_id_raises(mock_project, mock_api):
    mock_project.return_value = 'ourdna'
    # A malformed response without a cohort ID cannot be expressed as a NewCohort
    # (cohort_id is required there), so a bare dict stands in for it.
    mock_api.return_value.create_cohort_from_criteria.return_value = {}

    with pytest.raises(ValueError, match='no cohort ID'):
        create_custom_cohort(name='x', description='y', sg_ids=['CPG1'])


@patch('popgen_genotyping.metamist_utils.CohortApi')
@patch('popgen_genotyping.metamist_utils.metamist_project')
def test_create_custom_cohort_missing_sgs_raises(mock_project, mock_api):
    mock_project.return_value = 'ourdna'
    mock_api.return_value.create_cohort_from_criteria.return_value = NewCohort(
        cohort_id='COH42',
        sequencing_group_ids=['CPG1'],
        dry_run=False,
    )

    with pytest.raises(ValueError, match=r'missing 1 requested sequencing group.*CPG9'):
        create_custom_cohort(name='x', description='y', sg_ids=['CPG1', 'CPG9'])


@patch('popgen_genotyping.metamist_utils.time')
@patch('popgen_genotyping.metamist_utils.query_cohorts_with_analyses')
def test_wait_for_cohort_analyses_returns_when_registered(mock_query, mock_time):
    mock_time.monotonic.return_value = 0
    mock_query.return_value = [
        _cohort('COH1', ['CPG1'], analyses=[_analysis('array_cohort_bed', 'x'), _analysis('array_bafregress', 'y')]),
    ]

    wait_for_cohort_analyses(cohort_ids=['COH1'], analysis_types=['array_cohort_bed', 'array_bafregress'])

    mock_time.sleep.assert_not_called()


@patch('popgen_genotyping.metamist_utils.time')
@patch('popgen_genotyping.metamist_utils.query_cohorts_with_analyses')
def test_wait_for_cohort_analyses_polls_until_registered(mock_query, mock_time):
    mock_time.monotonic.side_effect = [0, 10]
    mock_query.side_effect = [
        [_cohort('COH1', ['CPG1'], analyses=[_analysis('array_cohort_bed', 'x')])],
        [_cohort('COH1', ['CPG1'], analyses=[_analysis('array_cohort_bed', 'x'), _analysis('array_bafregress', 'y')])],
    ]

    wait_for_cohort_analyses(cohort_ids=['COH1'], analysis_types=['array_cohort_bed', 'array_bafregress'])

    mock_time.sleep.assert_called_once()


@patch('popgen_genotyping.metamist_utils.time')
@patch('popgen_genotyping.metamist_utils.query_cohorts_with_analyses')
def test_wait_for_cohort_analyses_times_out(mock_query, mock_time):
    mock_time.monotonic.side_effect = [0, 1000]
    mock_query.return_value = [_cohort('COH1', ['CPG1'], analyses=[_analysis('array_cohort_bed', 'x')])]

    with pytest.raises(TimeoutError, match=r'COH1.*array_bafregress'):
        wait_for_cohort_analyses(cohort_ids=['COH1'], analysis_types=['array_cohort_bed', 'array_bafregress'])

    mock_time.sleep.assert_not_called()
