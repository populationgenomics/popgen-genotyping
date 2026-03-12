"""
Metamist GraphQL query utilities for the genotyping pipeline.
"""

import csv
import functools
from typing import TYPE_CHECKING, Any

from metamist.graphql import gql, query
from cpg_utils import to_path
from cpg_utils.config import config_retrieve
from popgen_genotyping.utils import get_sequencing_group_cohort

if TYPE_CHECKING:
    from cpg_flow.targets import SequencingGroup

# GQL to retrieve reported sex for active sequencing groups
QUERY_REPORTED_SEX = gql(
    """
    query ReportedSexQuery($project: String!) {
      project(name: $project) {
        sequencingGroups(activeOnly: {eq: true}) {
          id
          sample {
            participant {
              reportedSex
            }
          }
        }
      }
    }
    """
)

# GQL to retrieve manifest analysis entries
QUERY_GENOTYPING_MANIFESTS = gql(
    """
    query GenotypingManifestQuery($project: String!) {
      project(name: $project) {
        analyses(type: {eq: "manifest"}) {
          id
          type
          outputs
        }
      }
    }
    """
)

# GQL to retrieve previous aggregate metadata and active SGs
QUERY_PREVIOUS_AGGREGATE = gql(
    """
    query PreviousAggregateQuery($id: Int!) {
      analyses(id: {eq: $id}) {
        outputs
        project {
          sequencingGroups(activeOnly: {eq: true}) {
            id
          }
        }
      }
    }
    """
)


def query_genotyping_manifests(project: str | None = None) -> list[dict[str, Any]]:
    """
    Query Metamist for genotyping manifest analyses.

    Args:
        project (str, optional): Metamist project name. Defaults to the 'dataset' from config.

    Returns:
        list[dict[str, Any]]: List of 'outputs' dictionaries from manifest analyses.
    """
    if project is None:
        project = config_retrieve(['workflow', 'dataset'])

    # Execute the query
    query_result: dict[str, Any] = query(QUERY_GENOTYPING_MANIFESTS, {'project': project})

    if not query_result.get('project') or not query_result['project'].get('analyses'):
        return []

    # Filter for genotyping array manifests
    manifest_outputs: list[dict[str, Any]] = []
    for analysis in query_result['project']['analyses']:
        outputs = analysis.get('outputs')
        if not outputs or not isinstance(outputs, dict):
            continue

        basename: str = outputs.get('basename', '')
        if 'genotyping_array' in basename:
            # Include the analysis ID in the outputs for later sorting
            outputs['id'] = analysis.get('id')
            manifest_outputs.append(outputs)

    return manifest_outputs


@functools.lru_cache
def parse_genotyping_manifest(manifest_path: str) -> dict[str, str]:
    """
    Read and parse a genotyping manifest CSV.

    Args:
        manifest_path (str): Cloud path to the manifest CSV.

    Returns:
        dict[str, str]: Mapping of Sequencing Group ID to GTC filepath.
    """
    sg_to_gtc: dict[str, str] = {}
    with to_path(manifest_path).open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            sg_id: str | None = row.get('cpg_sequencing_group_id')
            gtc_path: str | None = row.get('cpg_gcp_filepath')
            if sg_id and gtc_path:
                sg_to_gtc[sg_id] = gtc_path

    return sg_to_gtc


def parse_genotyping_manifest_for_reheader(manifest_path: str) -> dict[str, dict[str, str]]:
    """
    Parse manifest to get both GTC path and reheadering mapping.

    Args:
        manifest_path (str): Cloud path to the manifest CSV.

    Returns:
        dict[str, dict[str, str]]: Mapping of SG ID to {'gtc': path, 'old_name': barcode_pos}
    """
    mapping: dict[str, dict[str, str]] = {}
    with to_path(manifest_path).open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            sg_id: str | None = row.get('cpg_sequencing_group_id')
            gtc_path: str | None = row.get('cpg_gcp_filepath')
            barcode: str | None = row.get('sentrix_barcode_a')
            pos: str | None = row.get('sentrix_position_a')
            
            if sg_id and gtc_path and barcode and pos:
                mapping[sg_id] = {
                    'gtc': gtc_path,
                    'old_name': f'{barcode}_{pos}'
                }

    return mapping


def resolve_cohort_gtc_mapping(cohort: 'Cohort') -> dict[str, dict[str, str]]:
    """
    Resolve GTC paths and reheadering names for all sequencing groups in a cohort.

    Args:
        cohort (Cohort): The cohort to resolve.

    Returns:
        dict[str, dict[str, str]]: Mapping of SG ID to {'gtc': path, 'old_name': barcode_pos}
    """
    # 1. Query manifests for the project
    all_manifests: list[dict[str, Any]] = query_genotyping_manifests(cohort.analysis_dataset.name)

    # 2. Find manifest with the cohort ID in its basename
    matching_manifests: list[dict[str, Any]] = [
        m for m in all_manifests if cohort.id in str(m.get('basename', ''))
    ]

    if not matching_manifests:
        raise ValueError(f'No manifest found for cohort {cohort.id}')

    matching_manifests.sort(key=lambda x: int(x.get('id', 0)), reverse=True)
    manifest_path: str = matching_manifests[0]['path']

    # 3. Parse manifest
    mapping: dict[str, dict[str, str]] = parse_genotyping_manifest_for_reheader(manifest_path)
    
    # 4. Filter to only include SGs active in this cohort
    cohort_sg_ids = set(cohort.get_sequencing_group_ids())
    return {sg_id: data for sg_id, data in mapping.items() if sg_id in cohort_sg_ids}


def resolve_gtc_path(sequencing_group: 'SequencingGroup') -> str:
    """
    Resolve the GTC cloud path for a sequencing group from Metamist manifests.

    Args:
        sequencing_group (SequencingGroup): The sequencing group to resolve.

    Returns:
        str: The cloud path to the GTC file.

    Raises:
        ValueError: If no manifest is found or the SG is not in the manifest.
    """
    # 1. Query all possible genotyping manifests for the project
    all_manifests: list[dict[str, Any]] = query_genotyping_manifests(sequencing_group.dataset.name)

    # 2. Resolve the cohort for this sequencing group
    cohort = get_sequencing_group_cohort(sequencing_group)
    cohort_id: str = cohort.id  # Expected format: COH[0-9]+

    # 3. Find manifest with the cohort ID in its basename
    matching_manifests: list[dict[str, Any]] = [
        m for m in all_manifests if cohort_id in str(m.get('basename', ''))
    ]

    if not matching_manifests:
        available_basenames: list[str] = [str(m.get('basename', 'unknown')) for m in all_manifests]
        raise ValueError(
            f'No manifest found for cohort {cohort_id} in project {sequencing_group.dataset.name}. '
            f'Available manifest basenames: {available_basenames}'
        )

    # Sort by ID to get the latest if multiple exist for the same cohort
    matching_manifests.sort(key=lambda x: int(x.get('id', 0)), reverse=True)
    latest_manifest: dict[str, Any] = matching_manifests[0]
    manifest_path: str | None = latest_manifest.get('path')

    if not manifest_path:
        raise ValueError(
            f'Manifest for cohort {cohort_id} (ID: {latest_manifest.get("id")}) has no path output'
        )

    # 4. Parse the chosen manifest and retrieve the path
    mapping: dict[str, str] = parse_genotyping_manifest(manifest_path)

    if sequencing_group.id not in mapping:
        raise ValueError(
            f'Sequencing group {sequencing_group.id} not found in genotyping manifest '
            f'{manifest_path}. Manifest contains {len(mapping)} entries.'
        )

    return mapping[sequencing_group.id]


def query_previous_aggregate(analysis_id: int) -> tuple[dict[str, Any], list[str]]:
    """
    Query Metamist for a previous aggregate analysis and its project's active samples.

    Args:
        analysis_id (int): The Metamist analysis ID.

    Returns:
        tuple[dict[str, Any], list[str]]: (outputs_dict, active_sg_ids)

    Raises:
        ValueError: If the analysis ID is not found.
    """
    query_result: dict[str, Any] = query(QUERY_PREVIOUS_AGGREGATE, {'id': analysis_id})

    if not query_result.get('analyses'):
        raise ValueError(f'Analysis with ID {analysis_id} not found in Metamist')

    analysis: dict[str, Any] = query_result['analyses'][0]
    outputs: dict[str, Any] = analysis.get('outputs', {})
    project: dict[str, Any] = analysis.get('project', {})
    active_sgs: list[str] = [sg['id'] for sg in project.get('sequencingGroups', [])]

    return outputs, active_sgs


def query_reported_sex(project: str | None = None) -> dict[str, str]:
    """
    Query Metamist for reported sex of participants, mapped to Sequencing Group IDs.

    Args:
        project (str, optional): Metamist project name. Defaults to the 'dataset' from config.

    Returns:
        dict[str, str]: Mapping of Sequencing Group ID to reported sex.
    """
    if project is None:
        project = config_retrieve(['workflow', 'dataset'])

    # Execute the query
    query_result: dict[str, Any] = query(QUERY_REPORTED_SEX, {'project': project})

    if not query_result.get('project') or not query_result['project'].get('sequencingGroups'):
        return {}

    # Extract mapping
    dict_samples: dict[str, str] = {}
    for sg in query_result['project']['sequencingGroups']:
        sg_id: str | None = sg.get('id')
        sample: dict[str, Any] | None = sg.get('sample')
        if not sg_id or not sample:
            continue

        participant: dict[str, Any] | None = sample.get('participant')
        if not participant:
            continue

        reported_sex: str | None = participant.get('reportedSex')
        if reported_sex:
            dict_samples[sg_id] = reported_sex

    return dict_samples


def resolve_rolling_aggregate(prev_analysis_id: int | str) -> tuple[dict[str, str], list[str]]:
    """
    Resolve paths and sample delta for a rolling multi-cohort aggregate.

    Args:
        prev_analysis_id (int | str): The Metamist analysis ID of the previous aggregate.

    Returns:
        tuple[dict[str, str], list[str]]: (previous_aggregate_paths, samples_to_remove)
    """
    # Import here to avoid circular dependency
    from popgen_genotyping.utils import parse_psam  # noqa: PLC0415

    prev_outputs, active_sg_ids = query_previous_aggregate(int(prev_analysis_id))

    # Expecting PLINK 1.9 BED/BIM/FAM in outputs
    previous_aggregate_paths: dict[str, str] = {
        'bed': prev_outputs['bed'],
        'bim': prev_outputs['bim'],
        'fam': prev_outputs['fam'],
    }

    # Parse the previous .fam to find all samples that were in the aggregate
    prev_samples: list[str] = parse_psam(previous_aggregate_paths['fam'])

    # Find samples that are in the previous aggregate but no longer active
    samples_to_remove: list[str] = list(set(prev_samples) - set(active_sg_ids))

    return previous_aggregate_paths, samples_to_remove
