"""
Metamist GraphQL query utilities for the genotyping pipeline.
"""

import csv
import functools
from typing import TYPE_CHECKING

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


def query_genotyping_manifests(project: str | None = None) -> list[dict]:
    """
    Query Metamist for genotyping manifest analyses.

    Args:
        project (str, optional): Metamist project name. Defaults to the 'dataset' from config.

    Returns:
        list[dict]: List of 'outputs' dictionaries from manifest analyses containing 'genotyping_array'.
    """
    if project is None:
        project = config_retrieve(['workflow', 'dataset'])

    # Execute the query
    query_result = query(QUERY_GENOTYPING_MANIFESTS, {'project': project})

    if not query_result.get('project') or not query_result['project'].get('analyses'):
        return []

    # Filter for genotyping array manifests
    manifest_outputs = []
    for analysis in query_result['project']['analyses']:
        outputs = analysis.get('outputs')
        if not outputs or not isinstance(outputs, dict):
            continue

        basename = outputs.get('basename', '')
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
    sg_to_gtc = {}
    with to_path(manifest_path).open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            sg_id = row.get('cpg_sequencing_group_id')
            gtc_path = row.get('cpg_gcp_filepath')
            if sg_id and gtc_path:
                sg_to_gtc[sg_id] = gtc_path

    return sg_to_gtc


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
    all_manifests = query_genotyping_manifests(sequencing_group.dataset.name)

    # 2. Resolve the cohort for this sequencing group
    cohort = get_sequencing_group_cohort(sequencing_group)
    cohort_id = cohort.id  # Expected format: COH[0-9]+

    # 3. Find manifest with the cohort ID in its basename
    matching_manifests = [
        m for m in all_manifests if cohort_id in m.get('basename', '')
    ]

    if not matching_manifests:
        available_basenames = [m.get('basename', 'unknown') for m in all_manifests]
        raise ValueError(
            f'No manifest found for cohort {cohort_id} in project {sequencing_group.dataset.name}. '
            f'Available manifest basenames: {available_basenames}'
        )

    # Sort by ID to get the latest if multiple exist for the same cohort
    # Ensure ID is treated as integer for sorting
    matching_manifests.sort(key=lambda x: int(x.get('id', 0)), reverse=True)
    latest_manifest = matching_manifests[0]
    manifest_path = latest_manifest.get('path')

    if not manifest_path:
        raise ValueError(
            f'Manifest for cohort {cohort_id} (ID: {latest_manifest.get("id")}) has no path output'
        )

    # 4. Parse the chosen manifest and retrieve the path
    mapping = parse_genotyping_manifest(manifest_path)

    if sequencing_group.id not in mapping:
        raise ValueError(
            f'Sequencing group {sequencing_group.id} not found in genotyping manifest '
            f'{manifest_path}. Manifest contains {len(mapping)} entries.'
        )

    return mapping[sequencing_group.id]


def query_previous_aggregate(analysis_id: int) -> tuple[dict, list[str]]:
    """
    Query Metamist for a previous aggregate analysis and its project's active samples.

    Args:
        analysis_id (int): The Metamist analysis ID.

    Returns:
        tuple[dict, list[str]]: (outputs_dict, active_sg_ids)
    """
    query_result = query(QUERY_PREVIOUS_AGGREGATE, {'id': analysis_id})

    if not query_result.get('analyses'):
        raise ValueError(f'Analysis with ID {analysis_id} not found in Metamist')

    analysis = query_result['analyses'][0]
    outputs = analysis.get('outputs', {})
    project = analysis.get('project', {})
    active_sgs = [sg['id'] for sg in project.get('sequencingGroups', [])]

    return outputs, active_sgs


def query_reported_sex(project: str | None = None) -> dict[str, str]:
    """
    Query Metamist for reported sex of participants, mapped to Sequencing Group IDs.

    Args:
        project (str, optional): Metamist project name. Defaults to the 'dataset' from config.

    Returns:
        dict[str, str]: Mapping of Sequencing Group ID to reported sex (e.g. {'CPG123': 'Female'}).
    """
    if project is None:
        project = config_retrieve(['workflow', 'dataset'])

    # Execute the query
    query_result = query(QUERY_REPORTED_SEX, {'project': project})

    if not query_result.get('project') or not query_result['project'].get('sequencingGroups'):
        return {}

    # Extract mapping
    dict_samples = {}
    for sg in query_result['project']['sequencingGroups']:
        sg_id = sg.get('id')
        sample = sg.get('sample')
        if not sg_id or not sample:
            continue

        participant = sample.get('participant')
        if not participant:
            continue

        reported_sex = participant.get('reportedSex')
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
    from popgen_genotyping.utils import parse_psam

    prev_outputs, active_sg_ids = query_previous_aggregate(int(prev_analysis_id))

    # Expecting PLINK 1.9 BED/BIM/FAM in outputs
    previous_aggregate_paths = {
        'bed': prev_outputs['bed'],
        'bim': prev_outputs['bim'],
        'fam': prev_outputs['fam'],
    }

    # Parse the previous .fam to find all samples that were in the aggregate
    prev_samples = parse_psam(previous_aggregate_paths['fam'])

    # Find samples that are in the previous aggregate but no longer active
    samples_to_remove = list(set(prev_samples) - set(active_sg_ids))

    return previous_aggregate_paths, samples_to_remove
