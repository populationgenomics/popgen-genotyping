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
          type
          outputs
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
    cohort_id = cohort.id

    # 3. Find manifest with the cohort ID in its basename
    matching_manifests = [
        m for m in all_manifests if cohort_id in m.get('basename', '')
    ]

    if not matching_manifests:
        raise ValueError(
            f'No manifest found for cohort {cohort_id} in project {sequencing_group.dataset.name}'
        )

    # Sort by ID to get the latest if multiple exist for the same cohort
    matching_manifests.sort(key=lambda x: x.get('id', 0), reverse=True)
    latest_manifest = matching_manifests[0]
    manifest_path = latest_manifest.get('path')

    if not manifest_path:
        raise ValueError(f'Manifest for cohort {cohort_id} has no path output')

    # 4. Parse the chosen manifest and retrieve the path
    mapping = parse_genotyping_manifest(manifest_path)

    if sequencing_group.id not in mapping:
        raise ValueError(
            f'Sequencing group {sequencing_group.id} not found in genotyping manifest '
            f'{manifest_path}'
        )
        
    return mapping[sequencing_group.id]


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
