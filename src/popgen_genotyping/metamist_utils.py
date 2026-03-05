"""
Metamist GraphQL query utilities for the genotyping pipeline.
"""

from metamist.graphql import gql, query
from cpg_utils.config import config_retrieve

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
