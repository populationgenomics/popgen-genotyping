"""
Metamist GraphQL query utilities for the genotyping pipeline.
"""

import csv
import functools
import time
from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

from cpg_utils import to_path
from loguru import logger
from cpg_utils.config import config_retrieve
from metamist.apis import CohortApi
from metamist.graphql import gql, query
from metamist.models import BodyCreateCohortFromCriteria, CohortBody, CohortCriteria

from popgen_genotyping.utils import get_sequencing_group_cohort

if TYPE_CHECKING:
    from cpg_flow.targets import Cohort, SequencingGroup

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

# GQL to list every cohort in a project with its membership and registered analyses.
# Analyses register against the cohort (with an empty sequencing_group_ids list), not
# against individual SGs, so per-SG discovery is impossible — we enumerate cohorts and
# invert cohort membership to SGs. Used to resolve per-plate `array_cohort_bed` and
# `array_bafregress` outputs, and the previous `array_aggregate_pgen` aggregate.
QUERY_COHORTS_WITH_ANALYSES = gql(
    """
    query CohortsWithAnalysesQuery($project: String!) {
      project(name: $project) {
        cohorts {
          id
          name
          sequencingGroups {
            id
          }
          analyses {
            id
            type
            timestampCompleted
            outputs
          }
        }
      }
    }
    """
)


def metamist_project(project: str | None = None) -> str:
    """
    Resolve the namespaced Metamist project name for the current run.

    Args:
        project (str, optional): Metamist project name. Defaults to the 'dataset' from config.

    Returns:
        str: The project name, with a -test suffix at test access level.
    """
    if project is None:
        project = config_retrieve(['workflow', 'dataset'])

    # At test access level the namespaced Metamist project carries a -test suffix.
    if config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in project:
        project += '-test'

    return project


def query_genotyping_manifests(project: str | None = None) -> list[dict[str, Any]]:
    """
    Query Metamist for genotyping manifest analyses.

    Args:
        project (str, optional): Metamist project name. Defaults to the 'dataset' from config.

    Returns:
        list[dict[str, Any]]: List of 'outputs' dictionaries from manifest analyses.
    """
    project = metamist_project(project)

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
                mapping[sg_id] = {'gtc': gtc_path, 'old_name': f'{barcode}_{pos}'}

    return mapping


def resolve_cohort_gtc_mapping(cohort: 'Cohort') -> dict[str, dict[str, str]]:
    """
    Resolve GTC paths and reheadering names for all sequencing groups in a cohort.

    Args:
        cohort (Cohort): The cohort to resolve.

    Returns:
        dict[str, dict[str, str]]: Mapping of SG ID to {'gtc': path, 'old_name': barcode_pos}

    Raises:
        ValueError: If no manifest is found for the cohort.
    """
    # 1. Query manifests for the project
    all_manifests: list[dict[str, Any]] = query_genotyping_manifests(cohort.dataset.name)

    # 2. Find manifest with the cohort ID in its basename
    matching_manifests: list[dict[str, Any]] = [m for m in all_manifests if cohort.id in str(m.get('basename', ''))]

    if not matching_manifests:
        raise ValueError(f'No manifest found for cohort {cohort.id}')

    matching_manifests.sort(key=lambda x: int(x.get('id', 0)), reverse=True)
    manifest_path: str = matching_manifests[0]['path']

    # 3. Parse manifest
    mapping: dict[str, dict[str, str]] = parse_genotyping_manifest_for_reheader(manifest_path)

    # 4. Filter to only include SGs active in this cohort
    cohort_sg_ids: set[str] = set(cohort.get_sequencing_group_ids())
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
    cohort: Cohort = get_sequencing_group_cohort(sequencing_group)
    cohort_id: str = cohort.id  # Expected format: COH[0-9]+

    # 3. Find manifest with the cohort ID in its basename
    matching_manifests: list[dict[str, Any]] = [m for m in all_manifests if cohort_id in str(m.get('basename', ''))]

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
        raise ValueError(f'Manifest for cohort {cohort_id} (ID: {latest_manifest.get("id")}) has no path output')

    # 4. Parse the chosen manifest and retrieve the path
    mapping: dict[str, str] = parse_genotyping_manifest(manifest_path)

    if sequencing_group.id not in mapping:
        raise ValueError(
            f'Sequencing group {sequencing_group.id} not found in genotyping manifest '
            f'{manifest_path}. Manifest contains {len(mapping)} entries.'
        )

    return mapping[sequencing_group.id]


def query_reported_sex(project: str | None = None) -> dict[str, str]:
    """
    Query Metamist for reported sex of participants, mapped to Sequencing Group IDs.

    Args:
        project (str, optional): Metamist project name. Defaults to the 'dataset' from config.

    Returns:
        dict[str, str]: Mapping of Sequencing Group ID to reported sex.
    """
    project = metamist_project(project)

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


def _extract_output_path(outputs: Any) -> str | None:
    """
    Extract the registered file path from an analysis ``outputs`` field.

    Metamist analyses store ``outputs`` either as a bare path string or as a dict
    carrying a ``path`` key; this normalises both.

    Args:
        outputs (Any): The analysis ``outputs`` value.

    Returns:
        str | None: The output path, or None if absent.
    """
    if not outputs:
        return None
    if isinstance(outputs, str):
        return outputs
    if isinstance(outputs, dict):
        return outputs.get('path')
    return None


def query_cohorts_with_analyses(project: str | None = None) -> list[dict[str, Any]]:
    """
    Query Metamist for every cohort in a project with its membership and analyses.

    Args:
        project (str, optional): Metamist project name. Defaults to the 'dataset' from config.

    Returns:
        list[dict[str, Any]]: Cohort dicts, each with 'id', 'name', 'sequencingGroups'
            and 'analyses'. Empty list if the project has no cohorts.
    """
    project = metamist_project(project)

    query_result: dict[str, Any] = query(QUERY_COHORTS_WITH_ANALYSES, {'project': project})

    project_result = query_result.get('project') or {}
    return project_result.get('cohorts') or []


def wait_for_cohort_analyses(
    cohort_ids: Iterable[str],
    analysis_types: Iterable[str],
    project: str | None = None,
    timeout_seconds: float = 900,
    poll_seconds: float = 30,
) -> None:
    """
    Block until every cohort has a registered analysis of every requested type.

    cpg-flow's Metamist registration jobs are not stage dependencies (the status
    reporter creates them but never adds them to the stage's output jobs), so a job
    that consumes registered analyses can start while registration is still in flight.
    Polling closes that race; the deadline keeps a genuinely failed registration loud.

    Args:
        cohort_ids (Iterable[str]): Cohorts whose analyses must exist.
        analysis_types (Iterable[str]): Analysis types required on every cohort.
        project (str, optional): Metamist project name. Defaults to the 'dataset' from config.
        timeout_seconds (float): Deadline before giving up. Defaults to 900.
        poll_seconds (float): Delay between Metamist queries. Defaults to 30.

    Raises:
        TimeoutError: If any (cohort, analysis type) pair is still unregistered at the deadline.
    """
    ids = list(cohort_ids)
    types = list(analysis_types)
    deadline = time.monotonic() + timeout_seconds
    while True:
        by_id = {c.get('id'): c for c in query_cohorts_with_analyses(project)}
        missing = [
            (cohort_id, analysis_type)
            for cohort_id in ids
            for analysis_type in types
            if not any(a.get('type') == analysis_type for a in ((by_id.get(cohort_id) or {}).get('analyses') or []))
        ]
        if not missing:
            return
        if time.monotonic() >= deadline:
            raise TimeoutError(
                f'Timed out after {timeout_seconds}s waiting for Metamist analysis registration: {missing}'
            )
        logger.info(f'Waiting for {len(missing)} Metamist analysis registration(s): {missing}')
        time.sleep(poll_seconds)


def _invert_cohort_analyses(cohorts: list[dict[str, Any]], analysis_type: str) -> dict[str, dict[str, str]]:
    """
    Invert cohort-first records into an SG -> analysis-output map for one analysis type.

    Cohorts carrying no analysis of ``analysis_type`` are skipped; when a cohort has
    several (e.g. a re-run), the latest by ``timestampCompleted`` wins. Each SG in a
    matching cohort maps to that cohort's output.

    Args:
        cohorts (list[dict[str, Any]]): Records from :func:`query_cohorts_with_analyses`.
        analysis_type (str): The analysis type to select (e.g. 'array_cohort_bed').

    Returns:
        dict[str, dict[str, str]]: sg_id -> {'path': output_path, 'cohort_id': cohort_id}.

    Raises:
        ValueError: If an SG resolves to more than one cohort of this type (a
            re-genotyped sample always gets a new SG ID, so this should not happen),
            or if a selected analysis has no output path.
    """
    sg_map: dict[str, dict[str, str]] = {}
    for cohort in cohorts:
        cohort_id: str = cohort.get('id', '')
        matching = [a for a in (cohort.get('analyses') or []) if a.get('type') == analysis_type]
        if not matching:
            continue

        latest = max(matching, key=lambda a: a.get('timestampCompleted') or '')
        path = _extract_output_path(latest.get('outputs'))
        if not path:
            raise ValueError(f'{analysis_type} analysis {latest.get("id")} for cohort {cohort_id} has no output path')

        for sg in cohort.get('sequencingGroups') or []:
            sg_id = sg.get('id')
            if not sg_id:
                continue
            existing = sg_map.get(sg_id)
            if existing and existing['cohort_id'] != cohort_id:
                raise ValueError(
                    f'Sequencing group {sg_id} maps to multiple {analysis_type} cohorts: '
                    f'{existing["cohort_id"]} and {cohort_id}'
                )
            sg_map[sg_id] = {'path': path, 'cohort_id': cohort_id}

    return sg_map


def _plink_fileset_from_bed(bed_path: str, cohort_id: str) -> dict[str, str]:
    """
    Build a PLINK 1.9 fileset dict from a registered ``.bed`` anchor path.

    Args:
        bed_path (str): The registered ``.bed`` path.
        cohort_id (str): The owning cohort ID (for error messages / provenance).

    Returns:
        dict[str, str]: {'bed', 'bim', 'fam', 'cohort_id'}.

    Raises:
        ValueError: If ``bed_path`` does not end in '.bed'.
    """
    if not bed_path.endswith('.bed'):
        raise ValueError(f'array_cohort_bed output for cohort {cohort_id} is not a .bed path: {bed_path}')
    stem = bed_path.removesuffix('.bed')
    return {'bed': bed_path, 'bim': f'{stem}.bim', 'fam': f'{stem}.fam', 'cohort_id': cohort_id}


def resolve_cohort_bed_map(
    project: str | None = None,
    cohorts: list[dict[str, Any]] | None = None,
) -> dict[str, dict[str, str]]:
    """
    Map each sequencing group to its registered per-plate PLINK 1.9 fileset.

    Discovery is cohort-first over ``array_cohort_bed`` analyses (registered against the
    plate cohort). The ``.bed`` is the registered anchor; ``.bim``/``.fam`` are derived
    by suffix.

    Args:
        project (str, optional): Metamist project name. Defaults to the 'dataset' from config.
        cohorts (list[dict[str, Any]], optional): Pre-fetched cohort records to avoid a
            second query; fetched via :func:`query_cohorts_with_analyses` if omitted.

    Returns:
        dict[str, dict[str, str]]: sg_id -> {'bed', 'bim', 'fam', 'cohort_id'}.
    """
    if cohorts is None:
        cohorts = query_cohorts_with_analyses(project)

    inverted = _invert_cohort_analyses(cohorts, 'array_cohort_bed')
    return {sg_id: _plink_fileset_from_bed(rec['path'], rec['cohort_id']) for sg_id, rec in inverted.items()}


def resolve_bafregress_map(
    sg_ids: Iterable[str],
    project: str | None = None,
    cohorts: list[dict[str, Any]] | None = None,
) -> dict[str, str]:
    """
    Map each requested sequencing group to its registered BafRegress output path.

    Full-membership discovery for the QC report: enumerates every ``array_bafregress``
    analysis in the project (cohort-first, since BafRegress registers against the cohort
    with an empty sequencing_group_ids list), inverts to sg_id -> path, and selects the
    requested SGs.

    Args:
        sg_ids (Iterable[str]): Sequencing groups to resolve (the aggregate's membership).
        project (str, optional): Metamist project name. Defaults to the 'dataset' from config.
        cohorts (list[dict[str, Any]], optional): Pre-fetched cohort records.

    Returns:
        dict[str, str]: sg_id -> BafRegress output path.

    Raises:
        ValueError: If any requested SG has no BafRegress analysis.
    """
    if cohorts is None:
        cohorts = query_cohorts_with_analyses(project)

    inverted = _invert_cohort_analyses(cohorts, 'array_bafregress')

    resolved: dict[str, str] = {}
    missing: list[str] = []
    for sg_id in sg_ids:
        rec = inverted.get(sg_id)
        if rec is None:
            missing.append(sg_id)
        else:
            resolved[sg_id] = rec['path']

    if missing:
        raise ValueError(f'No BafRegress analysis found for sequencing groups: {sorted(missing)}')

    return resolved


def resolve_merge_inputs(
    super_cohort_sg_ids: Iterable[str],
    previous_aggregate_cohort_id: str | None,
    project: str | None = None,
) -> dict[str, Any]:
    """
    Resolve phase-2 rolling-merge inputs from Metamist (no path reconstruction).

    The super cohort's membership is the source of truth. Given the previous aggregate
    *cohort* ID, its membership (AGG) and its ``array_aggregate_pgen`` PGEN are carried
    forward. Plate cohorts owning any delta SG (``NEW = ALL - AGG``) are pulled whole from
    their ``array_cohort_bed`` filesets; the merge job then applies a final ``--keep`` to
    super-cohort membership, trimming any withdrawn/excluded SGs. A plate is pulled only if
    it contributes a new SG, so already-aggregated plates are never re-merged
    (aggregate-priority). Withdrawn SGs (``AGG - ALL``) are flagged for removal from the
    carried aggregate.

    Args:
        super_cohort_sg_ids (Iterable[str]): Active membership of the super cohort.
        previous_aggregate_cohort_id (str, optional): Cohort ID of the previous aggregate,
            or None for a from-scratch (bootstrap) build.
        project (str, optional): Metamist project name. Defaults to the 'dataset' from config.

    Returns:
        dict[str, Any]: {
            'previous_aggregate_paths': {'pgen', 'pvar', 'psam'} | None,
            'samples_to_remove': list[str],   # withdrawn SGs to drop from the aggregate
            'plate_merge_list': list[dict],   # {'bed','bim','fam','cohort_id','new_count'}
            'super_cohort_size': int,
        }

    Raises:
        ValueError: If the previous aggregate cohort is missing / has no aggregate analysis,
            or if a new SG has no registered array_cohort_bed fileset.
    """
    all_sgs: set[str] = set(super_cohort_sg_ids)
    cohorts = query_cohorts_with_analyses(project)

    agg_sgs: set[str] = set()
    previous_aggregate_paths: dict[str, str] | None = None

    if previous_aggregate_cohort_id:
        agg_cohort = next((c for c in cohorts if c.get('id') == previous_aggregate_cohort_id), None)
        if agg_cohort is None:
            raise ValueError(f'Previous aggregate cohort {previous_aggregate_cohort_id} not found in project')

        agg_sgs = {sg['id'] for sg in (agg_cohort.get('sequencingGroups') or []) if sg.get('id')}

        agg_analyses = [a for a in (agg_cohort.get('analyses') or []) if a.get('type') == 'array_aggregate_pgen']
        if not agg_analyses:
            raise ValueError(
                f'Cohort {previous_aggregate_cohort_id} has no array_aggregate_pgen analysis to roll forward'
            )
        latest_agg = max(agg_analyses, key=lambda a: a.get('timestampCompleted') or '')
        pgen = _extract_output_path(latest_agg.get('outputs'))
        if not pgen or not pgen.endswith('.pgen'):
            raise ValueError(
                f'array_aggregate_pgen analysis {latest_agg.get("id")} for cohort '
                f'{previous_aggregate_cohort_id} has no valid .pgen output path'
            )
        stem = pgen.removesuffix('.pgen')
        previous_aggregate_paths = {'pgen': pgen, 'pvar': f'{stem}.pvar', 'psam': f'{stem}.psam'}

    new_sgs = all_sgs - agg_sgs
    samples_to_remove = sorted(agg_sgs - all_sgs)

    # Resolve the delta to the plate cohorts that own the new SGs. Whole filesets are
    # pulled; the merge job applies a final --keep to super-cohort membership, which trims
    # any withdrawn/excluded SGs a plate happens to carry. A plate is only pulled if it
    # contributes a new SG, so already-aggregated plates are never re-merged
    # (aggregate-priority). New plates never overlap AGG (a re-genotyped sample gets a new
    # SG ID), so the merge cannot hit a duplicate-sample collision with the carried aggregate.
    bed_map = resolve_cohort_bed_map(cohorts=cohorts)
    plate_groups: dict[str, dict[str, Any]] = {}
    missing: list[str] = []
    for sg_id in new_sgs:
        fileset = bed_map.get(sg_id)
        if fileset is None:
            missing.append(sg_id)
            continue
        group = plate_groups.setdefault(fileset['cohort_id'], {'fileset': fileset, 'new_count': 0})
        group['new_count'] += 1

    if missing:
        raise ValueError(f'No array_cohort_bed analysis found for new sequencing groups: {sorted(missing)}')

    plate_merge_list: list[dict[str, Any]] = [
        {
            'bed': group['fileset']['bed'],
            'bim': group['fileset']['bim'],
            'fam': group['fileset']['fam'],
            'cohort_id': group['fileset']['cohort_id'],
            'new_count': group['new_count'],
        }
        for group in plate_groups.values()
    ]

    return {
        'previous_aggregate_paths': previous_aggregate_paths,
        'samples_to_remove': samples_to_remove,
        'plate_merge_list': plate_merge_list,
        'super_cohort_size': len(all_sgs),
    }


def format_merge_plan(resolved: dict[str, Any], previous_aggregate_cohort_id: str | None = None) -> str:
    """
    Render a human-readable summary of a resolved rolling-merge plan for the driver log.

    Lets an operator confirm at a glance that phase 2 picked up the phase-1 plates: it
    lists each contributing plate cohort with its new-SG count. The printed totals are
    all derived from the super cohort's membership, so they are informational, not a
    cross-check — the expected merged total equals the super cohort size by construction.
    The check that the merged data actually matches is the in-job kept-sample-count
    assert after the final ``--keep``.

    Args:
        resolved (dict[str, Any]): The dict returned by :func:`resolve_merge_inputs`.
        previous_aggregate_cohort_id (str, optional): Cohort ID rolled forward, for the header.

    Returns:
        str: A multi-line summary suitable for ``logger.info``.
    """
    plate_merge_list: list[dict[str, Any]] = resolved['plate_merge_list']
    samples_to_remove: list[str] = resolved['samples_to_remove']
    super_size: int = resolved['super_cohort_size']

    new_count: int = sum(p['new_count'] for p in plate_merge_list)
    carried_forward: int = super_size - new_count
    expected_total: int = carried_forward + new_count

    header = (
        f'Rolling merge plan (previous aggregate cohort {previous_aggregate_cohort_id}):'
        if previous_aggregate_cohort_id
        else 'Rolling merge plan (bootstrap — no previous aggregate):'
    )
    lines: list[str] = [
        header,
        f'  super cohort:          {super_size} SGs',
        f'  carried forward:       {carried_forward} SGs (from previous aggregate)',
        f'  withdrawn:             {len(samples_to_remove)} SGs'
        + (f' {samples_to_remove}' if samples_to_remove else ''),
        f'  new:                   {new_count} SGs from {len(plate_merge_list)} plate cohort(s):',
    ]
    for plate in sorted(plate_merge_list, key=lambda p: p['cohort_id']):
        lines.append(f'       {plate["cohort_id"]}:  {plate["new_count"]} new SGs')
    lines.append(f'  expected merged total: {expected_total} SGs')
    return '\n'.join(lines)


def resolve_super_cohort_membership(
    plate_sg_ids: Iterable[str],
    previous_aggregate_cohort_id: str | None,
    project: str | None = None,
    cohorts: list[dict[str, Any]] | None = None,
) -> list[str]:
    """
    Compute the target super-cohort membership for a phase-2 run.

    The super cohort is the union of this phase-1 run's plate SGs and the previous
    aggregate cohort's current membership (so SGs withdrawn from the previous aggregate
    cohort since it was built are excluded automatically).

    Args:
        plate_sg_ids (Iterable[str]): SGs of the phase-1 plate cohorts.
        previous_aggregate_cohort_id (str, optional): Cohort ID of the previous aggregate,
            or None for a from-scratch (bootstrap) build.
        project (str, optional): Metamist project name. Defaults to the 'dataset' from config.
        cohorts (list[dict[str, Any]], optional): Pre-fetched cohort records.

    Returns:
        list[str]: Sorted super-cohort SG IDs.

    Raises:
        ValueError: If the previous aggregate cohort is not found in the project, or if
            every plate SG is already a member of it (nothing new to aggregate).
    """
    sg_ids: set[str] = set(plate_sg_ids)

    if previous_aggregate_cohort_id:
        if cohorts is None:
            cohorts = query_cohorts_with_analyses(project)
        agg_cohort = next((c for c in cohorts if c.get('id') == previous_aggregate_cohort_id), None)
        if agg_cohort is None:
            raise ValueError(f'Previous aggregate cohort {previous_aggregate_cohort_id} not found in project')
        agg_sg_ids = {sg['id'] for sg in (agg_cohort.get('sequencingGroups') or []) if sg.get('id')}
        # A super cohort identical to the previous aggregate would send phase 2 chasing
        # an empty new-SG set; fail here with the real reason instead.
        if sg_ids <= agg_sg_ids:
            raise ValueError(
                f'All {len(sg_ids)} plate SGs are already members of previous aggregate cohort '
                f'{previous_aggregate_cohort_id}; nothing new to aggregate'
            )
        sg_ids |= agg_sg_ids

    return sorted(sg_ids)


def find_cohort_by_membership(
    sg_ids: Iterable[str],
    project: str | None = None,
    cohorts: list[dict[str, Any]] | None = None,
) -> dict[str, Any] | None:
    """
    Find an existing cohort whose membership is exactly ``sg_ids``.

    Used to make super-cohort creation idempotent: a re-run of phase 1 (e.g. after a
    failed phase-2 submission) reuses the cohort it created last time instead of
    registering a duplicate.

    Args:
        sg_ids (Iterable[str]): The target membership.
        project (str, optional): Metamist project name. Defaults to the 'dataset' from config.
        cohorts (list[dict[str, Any]], optional): Pre-fetched cohort records.

    Returns:
        dict[str, Any] | None: The matching cohort record (latest by ID if several
            share the membership), or None.
    """
    if cohorts is None:
        cohorts = query_cohorts_with_analyses(project)

    target: set[str] = set(sg_ids)
    matches = [c for c in cohorts if {sg['id'] for sg in (c.get('sequencingGroups') or []) if sg.get('id')} == target]
    if not matches:
        return None
    # Cohort IDs are 'COH' + an unpadded integer (+ check digit), so 'latest' must
    # compare numerically: lexicographic max would rank COH99 above COH100.
    return max(matches, key=lambda c: int(str(c.get('id', '')).removeprefix('COH')))


def create_custom_cohort(name: str, description: str, sg_ids: Iterable[str], project: str | None = None) -> str:
    """
    Create a custom Metamist cohort from explicit sequencing group IDs.

    Args:
        name (str): Cohort name (must not collide with an existing cohort).
        description (str): Cohort description (provenance: source plates, previous aggregate).
        sg_ids (Iterable[str]): The cohort membership.
        project (str, optional): Metamist project name. Defaults to the 'dataset' from config.

    Returns:
        str: The new cohort ID.

    Raises:
        ValueError: If the created cohort is missing any requested SG (the cohort would
            be silently smaller than intended), or Metamist did not return a cohort ID.
    """
    project = metamist_project(project)

    requested = sorted(sg_ids)
    body = BodyCreateCohortFromCriteria(
        cohort_spec=CohortBody(name=name, description=description),
        # sg_ids_internal must be the sole criterion: the server rejects an explicit SG
        # list combined with any other criterion, including projects (metamist SET-839).
        cohort_criteria=CohortCriteria(sg_ids_internal=requested),
    )
    result = CohortApi().create_cohort_from_criteria(project=project, body_create_cohort_from_criteria=body)

    def _field(key: str) -> Any:
        return result.get(key) if isinstance(result, dict) else getattr(result, key, None)

    cohort_id = _field('cohort_id')
    if not cohort_id:
        raise ValueError(f'Cohort creation for {name!r} returned no cohort ID: {result}')

    # Current Metamist raises on inactive SGs, but a cohort whose membership differs
    # from the request would ship a wrong aggregate, so verify the created membership
    # regardless of server version.
    created = set(_field('sequencing_group_ids') or [])
    missing = sorted(set(requested) - created)
    if missing:
        raise ValueError(
            f'Cohort {name!r} ({cohort_id}) is missing {len(missing)} requested sequencing group(s): {missing}'
        )
    return str(cohort_id)
