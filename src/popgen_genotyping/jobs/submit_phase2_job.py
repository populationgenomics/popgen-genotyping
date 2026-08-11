"""
Create the phase-2 super cohort in Metamist and submit the second-phase workflow.

Phase 2 must run against a cohort that does not exist when phase 1 builds its DAG
(cpg-flow cannot create or validate a cohort mid-run). This job sidesteps that by doing
both steps at batch runtime, after every plate cohort has registered its outputs: it
creates the super cohort, then submits a fresh analysis-runner run whose driver sees the
cohort already existing. The submission POSTs to the analysis-runner server directly
rather than going through ``run_analysis_runner``: that helper prompts interactively for
full access and, worse, catches HTTP errors and returns normally, which would leave this
job green with phase 2 never submitted.
"""

from typing import TYPE_CHECKING

from cpg_utils import config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import PythonJob

# Run-specific workflow config keys that must not leak into the phase-2 submission:
# the ar-guid is unique per run, and any phase-1 stage selection would mask phase-2 stages.
NON_PORTABLE_WORKFLOW_KEYS = (config.AR_GUID_NAME, 'first_stages', 'last_stages', 'only_stages', 'skip_stages')


def create_super_cohort_and_submit(
    plate_sg_ids: list[str],
    plate_cohort_ids: list[str],
    previous_aggregate_cohort_id: str | None,
    super_cohort_name: str,
    output_path: str,
) -> str:
    """
    Create (or reuse) the super cohort, then submit the phase-2 workflow via analysis-runner.

    Runs inside the driver image as a PythonJob. The super cohort membership is the union
    of the phase-1 plate SGs and the previous aggregate cohort's current membership. An
    existing cohort with exactly this membership is reused, so a re-run cannot register a
    duplicate; a different cohort already using ``super_cohort_name`` is an error.

    Args:
        plate_sg_ids (list[str]): SGs of the phase-1 plate cohorts.
        plate_cohort_ids (list[str]): The phase-1 plate cohort IDs, whose durable outputs
            must be registered in Metamist before phase 2 is submitted.
        previous_aggregate_cohort_id (str, optional): Cohort ID of the previous aggregate,
            or None for a from-scratch (bootstrap) build.
        super_cohort_name (str): Name for the new super cohort.
        output_path (str): Cloud path for the submission-record sentinel TOML.

    Returns:
        str: The super cohort ID.

    Raises:
        ValueError: If ``super_cohort_name`` is taken by a cohort with different membership.
    """
    import copy  # noqa: PLC0415

    import requests  # noqa: PLC0415
    import toml  # noqa: PLC0415
    from analysis_runner.util import get_server_endpoint  # noqa: PLC0415
    from cpg_utils import to_path  # noqa: PLC0415
    from cpg_utils.cloud import get_google_identity_token  # noqa: PLC0415
    from loguru import logger  # noqa: PLC0415

    from popgen_genotyping import metamist_utils  # noqa: PLC0415

    config_dict = copy.deepcopy(dict(config.get_config()))
    phase1_ar_guid = config_dict['workflow'].get(config.AR_GUID_NAME)

    # 0. Refuse to double-submit. The stage-level skip only gates a later run's DAG
    # construction, so a forced re-run (e.g. check_expected_outputs = false) or a
    # rescheduled attempt of this job would otherwise launch a second phase-2 driver
    # racing the first on the same outputs.
    if to_path(output_path).exists():
        raise ValueError(
            f'Phase-2 submission record already exists at {output_path}; refusing to submit again. '
            f'If the recorded submission actually failed, delete the sentinel and re-run phase 1.'
        )

    # 1. Wait until every plate cohort has registered its durable outputs: cpg-flow's
    # registration jobs are not dependencies of this job (see wait_for_cohort_analyses),
    # and phase 2 resolves those analyses at driver startup.
    metamist_utils.wait_for_cohort_analyses(
        cohort_ids=plate_cohort_ids,
        analysis_types=('array_cohort_bed', 'array_bafregress'),
    )

    # 2. Resolve the target membership and check for an existing identical cohort.
    cohorts = metamist_utils.query_cohorts_with_analyses()
    membership = metamist_utils.resolve_super_cohort_membership(
        plate_sg_ids=plate_sg_ids,
        previous_aggregate_cohort_id=previous_aggregate_cohort_id,
        cohorts=cohorts,
    )

    existing = metamist_utils.find_cohort_by_membership(membership, cohorts=cohorts)
    if existing:
        cohort_id = str(existing['id'])
        logger.info(f'Reusing existing cohort {cohort_id} ({existing.get("name")}) with identical membership')
    else:
        name_clash = next((c for c in cohorts if c.get('name') == super_cohort_name), None)
        if name_clash:
            raise ValueError(
                f'Cohort name {super_cohort_name!r} is already used by {name_clash.get("id")} with different '
                f'membership; choose a new super_cohort_name'
            )
        description = (
            f'popgen-genotyping aggregate: {len(membership)} SGs = '
            f'previous aggregate {previous_aggregate_cohort_id or "none (bootstrap)"} + '
            f'{len(plate_sg_ids)} plate SGs (phase-1 ar-guid {phase1_ar_guid or "unknown"})'
        )
        cohort_id = metamist_utils.create_custom_cohort(
            name=super_cohort_name,
            description=description,
            sg_ids=membership,
        )
        logger.info(f'Created super cohort {cohort_id} ({super_cohort_name}) with {len(membership)} SGs')

    # 3. Rewrite the run config for phase 2.
    for key in NON_PORTABLE_WORKFLOW_KEYS:
        config_dict['workflow'].pop(key, None)
    config_dict['workflow']['input_cohorts'] = [cohort_id]

    # 4. Record the hand-off BEFORE submitting: a crash between submission and record
    # would otherwise let a phase-1 re-run submit phase 2 twice, and two concurrent
    # phase-2 runs race on the same outputs. The inverse failure (record written,
    # submission failed) is recoverable: delete this sentinel and re-run phase 1.
    record = {
        'super_cohort_id': cohort_id,
        'super_cohort_size': len(membership),
        'previous_aggregate_cohort_id': previous_aggregate_cohort_id or '',
        'phase1_ar_guid': phase1_ar_guid or '',
    }
    with to_path(output_path).open('w') as output_file:
        output_file.write(toml.dumps(record))

    # POST to the analysis-runner server directly and raise on any HTTP error; see the
    # module docstring for why run_analysis_runner is not used. No repo/commit fields
    # means no repo checkout: the driver image alone carries the code.
    server_endpoint = get_server_endpoint()
    response = requests.post(
        server_endpoint,
        json={
            'dataset': config_dict['workflow']['dataset'],
            # Phase 1's output prefix, so the analysis-runner audit record carries a
            # meaningful output field (the forwarded config sets the same value anyway).
            'output': config_dict['workflow']['output_prefix'],
            'accessLevel': config_dict['workflow']['access_level'],
            'script': ['second_workflow'],
            'description': f'popgen-genotyping phase 2: aggregate cohort {cohort_id}',
            'image': config_dict['workflow']['driver_image'],
            'config': config_dict,
        },
        headers={'Authorization': f'Bearer {get_google_identity_token(server_endpoint)}'},
        timeout=60,
    )
    response.raise_for_status()
    logger.info(f'Phase-2 submission accepted: {response.text}')

    return cohort_id


def run_submit_phase2(
    plate_sg_ids: list[str],
    plate_cohort_ids: list[str],
    previous_aggregate_cohort_id: str | None,
    super_cohort_name: str,
    output_path: str,
    job_name: str = 'SubmitPhase2',
) -> 'PythonJob':
    """
    Queue the super-cohort creation + phase-2 submission as a PythonJob in the driver image.

    Args:
        plate_sg_ids (list[str]): SGs of the phase-1 plate cohorts.
        plate_cohort_ids (list[str]): The phase-1 plate cohort IDs.
        previous_aggregate_cohort_id (str, optional): Cohort ID of the previous aggregate.
        super_cohort_name (str): Name for the new super cohort.
        output_path (str): Cloud path for the submission-record sentinel TOML.
        job_name (str): Name for the job.

    Returns:
        PythonJob: The queued job.
    """
    batch = hail_batch.get_batch()
    j: PythonJob = batch.new_python_job(job_name)
    j.image(config.config_retrieve(['workflow', 'driver_image']))
    # Not exactly-once on spot: a preemption between the submission POST and job
    # completion would rerun the function (and trip the sentinel check as a red job
    # even though phase 2 was submitted), so keep this job off spot instances.
    j.spot(is_spot=False)
    j.call(
        create_super_cohort_and_submit,
        plate_sg_ids,
        plate_cohort_ids,
        previous_aggregate_cohort_id,
        super_cohort_name,
        output_path,
    )
    return j
