"""
Create the phase-2 super cohort in Metamist and submit the second-phase workflow.

Phase 2 must run against a cohort that does not exist when phase 1 builds its DAG
(cpg-flow cannot create or validate a cohort mid-run). This job sidesteps that by doing
both steps at batch runtime, after every plate cohort has registered its outputs: it
creates the super cohort, then submits a fresh analysis-runner run whose driver sees the
cohort already existing. Submission from a batch job requires analysis-runner >= 3.3.0,
which auto-confirms the full-access prompt when stdin is not a TTY.
"""

from typing import TYPE_CHECKING

from cpg_utils import config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import PythonJob

# Run-specific workflow config keys that must not leak into the phase-2 submission:
# ar-guid is unique per run, and any phase-1 stage selection would mask phase-2 stages.
NON_PORTABLE_WORKFLOW_KEYS = ('ar-guid', 'first_stages', 'last_stages', 'only_stages', 'skip_stages')


def create_super_cohort_and_submit(
    plate_sg_ids: list[str],
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
    import logging  # noqa: PLC0415

    import toml  # noqa: PLC0415
    from analysis_runner.cli_analysisrunner import run_analysis_runner  # noqa: PLC0415
    from cpg_utils import to_path  # noqa: PLC0415

    from popgen_genotyping import metamist_utils  # noqa: PLC0415

    # 1. Resolve the target membership and check for an existing identical cohort.
    cohorts = metamist_utils.query_cohorts_with_analyses()
    membership = metamist_utils.resolve_super_cohort_membership(
        plate_sg_ids=plate_sg_ids,
        previous_aggregate_cohort_id=previous_aggregate_cohort_id,
        cohorts=cohorts,
    )

    existing = metamist_utils.find_cohort_by_membership(membership, cohorts=cohorts)
    if existing:
        cohort_id = str(existing['id'])
        logging.info(f'Reusing existing cohort {cohort_id} ({existing.get("name")}) with identical membership')
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
            f'{len(plate_sg_ids)} plate SGs (phase-1 ar-guid {config.config_retrieve(["workflow", "ar-guid"])})'
        )
        cohort_id = metamist_utils.create_custom_cohort(
            name=super_cohort_name,
            description=description,
            sg_ids=membership,
        )
        logging.info(f'Created super cohort {cohort_id} ({super_cohort_name}) with {len(membership)} SGs')

    # 2. Rewrite the run config for phase 2.
    config_dict = copy.deepcopy(dict(config._config))  # noqa: SLF001
    phase1_ar_guid = config_dict['workflow'].get('ar-guid')
    for key in NON_PORTABLE_WORKFLOW_KEYS:
        config_dict['workflow'].pop(key, None)
    config_dict['workflow']['input_cohorts'] = [cohort_id]

    with open('config.toml', 'w') as config_file:
        config_file.write(toml.dumps(config_dict))

    # 3. Record the hand-off BEFORE submitting: a crash between submission and record
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

    run_analysis_runner(
        dataset=config_dict['workflow']['dataset'],
        image=config_dict['workflow']['driver_image'],
        output_dir='',
        script=['second_workflow'],
        description=f'popgen-genotyping phase 2: aggregate cohort {cohort_id}',
        access_level=config_dict['workflow']['access_level'],
        config=['config.toml'],
        skip_repo_checkout=True,
    )

    return cohort_id


def run_submit_phase2(
    plate_sg_ids: list[str],
    previous_aggregate_cohort_id: str | None,
    super_cohort_name: str,
    output_path: str,
    job_name: str = 'SubmitPhase2',
) -> 'PythonJob':
    """
    Queue the super-cohort creation + phase-2 submission as a PythonJob in the driver image.

    Args:
        plate_sg_ids (list[str]): SGs of the phase-1 plate cohorts.
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
    j.call(
        create_super_cohort_and_submit,
        plate_sg_ids,
        previous_aggregate_cohort_id,
        super_cohort_name,
        output_path,
    )
    return j
