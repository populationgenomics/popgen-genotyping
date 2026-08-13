"""
Standard utilities and constants for the genotyping pipeline.
"""

from __future__ import annotations

import csv
from typing import TYPE_CHECKING

from cpg_flow.inputs import get_multicohort
from cpg_flow.workflow import get_workflow
from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve
from hailtop.batch.job import BashJob

if TYPE_CHECKING:
    from cpg_flow.targets import Cohort, Dataset, SequencingGroup
    from hailtop.batch import Batch

FAM_COLUMN_COUNT: int = 6


def get_output_prefix(dataset: Dataset, stage_name: str, tmp: bool = False, versioned: bool = True) -> Path:
    """
    Standardised output prefix for all stages.

    Args:
        dataset (Dataset): The CPG Dataset object.
        stage_name (str): Name of the current pipeline stage.
        tmp (bool): If True, returns a path in the 'tmp' category. Defaults to False.
        versioned (bool): If True, appends the ``workflow.version`` segment. Set
            False for durable per-cohort artifacts that should live at a single
            stable location, be processed once, and be reused across pipeline
            versions. With no version segment the path is overwritten in place if a
            cohort is reprocessed, so these outputs must be treated as immutable once
            an aggregate references them (cpg-flow skips regeneration while they exist).
            Defaults to True.

    Returns:
        Path: The resolved cloud path prefix.
    """
    prefix: Path = dataset.prefix(category='tmp' if tmp else 'default')
    stage_prefix: Path = prefix / get_workflow().name / stage_name
    if not versioned:
        return stage_prefix
    version: int = config_retrieve(['workflow', 'version'], 1)
    return stage_prefix / str(version)


def get_previous_aggregate_cohort_id() -> str | None:
    """
    Read the previous aggregate cohort ID from config, with bootstrap made explicit.

    The key is required so a forgotten config entry cannot silently build a
    new-plates-only aggregate: a from-scratch build must be declared with the literal
    'bootstrap' (returns None); any other value must be a Metamist cohort ID.

    Returns:
        str | None: The cohort ID to roll forward, or None for a declared bootstrap.

    Raises:
        ConfigError: If the key is missing from config.
        ValueError: If the value is neither a cohort ID (COH...) nor 'bootstrap'.
    """
    value: str = config_retrieve(['popgen_genotyping', 'merge_cohort_plink', 'previous_aggregate_cohort_id'])
    if value == 'bootstrap':
        return None
    if not isinstance(value, str) or not value.startswith('COH'):
        raise ValueError(
            f"previous_aggregate_cohort_id must be a cohort ID (COH...) or the literal 'bootstrap', got {value!r}"
        )
    return value


def validate_only_stages(only_stages: list[str], phase_stages: list, entry_point: str) -> None:
    """
    Reject a ``workflow.only_stages`` selection naming stages this entry point does not run.

    Each phase has its own entry point with a fixed stage list, so ``only_stages`` is only
    ever a within-phase subset (e.g. re-running just the QC report). cpg-flow skips stages
    by exact name match, so a typo, a wrong-cased name, or a stage from the other phase
    would otherwise be silently skipped while the rest of the selection runs.

    Args:
        only_stages (list[str]): The ``workflow.only_stages`` config value; empty runs
            every stage of the phase.
        phase_stages (list): The stage classes this entry point submits.
        entry_point (str): The entry-point name, for the error message.

    Raises:
        ValueError: If ``only_stages`` names a stage outside this entry point's phase.
    """
    known = {cls.__name__ for cls in phase_stages}
    unknown = set(only_stages) - known
    if unknown:
        raise ValueError(
            f'workflow.only_stages names stages {sorted(unknown)} that {entry_point} does not run; '
            f'its stages are {sorted(known)} (exact case). The phases run against different cohorts '
            '(new plates vs the super cohort) and each has its own entry point — see the README.'
        )


def get_sequencing_group_cohort(sequencing_group: SequencingGroup) -> Cohort:
    """
    Resolve the cohort a sequencing group belongs to by searching the multi-cohort.

    Args:
        sequencing_group (SequencingGroup): The sequencing group to search for.

    Returns:
        Cohort: The resolved Cohort object.

    Raises:
        ValueError: If the sequencing group ID is not found in any cohort.
    """
    multicohort = get_multicohort()
    for cohort in multicohort.get_cohorts():
        if sequencing_group.id in cohort.get_sequencing_group_ids():
            return cohort

    raise ValueError(f'Sequencing group {sequencing_group.id} not found in any cohort')


def parse_psam(psam_path: str | Path) -> list[str]:
    """
    Extract sequencing group IDs from a PLINK2 .psam or PLINK 1.9 .fam file.

    Args:
        psam_path (str | Path): Path to the sample metadata file.

    Returns:
        list[str]: List of sample (IID) strings found in the file.
    """
    ids: list[str] = []
    with to_path(psam_path).open() as f:
        reader = csv.reader(f, delimiter='\t')
        header: list[str] | None = None
        iid_idx: int = 0

        for row in reader:
            if not row or not row[0].strip():
                continue

            # Handle PLINK2 .psam header
            if row[0].startswith('#'):
                header = [c.lstrip('#') for c in row]
                if 'IID' in header:
                    iid_idx = header.index('IID')
                continue

            # PLINK 1.9 .fam files have exactly 6 columns, Sample ID is the 2nd (index 1)
            # FamilyID SampleID PatID MatID Sex Pheno
            if not header and len(row) == FAM_COLUMN_COUNT:
                ids.append(row[1])
            else:
                ids.append(row[iid_idx])

    return ids


def register_job(
    batch: Batch,
    job_name: str,
    config_path: list[str],
    image: str | None = None,
    default_cpu: int = 1,
    default_memory: str = 'standard',
    default_storage: str = '10G',
) -> BashJob:
    """
    Initialize a Hail Batch job with standard configuration from the project config.

    Args:
        batch (Batch): The Hail Batch instance.
        job_name (str): Name for the job.
        config_path (list[str]): Path in the config TOML to retrieve resources from.
        image (str, optional): Docker image path. Defaults to driver_image from config.
        default_cpu (int): Default CPU count if not found in config. Defaults to 1.
        default_memory (str): Default memory if not found in config. Defaults to 'standard'.
        default_storage (str): Default storage if not found in config. Defaults to '10G'.

    Returns:
        BashJob: The initialized Hail Batch BashJob object.
    """
    j: BashJob = batch.new_job(name=job_name)
    assert isinstance(j, BashJob)

    if image:
        j.image(image)
    else:
        j.image(config_retrieve(['workflow', 'driver_image']))

    j.cpu(config_retrieve([*config_path, 'cpu'], default_cpu))
    j.memory(config_retrieve([*config_path, 'memory'], default_memory))
    j.storage(config_retrieve([*config_path, 'storage'], default_storage))

    return j
