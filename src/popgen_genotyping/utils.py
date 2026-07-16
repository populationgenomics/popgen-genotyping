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


def get_output_prefix(dataset: Dataset, stage_name: str, tmp: bool = False) -> Path:
    """
    Standardised output prefix for all stages.

    Args:
        dataset (Dataset): The CPG Dataset object.
        stage_name (str): Name of the current pipeline stage.
        tmp (bool): If True, returns a path in the 'tmp' category. Defaults to False.

    Returns:
        Path: The resolved cloud path prefix.
    """
    version: int = config_retrieve(['workflow', 'version'], 1)
    prefix: Path = dataset.prefix(category='tmp' if tmp else 'default')
    return prefix / get_workflow().name / stage_name / str(version)


def reconstruct_cohort_output_paths(
    dataset: Dataset,
    stage_name: str,
    cohort_ids: list[str],
    suffixes: dict[str, str],
    tmp: bool = False,
) -> list[dict[str, str]]:
    """
    Rebuild a prior CohortStage's per-cohort output paths deterministically.

    Used by the phase-2 (super-cohort) stages to locate phase-1 per-plate outputs
    without cross-cohort stage wiring or a Metamist lookup. The reconstruction is
    kept in lockstep with production output naming by reusing ``get_output_prefix``;
    it therefore depends on both phases sharing the same ``workflow.version`` and
    dataset.

    Args:
        dataset (Dataset): The CPG Dataset the phase-1 outputs were written under.
        stage_name (str): Name of the phase-1 stage that produced the outputs.
        cohort_ids (list[str]): Component cohort IDs to reconstruct paths for.
        suffixes (dict[str, str]): Mapping of output key to filename suffix, e.g.
            ``{'bed': '.bed', 'bim': '.bim', 'fam': '.fam'}`` or
            ``{'txt': '.BAFRegress.txt'}``.
        tmp (bool): Whether the phase-1 stage wrote to 'tmp' storage. Defaults to False.

    Returns:
        list[dict[str, str]]: One dict per cohort ID mapping each key to its cloud path.
    """
    prefix: Path = get_output_prefix(dataset=dataset, stage_name=stage_name, tmp=tmp)
    return [{key: str(prefix / f'{cohort_id}{sfx}') for key, sfx in suffixes.items()} for cohort_id in cohort_ids]


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
