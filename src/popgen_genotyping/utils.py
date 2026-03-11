"""
Standard utilities and constants for the genotyping pipeline.
"""

import csv
from pathlib import Path
from typing import TYPE_CHECKING, Any

from cpg_utils import to_path, Path as CPGPath
from cpg_utils.config import config_retrieve
from cpg_flow.workflow import get_workflow
from cpg_flow.inputs import get_multicohort

if TYPE_CHECKING:
    from hailtop.batch.job import Job
    from hailtop.batch import Batch
    from cpg_flow.targets import Dataset, SequencingGroup, Cohort


def get_output_prefix(dataset: 'Dataset', stage_name: str, tmp: bool = False) -> 'CPGPath':
    """
    Standardised output prefix for all stages.
    Format: dataset.[tmp_]prefix() / workflow.name / stage_name / version
    """
    version = config_retrieve(['workflow', 'version'], 'v1')
    prefix = dataset.prefix(category='tmp' if tmp else 'main')
    return prefix / get_workflow().name / stage_name / version


def get_sequencing_group_cohort(sequencing_group: 'SequencingGroup') -> 'Cohort':
    """
    Resolve the cohort a sequencing group belongs to by searching the multi-cohort.
    """
    multicohort = get_multicohort()
    for cohort in multicohort.get_cohorts():
        if sequencing_group.id in cohort.get_sequencing_group_ids():
            return cohort

    raise ValueError(f'Sequencing group {sequencing_group.id} not found in any cohort')


def parse_psam(psam_path: str | Path | CPGPath) -> list[str]:
    """
    Extract sequencing group IDs from a PLINK2 .psam or PLINK 1.9 .fam file.
    Sample IDs are in the first column or identified by the #IID/IID header.
    """
    ids = []
    with to_path(psam_path).open() as f:
        reader = csv.reader(f, delimiter='\t')
        header = None
        iid_idx = 0

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
            if not header and len(row) == 6:
                ids.append(row[1])
            else:
                ids.append(row[iid_idx])

    return ids


def register_job(
    batch: 'Batch',
    job_name: str,
    config_path: list[str],
    image: str | None = None,
    default_cpu: int = 1,
    default_memory: str = 'standard',
    default_storage: str = '10G',
) -> 'Job':
    """
    Initialize a Hail Batch job with standard configuration from the project config.

    Args:
        batch (Batch): The Hail Batch instance.
        job_name (str): Name for the job.
        config_path (list[str]): Path in the config TOML to retrieve resources from.
        image (str, optional): Docker image. Defaults to driver_image.
        default_cpu (int): Default CPU count if not in config.
        default_memory (str): Default memory if not in config.
        default_storage (str): Default storage if not in config.

    Returns:
        Job: The initialized job.
    """
    j = batch.new_job(name=job_name)

    if image:
        j.image(image)
    else:
        j.image(config_retrieve(['workflow', 'driver_image']))

    j.cpu(config_retrieve(config_path + ['cpu'], default_cpu))
    j.memory(config_retrieve(config_path + ['memory'], default_memory))
    j.storage(config_retrieve(config_path + ['storage'], default_storage))

    return j
