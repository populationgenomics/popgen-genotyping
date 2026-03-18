"""
Job logic for converting a PLINK2 fileset to PLINK 1.9 format.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob
    from hailtop.batch.resource import ResourceGroup


def run_plink2_to_plink1(
    pfile_prefix: dict[str, str],
    output_prefix: str,
    job_name: str = 'plink2_to_plink1',
) -> tuple[BashJob, ResourceGroup]:
    """
    Convert a PLINK2 fileset (.pgen, .pvar, .psam) to PLINK 1.9 format (.bed, .bim, .fam).

    Args:
        pfile_prefix (dict[str, str]): Cloud paths for the input PLINK2 fileset.
        output_prefix (str): Cloud prefix for the output PLINK 1.9 files.
        job_name (str): Name for the Hail Batch job.

    Returns:
        tuple[BashJob, ResourceGroup]: The job object and the output resource group.
    """
    b = get_batch()
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'merge_cohort_plink'],
        image=config_retrieve(['workflow', 'plink_image']),
        default_cpu=2,
        default_storage='50G',
    )

    # Stage the input PLINK2 files
    pfile = b.read_input_group(
        pgen=pfile_prefix['pgen'],
        pvar=pfile_prefix['pvar'],
        psam=pfile_prefix['psam'],
    )

    # Define the output resource group
    j.declare_resource_group(
        output_plink1={
            'bed': '{root}.bed',
            'bim': '{root}.bim',
            'fam': '{root}.fam',
        }
    )

    # Construct the command
    j.command(
        f"""
        set -ex
        plink2 --pfile {pfile} --make-bed --out {j.output_plink1}
        """
    )

    # Write the output to the cloud
    b.write_output(j.output_plink1, output_prefix)

    return j, j.output_plink1
