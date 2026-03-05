"""
Job logic for merging individual BCFs into a cohort PLINK2 dataset.
"""

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

if TYPE_CHECKING:
    from hailtop.batch.job import Job


def run_cohort_bcf_to_plink(
    manifest_path: str,
    output_prefix: str,
    job_name: str = 'cohort_bcf_to_plink',
) -> 'Job':
    """
    Run the PLINK2 conversion and merging orchestration script.

    Args:
        manifest_path (str): Cloud path to the JSON manifest.
        output_prefix (str): Cloud prefix for the merged PLINK2 files.
        job_name (str): Name for the Hail Batch job.

    Returns:
        Job: A Hail Batch job object.
    """
    b = get_batch()
    j = b.new_job(name=job_name)

    j.image(config_retrieve(['workflow', 'driver_image']))
    j.cpu(8)
    j.memory('highmem')
    j.storage('50G')

    j.command(
        f"""
        set -ex

        # Use the orchestration script
        python3 -m popgen_genotyping.scripts.vcf_to_plink \\
            --manifest {manifest_path} \\
            --out-prefix {output_prefix} \\
            --threads 8
        """
    )

    return j
