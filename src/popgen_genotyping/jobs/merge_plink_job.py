"""
Job logic for merging multiple cohort PLINK2 datasets into a single unified dataset.
"""

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

if TYPE_CHECKING:
    from hailtop.batch.job import Job


def run_merge_plink(
    cohort_plink_paths: list[dict[str, str]],
    output_prefix: str,
    job_name: str = 'merge_cohort_plink',
) -> 'Job':
    """
    Merge multiple PLINK2 datasets into a single unified dataset.

    Args:
        cohort_plink_paths (list[dict[str, str]]): List of dicts, each with 'pgen', 'pvar', 'psam' cloud paths.
        output_prefix (str): Cloud prefix for the final merged PLINK2 files.
        job_name (str): Name for the Hail Batch job.

    Returns:
        Job: A Hail Batch job object.
    """
    b = get_batch()
    j = b.new_job(name=job_name)

    j.image(config_retrieve(['workflow', 'BcfToPlink_image']))
    j.cpu(4)
    j.memory('highmem')
    j.storage('100G')

    # 1. Stage all input datasets
    staged_prefixes = []
    for _i, paths in enumerate(cohort_plink_paths):
        # We read as an input group to keep them together
        # Note: PLINK2 --pmerge-list expects prefixes
        resource = b.read_input_group(
            pgen=paths['pgen'],
            pvar=paths['pvar'],
            psam=paths['psam']
        )
        # The prefix for plink2 is just the base path without extension
        # We can use the resource root
        staged_prefixes.append(str(resource))

    # 2. Define output resource group
    j.declare_resource_group(
        output_plink={
            'pgen': '{root}.pgen',
            'pvar': '{root}.pvar',
            'psam': '{root}.psam',
        }
    )

    # 3. Construct merge list and execute
    # We use a python loop to build the shell command
    merge_list_content = '\n'.join(staged_prefixes)

    j.command(
        f"""
        set -ex

        echo "{merge_list_content}" > mergelist.txt

        plink2 --pmerge-list mergelist.txt --make-pgen --out {j.output_plink}
        """
    )

    # 4. Write outputs back to cloud
    b.write_output(j.output_plink, output_prefix)

    return j
