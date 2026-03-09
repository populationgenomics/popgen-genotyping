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
    previous_aggregate_paths: dict[str, str] | None = None,
    samples_to_remove: list[str] | None = None,
    job_name: str = 'merge_cohort_plink',
) -> 'Job':
    """
    Merge multiple PLINK2 datasets into a single unified dataset, with rolling aggregate support.

    Args:
        cohort_plink_paths (list[dict[str, str]]): List of dicts, each with 'pgen', 'pvar', 'psam' cloud paths.
        output_prefix (str): Cloud prefix for the final merged PLINK2 files.
        previous_aggregate_paths (dict[str, str], optional): Paths for the previous rolling aggregate.
        samples_to_remove (list[str], optional): List of SG IDs to remove from the previous aggregate.
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

    staged_prefixes = []

    # 1. Stage the previous aggregate and filter if necessary
    if previous_aggregate_paths:
        prev_resource = b.read_input_group(
            pgen=previous_aggregate_paths['pgen'],
            pvar=previous_aggregate_paths['pvar'],
            psam=previous_aggregate_paths['psam']
        )

        if samples_to_remove:
            # We must filter the previous aggregate before merging
            # We declare a resource group for the filtered version
            j.declare_resource_group(
                filtered_prev={
                    'pgen': '{root}.pgen',
                    'pvar': '{root}.pvar',
                    'psam': '{root}.psam',
                }
            )

            remove_list_content = '\n'.join(samples_to_remove)
            j.command(
                f"""
                echo "{remove_list_content}" > samples_to_remove.txt
                plink2 --pfile {prev_resource} --remove samples_to_remove.txt --make-pgen --out {j.filtered_prev}
                """
            )
            staged_prefixes.append(str(j.filtered_prev))
        else:
            staged_prefixes.append(str(prev_resource))

    # 2. Stage all new input datasets
    for _i, paths in enumerate(cohort_plink_paths):
        resource = b.read_input_group(
            pgen=paths['pgen'],
            pvar=paths['pvar'],
            psam=paths['psam']
        )
        staged_prefixes.append(str(resource))

    # 3. Define output resource group
    j.declare_resource_group(
        output_plink={
            'pgen': '{root}.pgen',
            'pvar': '{root}.pvar',
            'psam': '{root}.psam',
        }
    )

    # 4. Construct merge list and execute
    merge_list_content = '\n'.join(staged_prefixes)

    j.command(
        f"""
        set -ex

        echo "{merge_list_content}" > mergelist.txt

        plink2 --pmerge-list mergelist.txt --make-pgen --out {j.output_plink}
        """
    )

    # 5. Write outputs back to cloud
    b.write_output(j.output_plink, output_prefix)

    return j
