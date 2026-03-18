"""
Job logic for merging multiple cohort PLINK 1.9 datasets into a single unified dataset.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from cpg_utils import to_path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def run_merge_plink(
    cohort_plink_paths: list[dict[str, str]],
    output_prefix: str,
    previous_aggregate_paths: dict[str, str] | None = None,
    samples_to_remove: list[str] | None = None,
    job_name: str = 'merge_cohort_plink',
) -> BashJob:
    """
    Merge multiple PLINK 1.9 datasets into a single unified dataset, with rolling aggregate support.

    Args:
        cohort_plink_paths (list[dict[str, str]]): List of dicts, each with 'bed', 'bim', 'fam' cloud paths.
        output_prefix (str): Cloud prefix for the final merged PLINK 1.9 files.
        previous_aggregate_paths (dict[str, str], optional): Paths for the previous rolling aggregate.
        samples_to_remove (list[str], optional): List of SG IDs to remove from the previous aggregate.
        job_name (str): Name for the Hail Batch job.

    Returns:
        Job: A Hail Batch job object.
    """
    b = get_batch()
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'merge_cohort_plink'],
        image=config_retrieve(['workflow', 'plink_image']),
        default_cpu=4,
        default_storage='100G',
    )

    staged_prefixes = []

    # 1. Stage the previous aggregate and filter if necessary
    if previous_aggregate_paths:
        prev_resource = b.read_input_group(
            bed=previous_aggregate_paths['bed'],
            bim=previous_aggregate_paths['bim'],
            fam=previous_aggregate_paths['fam'],
        )

        if samples_to_remove:
            # We must filter the previous aggregate before merging
            j.declare_resource_group(
                filtered_prev={
                    'bed': '{root}.bed',
                    'bim': '{root}.bim',
                    'fam': '{root}.fam',
                }
            )

            # Harden the removal list by writing it to a file via to_path
            remove_list_path = f'{output_prefix}_samples_to_remove.txt'
            to_path(remove_list_path).write_text('\n'.join([f'0\t{s}' for s in samples_to_remove]))
            samples_to_remove_resource = b.read_input(remove_list_path)

            j.command(
                f"""
                set -ex
                plink --bfile {prev_resource} --allow-extra-chr --remove {samples_to_remove_resource} \\
                    --make-bed --out {j.filtered_prev}
                """
            )
            staged_prefixes.append(str(j.filtered_prev))
        else:
            staged_prefixes.append(str(prev_resource))

    # 2. Stage all new input datasets
    for _i, paths in enumerate(cohort_plink_paths):
        resource = b.read_input_group(bed=paths['bed'], bim=paths['bim'], fam=paths['fam'])
        staged_prefixes.append(str(resource))

    # 3. Define output resource group
    j.declare_resource_group(
        output_plink={
            'bed': '{root}.bed',
            'bim': '{root}.bim',
            'fam': '{root}.fam',
        }
    )

    # 4. Construct merge list and execute
    # Note: PLINK 1.9 --merge-list expects prefixes of datasets to merge
    # excluding the one specified by --bfile.
    if not staged_prefixes:
        raise ValueError('No datasets to merge')

    first_prefix = staged_prefixes[0]
    rest_prefixes = staged_prefixes[1:]

    if not rest_prefixes:
        j.command(f'plink --bfile {first_prefix} --allow-extra-chr --make-bed --out {j.output_plink}')
    else:
        merge_list_content = '\n'.join(rest_prefixes)
        j.command(
            f"""
            set -ex

            echo "{merge_list_content}" > mergelist.txt

            plink --bfile {first_prefix} --merge-list mergelist.txt --allow-extra-chr --make-bed --out {j.output_plink}
            """
        )

    # 5. Write outputs back to cloud
    b.write_output(j.output_plink, output_prefix)

    return j
