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
    from hailtop.batch.resource import ResourceGroup


def run_merge_plink(
    cohort_plink_paths: list[dict[str, str]],
    output_prefix: str,
    previous_aggregate_resource: ResourceGroup | None = None,
    samples_to_remove: list[str] | None = None,
    keep_samples: list[str] | None = None,
    job_name: str = 'merge_cohort_plink',
) -> BashJob:
    """
    Merge multiple PLINK 1.9 datasets into a single unified dataset, with rolling aggregate support.

    Args:
        cohort_plink_paths (list[dict[str, str]]): List of dicts, each with 'bed', 'bim', 'fam' cloud paths.
        output_prefix (str): Cloud prefix for the final merged PLINK 1.9 files.
        previous_aggregate_resource (ResourceGroup, optional): Resource group for the previous rolling aggregate.
        samples_to_remove (list[str], optional): List of SG IDs to remove from the previous aggregate.
        keep_samples (list[str], optional): SG IDs to retain in the final fileset. Whole plate filesets are
            merged (a plate may carry withdrawn/excluded SGs), so a final ``plink --keep`` pass trims the
            merged result to exactly this membership (the super cohort). When None (the default) no trim is
            applied and the merged fileset is written as-is — the pre-two-phase behaviour.
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
    if previous_aggregate_resource:
        prev_resource = previous_aggregate_resource

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

            # --keep-allele-order: PLINK 1.9 otherwise resets A1 to the minor allele,
            # diverging from the A1=ALT/A2=REF orientation that every other step
            # (and any future reference-panel merge) relies on.
            j.command(
                f"""
                set -ex
                plink --bfile {prev_resource} --allow-extra-chr --output-chr chrM \\
                    --remove {samples_to_remove_resource} \\
                    --keep-allele-order --make-bed --out {j.filtered_prev}
                """
            )
            staged_prefixes.append(str(j.filtered_prev))
        else:
            staged_prefixes.append(str(prev_resource))

    # 2. Stage all new input datasets
    for _i, paths in enumerate(cohort_plink_paths):
        resource = b.read_input_group(bed=paths['bed'], bim=paths['bim'], fam=paths['fam'])
        staged_prefixes.append(str(resource))

    # 3. Define output resource groups. When keep_samples is set the merge lands in an
    #    intermediate fileset and a final --keep pass (step 5) trims it to super-cohort
    #    membership; otherwise the merge writes straight to the final output.
    j.declare_resource_group(
        output_plink={
            'bed': '{root}.bed',
            'bim': '{root}.bim',
            'fam': '{root}.fam',
        }
    )
    if keep_samples:
        j.declare_resource_group(
            merged_untrimmed={
                'bed': '{root}.bed',
                'bim': '{root}.bim',
                'fam': '{root}.fam',
            }
        )
        merge_target = j.merged_untrimmed
    else:
        merge_target = j.output_plink

    # 4. Construct merge list and execute
    # Note: PLINK 1.9 --merge-list expects prefixes of datasets to merge
    # excluding the one specified by --bfile.
    if not staged_prefixes:
        raise ValueError('No datasets to merge')

    first_prefix = staged_prefixes[0]
    rest_prefixes = staged_prefixes[1:]

    if not rest_prefixes:
        j.command(
            f"""
            plink --bfile {first_prefix} \
                --allow-extra-chr \
                --output-chr chrM \
                --make-bed \
                --keep-allele-order \
                --out {merge_target}
            """
        )
    else:
        merge_list_content = '\n'.join(rest_prefixes)
        j.command(
            f"""
            set -ex

            echo "{merge_list_content}" > mergelist.txt

            plink --bfile {first_prefix} \
                --merge-list mergelist.txt \
                --allow-extra-chr \
                --output-chr chrM \
                --make-bed \
                --keep-allele-order \
                --out {merge_target}
            """
        )

    # 5. Trim the merged fileset to the super-cohort membership. Whole plate filesets are
    #    merged (a plate may carry withdrawn/excluded SGs), so this final --keep is what
    #    guarantees merged ⊆ super cohort. Uses the same FID=0 / IID convention as the
    #    --remove list above.
    if keep_samples:
        keep_list_path = f'{output_prefix}_samples_to_keep.txt'
        to_path(keep_list_path).write_text('\n'.join([f'0\t{s}' for s in keep_samples]))
        keep_samples_resource = b.read_input(keep_list_path)

        # --keep-allele-order: preserve the A1=ALT/A2=REF orientation through the trim
        # (PLINK 1.9 otherwise resets A1 to the minor allele).
        j.command(
            f"""
            set -ex
            plink --bfile {merge_target} --allow-extra-chr --output-chr chrM \\
                --keep {keep_samples_resource} \\
                --keep-allele-order --make-bed --out {j.output_plink}

            kept=$(wc -l < {j.output_plink.fam})
            if [ "$kept" -ne {len(keep_samples)} ]; then
                echo "expected {len(keep_samples)} samples after --keep, got $kept" >&2
                exit 1
            fi
            """
        )

    # 6. Write outputs back to cloud
    b.write_output(j.output_plink, output_prefix)

    return j
