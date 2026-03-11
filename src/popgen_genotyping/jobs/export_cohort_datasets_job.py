"""
Job logic for exporting merged cohort PLINK 1.9 data to PLINK2 and BCF formats.
"""

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch.job import Job


def run_export_cohort_datasets(
    input_plink_prefix: dict[str, str],
    output_prefix: str,
    job_name: str = 'export_cohort_datasets',
) -> 'Job':
    """
    Convert PLINK 1.9 BED/BIM/FAM to PLINK2 (PGEN/PVAR/PSAM) and BCF.

    Args:
        input_plink_prefix (dict[str, str]): Dict with 'bed', 'bim', 'fam' cloud paths.
        output_prefix (str): Cloud prefix for the exported datasets.
        job_name (str): Name for the Hail Batch job.

    Returns:
        Job: A Hail Batch job object.
    """
    b = get_batch()
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'export_cohort_datasets'],
        image=config_retrieve(['workflow', 'plink_image']),
        default_cpu=4,
        default_storage='100G',
    )

    # 1. Stage inputs
    plink_input = b.read_input_group(
        bed=input_plink_prefix['bed'],
        bim=input_plink_prefix['bim'],
        fam=input_plink_prefix['fam']
    )

    # 2. Define output resource group
    j.declare_resource_group(
        output_data={
            'pgen': '{root}.pgen',
            'pvar': '{root}.pvar',
            'psam': '{root}.psam',
            'bcf': '{root}.bcf',
        }
    )

    # 3. Execute combined export command
    # Using plink2 to generate both PGEN and BCF in a single pass.
    # --allow-extra-chr is included to handle non-standard chromosomes.
    # --split-par is NOT included here because it was already handled during
    # the initial conversion from BCF to BED.
    j.command(
        f"""
        set -ex
        plink2 --bfile {plink_input} \\
            --allow-extra-chr \\
            --make-pgen \\
            --export bcf \\
            --out {j.output_data}
        """
    )

    # 4. Write outputs back to cloud
    b.write_output(j.output_data, output_prefix)

    return j
