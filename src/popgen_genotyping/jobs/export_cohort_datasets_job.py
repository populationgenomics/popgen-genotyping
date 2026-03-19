"""
Job logic for exporting merged cohort PLINK 1.9 data to PLINK2 and BCF formats.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def run_export_cohort_datasets(
    input_plink_prefix: dict[str, str],
    output_prefix: str,
    bcf_output_path: str,
    job_name: str = 'export_cohort_datasets',
) -> BashJob:
    """
    Convert PLINK 1.9 BED/BIM/FAM to PLINK2 (PGEN/PVAR/PSAM) and BCF.

    Args:
        input_plink_prefix (dict[str, str]): Dict with 'bed', 'bim', 'fam' cloud paths.
        output_prefix (str): Cloud prefix for the PLINK2 (pgen/pvar/psam) outputs.
        bcf_output_path (str): Cloud path for the BCF output (written separately to tmp storage).
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
        bed=input_plink_prefix['bed'], bim=input_plink_prefix['bim'], fam=input_plink_prefix['fam']
    )

    # 2. Define output resource groups with separate roots
    j.declare_resource_group(
        plink2_output={
            'pgen': '{root}.pgen',
            'pvar': '{root}.pvar',
            'psam': '{root}.psam',
        }
    )
    j.declare_resource_group(
        bcf_output={
            'bcf': '{root}.bcf',
        }
    )

    # 3. Execute combined export command
    # Using plink2 to generate both PGEN and BCF in a single pass.
    # --allow-extra-chr is included to handle non-standard chromosomes.
    # --split-par is NOT included here because it was already handled during
    # the initial conversion from BCF to BED.
    # Both resource groups share the same --out prefix so plink2 writes all files together.
    j.command(
        f"""
        set -ex
        plink2 --bfile {plink_input} \\
            --allow-extra-chr \\
            --make-pgen \\
            --export bcf ref-first \\
            --out {j.plink2_output}
        cp {j.plink2_output}.bcf {j.bcf_output.bcf}
        """
    )

    # 4. Write outputs back to cloud
    # PLINK2 files go to main storage, BCF goes to tmp storage separately
    b.write_output(j.plink2_output, output_prefix)
    b.write_output(j.bcf_output.bcf, bcf_output_path)

    return j
