"""
Job logic for converting a cohort-level multi-sample BCF into a PLINK 1.9 dataset.
"""

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def run_cohort_bcf_to_plink(
    bcf_path: str,
    output_prefix: str,
    sex_mapping: dict[str, str] | None = None,
    job_name: str = 'cohort_bcf_to_plink',
) -> 'BashJob':
    """
    Directly convert a multi-sample cohort BCF to PLINK 1.9 format using PLINK2.

    Args:
        bcf_path (str): Cloud path to the cohort-level multi-sample BCF.
        output_prefix (str): Cloud prefix for the output PLINK 1.9 files.
        sex_mapping (dict[str, str], optional): Mapping of SG ID to sex code (1 or 2).
        job_name (str): Name for the Hail Batch job.

    Returns:
        BashJob: A Hail Batch job object.
    """
    b = get_batch()
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'cohort_bcf_to_plink'],
        image=config_retrieve(['workflow', 'plink_image']),
        default_cpu=4,
        default_storage='50G',
    )

    # 1. Stage the cohort BCF and index
    bcf_file = b.read_input_group(bcf=bcf_path, csi=f'{bcf_path}.csi')

    # 2. Define output resource group
    j.declare_resource_group(
        output_plink={
            'bed': '{root}.bed',
            'bim': '{root}.bim',
            'fam': '{root}.fam',
        }
    )

    # 3. Handle sex metadata if provided
    update_sex_cmd = ''
    if sex_mapping:
        sex_tsv_lines: list[str] = []
        for sg_id, sex_code in sex_mapping.items():
            # PLINK format: FamilyID(0) SampleID Sex
            sex_tsv_lines.append(f'0\t{sg_id}\t{sex_code}')

        sex_tsv_content = '\\n'.join(sex_tsv_lines)
        j.command(f'echo -e "{sex_tsv_content}" > sex_metadata.tsv')
        update_sex_cmd = '--update-sex sex_metadata.tsv'

    # 4. Direct conversion using PLINK2
    j.command(
        f"""
        set -ex

        plink2 \\
            --bcf {bcf_file.bcf} \\
            --max-alleles 2 \\
            --split-par hg38 \\
            --set-all-var-ids '@:#:$r:$a' \\
            --allow-extra-chr \\
            {update_sex_cmd} \\
            --make-bed \\
            --out {j.output_plink}
        """
    )

    # 5. Write outputs back to cloud
    b.write_output(j.output_plink, output_prefix)

    return j
