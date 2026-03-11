from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch.job import Job


def run_bafregress(
    bcf_path: str,
    output_path: str,
    job_name: str = 'bafregress',
) -> 'Job':
    """Run BAFRegress on a BCF file to estimate sample contamination.

    Args:
        bcf_path (str): Cloud path to input BCF file.
        output_path (str): Cloud path to output BAFRegress.txt file.

    Returns:
        Job: A Hail Batch job object.

    """
    b = get_batch()
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'bafregress'],
        image=config_retrieve(['workflow', 'bcftools_image']),
        default_cpu=1,
        default_storage='10G',
    )

    # Read the input BCF file with index.
    bcf_file = b.read_input_group(bcf=bcf_path, csi=f'{bcf_path}.csi')

    j.command(
        f"""
        set -ex

        # Run BAFRegress and redirect stdout directly to the Hail resource
        bcftools +BAFregress \\
            "{bcf_file.bcf}" > {j.baf_regress_out}
        """
    )

    # Write the resulting text file to the cloud bucket
    b.write_output(j.baf_regress_out, output_path)

    return j
