from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

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
    j = b.new_job(name=job_name)

    j.image(config_retrieve(['workflow', 'driver_image']))

    j.cpu(config_retrieve(['popgen_genotyping', 'bafregress', 'cpu'], 1))
    j.memory(config_retrieve(['popgen_genotyping', 'bafregress', 'memory'], 'standard'))
    j.storage(config_retrieve(['popgen_genotyping', 'bafregress', 'storage'], '10G'))

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
