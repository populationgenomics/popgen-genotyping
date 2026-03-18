"""
Job logic for estimating sample contamination using BAFRegress.
"""

from __future__ import annotations

from hailtop.batch.job import BashJob
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from popgen_genotyping.utils import register_job


def run_bafregress(
    bcf_path: str,
    output_path: str,
    af_ref_path: str | None = None,
    job_name: str = 'bafregress',
) -> BashJob:
    """
    Run BAFRegress on a BCF file to estimate sample contamination.

    Args:
        bcf_path (str): Cloud path to input BCF file.
        output_path (str): Cloud path to output BAFRegress.txt file.
        af_ref_path (str, optional): Cloud path to population AF reference VCF.
            If None, BAFRegress is run without an external AF reference, and AF
            is estimated from the input BCF using fill-tags.
        job_name (str): Name for the Hail Batch job. Defaults to 'bafregress'.

    Returns:
        BashJob: A Hail Batch BashJob object.
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
    assert isinstance(j, BashJob)

    # Read the input BCF file with index.
    bcf_file = b.read_input_group(bcf=bcf_path, csi=f'{bcf_path}.csi')

    # Explicitly define the output resource to avoid dynamic attribute confusion
    j.declare_resource_group(output={'txt': '{root}.txt'})

    if af_ref_path:
        af_ref = b.read_input_group(vcf=af_ref_path, tbi=f'{af_ref_path}.tbi')
        j.command(
            f"""
            set -ex

            # Run BAFRegress and redirect stdout directly to the Hail resource
            bcftools +BAFregress \\
                -a "{af_ref.vcf}" \\
                --tag AF \\
                "{bcf_file.bcf}" > {j.output.txt}
            """
        )
    else:
        # Fallback to cohort-level AF estimation
        j.command(
            f"""
            set -ex

            # Estimate AF from the cohort itself and pipe to BAFRegress
            bcftools +fill-tags "{bcf_file.bcf}" -- -t AF | \\
                bcftools +BAFregress \\
                --tag AF > {j.output.txt}
            """
        )

    # Write the resulting text file to the cloud bucket
    b.write_output(j.output.txt, output_path)

    return j
