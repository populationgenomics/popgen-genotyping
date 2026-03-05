"""
Convert GTC to VCF.
"""

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

if TYPE_CHECKING:
    from hailtop.batch.job import Job


def run_gtc_to_heavy_vcf(
    gtc_path: str,
    sample_id: str,
    output_bcf_path: str,
    bpm_manifest_path: str,
    egt_cluster_path: str,
    fasta_ref_path: str,
    job_name: str = 'gtc_to_heavy_vcf',
) -> 'Job':
    """Convert a single GTC sample to VCF retaining array intensities.

    Args:
        gtc_path (str): Cloud path to the input GTC file.
        sample_id (str): Sample ID for the GTC file.
        output_bcf_path (str): Cloud path for the final output BCF file.
        bpm_manifest_path (str): Cloud path to the BPM manifest file.
        egt_cluster_path (str): Cloud path to the EGT cluster file.
        fasta_ref_path (str): Cloud path to the FASTA reference genome.
        job_name (str): Name for the Hail Batch job.

    Returns:
        Job: A Hail Batch job object.

    """
    b = get_batch()
    j = b.new_job(name=job_name)

    j.image(config_retrieve(['workflow', 'driver_image']))

    j.cpu(config_retrieve(['popgen_genotyping', 'gtc_to_heavy_vcf', 'cpu'], 2))
    j.memory(config_retrieve(['popgen_genotyping', 'gtc_to_heavy_vcf', 'memory'], 'standard'))
    j.storage(config_retrieve(['popgen_genotyping', 'gtc_to_heavy_vcf', 'storage'], '50G'))

    # Read inputs into local env
    gtc_file = b.read_input(gtc_path)
    bpm_manifest_file = b.read_input(bpm_manifest_path)
    egt_cluster_file = b.read_input(egt_cluster_path)

    # Ensure FASTA and index are in same location for bcftools
    ref_fasta = b.read_input_group(base=fasta_ref_path, fai=f'{fasta_ref_path}.fai')

    j.declare_resource_group(
        reheadered_bcf={
            'bcf': '{root}.bcf',
            'bcf.csi': '{root}.bcf.csi',
        }
    )

    j.command(
        f"""
        set -ex

        bcftools +gtc2vcf \\
            --no-version \\
            --bpm {bpm_manifest_file} \\
            --egt {egt_cluster_file} \\
            --fasta-ref {ref_fasta.base} \\
            --extra temp_vcf.tsv | \\
            {gtc_file} | \\
        bcftools sort -T ./bcftools-sort-tmp.dir | \\
        bcftools norm --no-version -c x -f {ref_fasta.base} | \\
        bcftools reheader \\
            -n {sample_id} \\
            raw.bcf \\
            -o {j.reheadered_bcf.bcf} --write-index
        """
    )

    b.write_output(j.reheadered_bcf, output_bcf_path.replace('.bcf', ''))

    return j
