"""
Combined job to convert GTC to both Heavy and Light BCFs.
"""

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
...
def run_gtc_to_bcfs(
    gtc_path: str,
    sample_id: str,
    output_heavy_bcf_path: str,
    output_light_bcf_path: str,
    output_metadata_path: str,
    bpm_manifest_path: str,
    egt_cluster_path: str,
    fasta_ref_path: str,
    job_name: str = 'gtc_to_bcfs',
) -> 'Job':
    """Convert a single GTC sample to both Heavy (with intensities) and Light BCFs.

    Args:
        gtc_path (str): Cloud path to the input GTC file.
        sample_id (str): Sample ID for the GTC file.
        output_heavy_bcf_path (str): Cloud path for the heavy BCF output.
        output_light_bcf_path (str): Cloud path for the light BCF output.
        output_metadata_path (str): Cloud path for the GTC metadata TSV output.
        bpm_manifest_path (str): Cloud path to the BPM manifest file.
        egt_cluster_path (str): Cloud path to the EGT cluster file.
        fasta_ref_path (str): Cloud path to the FASTA reference genome.
        job_name (str): Name for the Hail Batch job.

    Returns:
        Job: A Hail Batch job object.
    """
    b = get_batch()
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'gtc_to_bcfs'],
        image=config_retrieve(['workflow', 'bcftools_image']),
        default_cpu=2,
        default_storage='50G',
    )

    # Read inputs


    # Read inputs
    gtc_file = b.read_input(gtc_path)
    bpm_manifest_file = b.read_input(bpm_manifest_path)
    egt_cluster_file = b.read_input(egt_cluster_path)
    ref_fasta = b.read_input_group(base=fasta_ref_path, fai=f'{fasta_ref_path}.fai')

    j.declare_resource_group(
        heavy_bcf={
            'bcf': '{root}.bcf',
            'bcf.csi': '{root}.bcf.csi',
        },
        light_bcf={
            'bcf': '{root}.bcf',
            'bcf.csi': '{root}.bcf.csi',
        },
    )

    # Optimized command:
    # 1. gtc2vcf -> sort (in /dev/shm) -> norm -> reheader -> Heavy BCF
    # 2. Heavy BCF -> annotate (strip) -> Light BCF
    # Note: reheadering is done only once at the source of truth for the sample ID.
    metadata_filename = f'{sample_id}_gtc_metadata.tsv'
    j.command(
        f"""
        set -ex

        mkdir -p /dev/shm/bcftools-tmp

        bcftools +gtc2vcf \\
            --no-version \\
            --bpm {bpm_manifest_file} \\
            --egt {egt_cluster_file} \\
            --fasta-ref {ref_fasta.base} \\
            --extra {metadata_filename} \\
            {gtc_file} | \\
        bcftools norm -m -both --no-version -c x -f {ref_fasta.base} | \\
        bcftools sort -T /dev/shm/bcftools-tmp | \\
        bcftools reheader -n {sample_id} | \\
        bcftools view -O b -o {j.heavy_bcf.bcf} --write-index

        bcftools annotate --no-version -x ^FORMAT/GT,FORMAT/GQ {j.heavy_bcf.bcf} \\
        -O b -o {j.light_bcf.bcf} --write-index


        # Move metadata to a resource or just write it directly if it's small
        mv {metadata_filename} {j.metadata_tsv}
        """
    )

    b.write_output(j.heavy_bcf, output_heavy_bcf_path.replace('.bcf', ''))
    b.write_output(j.light_bcf, output_light_bcf_path.replace('.bcf', ''))
    b.write_output(j.metadata_tsv, output_metadata_path)

    return j
