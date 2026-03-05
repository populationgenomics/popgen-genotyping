"""
Convert Heavy BCF to Light BCF by stripping FORMAT fields.
"""

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

if TYPE_CHECKING:
    from hailtop.batch.job import Job


def run_heavy_to_light_vcf(
    input_bcf_path: str,
    id_mappings_path: str,
    output_bcf_path: str,
    fasta_ref_path: str,
    job_name: str = 'heavy_to_light_vcf',
) -> 'Job':
    """Convert a heavy BCF sample to a light BCF.

    Args:
        input_bcf_path (str): Cloud path to the input heavy BCF file.
        id_mappings_path (str): Cloud path to the id_mappings.txt file for reheadering.
        output_bcf_path (str): Cloud path for the final light BCF file.
        fasta_ref_path (str): Cloud path to the FASTA reference genome.
        job_name (str): Name for the Hail Batch job.

    Returns:
        Job: A Hail Batch job object.

    """
    b = get_batch()
    j = b.new_job(name=job_name)

    j.image(config_retrieve(['workflow', 'driver_image']))

    # Use smaller resources as this is a lighter job
    j.cpu(config_retrieve(['popgen_genotyping', 'heavy_to_light_vcf', 'cpu'], 1))
    j.memory(config_retrieve(['popgen_genotyping', 'heavy_to_light_vcf', 'memory'], 'standard'))
    j.storage(config_retrieve(['popgen_genotyping', 'heavy_to_light_vcf', 'storage'], '10G'))

    # Read inputs
    input_bcf = b.read_input_group(
        bcf=input_bcf_path,
        csi=f'{input_bcf_path}.csi'
    )
    id_map = b.read_input(id_mappings_path)
    ref_fasta = b.read_input_group(base=fasta_ref_path, fai=f'{fasta_ref_path}.fai')

    j.declare_resource_group(
        light_bcf={
            'bcf': '{root}.bcf',
            'bcf.csi': '{root}.bcf.csi',
        }
    )

    # Optimized command:
    # 1. Annotate to strip FORMAT fields (except GT, GQ)
    # 2. Sort using /dev/shm for speed
    # 3. Norm with reference
    # 4. Reheader with sample ID and write index
    j.command(
        f"""
        set -ex

        mkdir -p /dev/shm/bcftools-tmp

        bcftools annotate -x ^FORMAT/GT,FORMAT/GQ -Ou {input_bcf.bcf} | \\
        bcftools sort -Ou -T /dev/shm/bcftools-tmp | \\
        bcftools norm --no-version -c x -f {ref_fasta.base} -Ou | \\
        bcftools reheader -s {id_map} -o {j.light_bcf.bcf} --write-index

        # Ensure index has the correct extension for Hail Batch to recognize it
        # bcftools reheader --write-index creates .bcf.csi
        """
    )

    b.write_output(j.light_bcf, output_bcf_path.replace('.bcf', ''))

    return j
