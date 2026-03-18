"""
Combined job to convert multiple GTCs to cohort-level Heavy and Light BCFs.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def run_gtc_to_bcfs(
    gtc_paths: list[str],
    sample_mapping: dict[str, str],
    output_heavy_bcf_path: str,
    output_heavy_bcf_index_path: str,
    output_light_bcf_path: str,
    output_light_bcf_index_path: str,
    output_metadata_path: str,
    bpm_manifest_path: str,
    egt_cluster_path: str,
    fasta_ref_path: str,
    job_name: str = 'GtcToBcfs',
) -> BashJob:
    """
    Queue a job to convert multiple Illumina GTC files to a multi-sample cohort BCF.

    Args:
        gtc_paths (list[str]): List of cloud paths to the GTC files.
        sample_mapping (dict[str, str]): Mapping of old_name (barcode_pos) to new_name (SG_ID).
        output_heavy_bcf_path (str): Cloud path to save the heavy BCF.
        output_light_bcf_path (str): Cloud path to save the light BCF.
        output_metadata_path (str): Cloud path to save the metadata TSV.
        bpm_manifest_path (str): Cloud path to the BPM manifest file.
        egt_cluster_path (str): Cloud path to the EGT cluster file.
        fasta_ref_path (str): Cloud path to the FASTA reference file.
        job_name (str): Name for the Batch job.

    Returns:
        BashJob: The queued Hail Batch job.
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

    # Read reference files into the job's local storage
    bpm_file = b.read_input(bpm_manifest_path)
    egt_file = b.read_input(egt_cluster_path)
    fasta_file = b.read_input_group(
        base=fasta_ref_path,
        fai=fasta_ref_path + '.fai',
    )

    # Read all GTC files
    gtc_files = [b.read_input(p) for p in gtc_paths]
    gtc_arg = ' '.join([str(f) for f in gtc_files])

    # Create reheader mapping file content
    mapping_content = '\n'.join([f'{old} {new}' for old, new in sample_mapping.items()])

    # Index outputs
    j.declare_resource_group(
        heavy_bcf={'heavy_bcf': '{root}.bcf', 'heavy_bcf_index': '{root}_index.csi'},
        light_bcf={'light_bcf': '{root}.bcf', 'light_bcf_index': '{root}_index.csi'},
        metadata_tsv={'metadata_tsv': '{root}.tsv'},
    )

    # Building the command
    j.command(
        f"""
        set -ex

        mkdir -p $BATCH_TMPDIR/bcftools-tmp

        # Create reheader mapping file inside the job
        cat <<EOF > reheader_map.txt
{mapping_content}
EOF

        bcftools +gtc2vcf \\
            --no-version \\
            --do-not-check-bpm \\
            --bpm {bpm_file} \\
            --egt {egt_file} \\
            --fasta-ref {fasta_file.base} \\
            --extra metadata_raw.tsv \\
            {gtc_arg} | \\
        bcftools norm -m -both --no-version -c x -f {fasta_file.base} | \\
        bcftools sort -T $BATCH_TMPDIR/bcftools-tmp | \\
        bcftools reheader -s reheader_map.txt | \\
        bcftools view -O b -o {j.heavy_bcf} --write-index=csi

        bcftools annotate --no-version -x ^FORMAT/GT,FORMAT/GQ {j.heavy_bcf} \\
        -O b -o {j.light_bcf} --write-index=csi

        mv metadata_raw.tsv {j.metadata_tsv}
        """
    )

    b.write_output(j.heavy_bcf, output_heavy_bcf_path)
    b.write_output(j.heavy_bcf_index, output_heavy_bcf_index_path)
    b.write_output(j.light_bcf, output_light_bcf_path)
    b.write_output(j.light_bcf_index, output_light_bcf_index_path)
    b.write_output(j.metadata_tsv, output_metadata_path)

    return j
