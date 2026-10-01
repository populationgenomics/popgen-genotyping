"""
Combined job to convert multiple GTCs to cohort-level Heavy and Light BCFs.
"""

from __future__ import annotations

from collections import Counter
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
    output_light_bcf_path: str,
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

    Raises:
        ValueError: If sample_mapping is empty or maps two samples to the same SG ID.
    """
    if not sample_mapping:
        raise ValueError('sample_mapping is empty: no expected SG IDs to reheader the cohort BCF to')
    sg_id_counts: Counter[str] = Counter(sample_mapping.values())
    duplicated: list[str] = sorted(sg_id for sg_id, n in sg_id_counts.items() if n > 1)
    if duplicated:
        raise ValueError(f'sample_mapping maps several names to the same SG ID, duplicated: {duplicated}')

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
    # Expected cohort sample set, asserted against the reheadered BCF below
    expected_samples_content = '\n'.join(sorted(sg_id_counts))

    # Outputs
    j.declare_resource_group(
        heavy_bcf={'bcf': '{root}.bcf', 'bcf.csi': '{root}.bcf.csi'},
        light_bcf={'bcf': '{root}.bcf', 'bcf.csi': '{root}.bcf.csi'},
        metadata_tsv={'tsv': '{root}.tsv'},
    )

    # Building the command
    j.command(
        f"""
        set -exo pipefail

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
        bcftools view -O b -o {j.heavy_bcf.bcf} --write-index=csi

        # bcftools reheader silently keeps any sample name it cannot match (e.g. a manifest barcode
        # mangled to scientific notation), so assert the BCF holds exactly the cohort's SG IDs here
        # instead of letting raw barcode IDs surface only in a much later stage.
        cat <<EOF > expected_samples.txt
{expected_samples_content}
EOF
        LC_ALL=C sort expected_samples.txt -o expected_samples.txt
        bcftools query -l {j.heavy_bcf.bcf} | LC_ALL=C sort > actual_samples.txt
        if ! diff expected_samples.txt actual_samples.txt > sample_diff.txt; then
            echo "ERROR: heavy BCF samples are not exactly the expected SG IDs:" \\
                "expected $(wc -l < expected_samples.txt | tr -d ' ')," \\
                "found $(wc -l < actual_samples.txt | tr -d ' ')," \\
                "$(grep -c '^[<>]' sample_diff.txt) differing lines" \\
                "('<' = expected SG ID missing from BCF, '>' = unexpected name in BCF; first 10 shown)" >&2
            head -n 10 sample_diff.txt >&2
            exit 1
        fi

        bcftools annotate --no-version -x ^FORMAT/GT,FORMAT/GQ {j.heavy_bcf.bcf} \\
        -O b -o {j.light_bcf.bcf} --write-index=csi

        mv metadata_raw.tsv {j.metadata_tsv.tsv}
        """
    )

    b.write_output(j.heavy_bcf, str(output_heavy_bcf_path).replace('.bcf', ''))
    b.write_output(j.light_bcf, str(output_light_bcf_path).replace('.bcf', ''))
    b.write_output(j.metadata_tsv, str(output_metadata_path).replace('.tsv', ''))

    return j
