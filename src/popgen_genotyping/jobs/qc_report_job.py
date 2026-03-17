"""Job to merge the QC data files into a CSV."""

from __future__ import annotations

from importlib.resources import files
from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def run_qc_report(
    plink_qc_prefix: str,
    bafregress_paths: list[str],
    output_path: str,
    job_name: str = 'qc_report',
) -> BashJob:
    """Merge the QC files into a single CSV QC summary file.

    Args:
        plink_qc_prefix: Cloud path prefix for PLINK2 QC files.
        bafregress_paths: Cloud paths to bafregress files for merging.
        output_path: Cloud path to output QC summary CSV.
        job_name: Name for the Hail Batch Job.

    Returns:
        A Hail Batch BashJob object.
    """
    b = get_batch()
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'merge_qc'],
        image=config_retrieve(['workflow', 'plink_image']),
        default_cpu=2,
        default_storage='20G',
    )

    # Stage the merge_qc script into the container
    script_path = str(files('popgen_genotyping.scripts').joinpath('merge_qc.py'))
    script = b.read_input(script_path)

    # Read in PLINK2 QC files
    sexcheck_file = b.read_input(f'{plink_qc_prefix}.sexcheck')
    het_file = b.read_input(f'{plink_qc_prefix}.het')
    smiss_file = b.read_input(f'{plink_qc_prefix}.smiss')
    kin0_file = b.read_input(f'{plink_qc_prefix}.kin0')

    # Read in bafregress files
    bafregress_files = [b.read_input(path) for path in bafregress_paths]
    bafregress_bash_args = ' '.join([str(f) for f in bafregress_files])

    # Define output file
    j.declare_resource_group(output_csv={'qc_report': '{root}.csv'})

    # Build the bafregress flag only if files are provided
    bafregress_flag = f'--bafregress {bafregress_bash_args}' if bafregress_files else ''

    j.command(
        f"""\
        set -ex
        python3 {script} \
            --sexcheck {sexcheck_file} \
            --het {het_file} \
            --smiss {smiss_file} \
            --kin0 {kin0_file} \
            --output {j.output_csv.qc_report} \
            {bafregress_flag}
        """,
    )

    b.write_output(j.output_csv.qc_report, output_path)

    return j
