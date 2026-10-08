"""
Job logic for selecting the samples left out of every relatedness estimator.
"""

from __future__ import annotations

from importlib.resources import files
from typing import TYPE_CHECKING

from cpg_utils.config import get_driver_image
from cpg_utils.hail_batch import get_batch

from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def run_relatedness_exclusions(
    bafregress_paths: list[str],
    smiss_path: str,
    contamination_max: float,
    fmiss_max: float,
    output_excluded_path: str,
    output_remove_list_path: str,
    job_name: str = 'relatedness_exclusions',
) -> BashJob:
    """
    Select samples to exclude from relatedness estimation for contamination or missingness.

    Runs ``scripts/relatedness_exclusions.py`` against the BAFRegress output(s) and the
    merged-cohort ``.smiss`` to select samples whose BAFRegress estimate exceeds
    ``contamination_max`` or whose per-sample ``F_MISS`` exceeds ``fmiss_max``. Either condition
    generates thousands of spurious pairwise relationships under KING, KING-robust and IBIS,
    which would otherwise contaminate every clean sample's relatedness list. Writes an audit
    TSV (``IID``, ``REASON``, ``BAF_REGRESS``, ``F_MISS``; header-only when nothing is
    excluded) and a plink2 ``--remove`` sample-ID list (``#IID`` header, one IID per line;
    header-only when nothing is excluded, which removes nothing).

    Args:
        bafregress_paths: Cloud paths to the cohort's BAFRegress output file(s).
        smiss_path: Cloud path to the merged-cohort PLINK2 ``.smiss``.
        contamination_max: BAFRegress exclusion threshold (strictly greater than).
        fmiss_max: Per-sample ``F_MISS`` exclusion threshold (strictly greater than).
        output_excluded_path: Cloud path for the excluded-samples audit TSV.
        output_remove_list_path: Cloud path for the plink2 ``--remove`` list.
        job_name: Name for the Hail Batch job.

    Returns:
        The queued Hail Batch job.
    """
    b = get_batch()
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'relatedness_exclusions'],
        image=get_driver_image(),
        default_cpu=1,
        default_storage='10G',
    )

    script_path = str(files('popgen_genotyping.scripts').joinpath('relatedness_exclusions.py'))
    script = b.read_input(script_path)
    bafregress_files = [b.read_input(path) for path in bafregress_paths]
    bafregress_bash_args = ' '.join(str(f) for f in bafregress_files)
    smiss_file = b.read_input(smiss_path)

    j.command(
        f"""
        set -ex
        python3 {script} \\
            --bafregress {bafregress_bash_args} \\
            --smiss {smiss_file} \\
            --contamination-max {contamination_max} \\
            --fmiss-max {fmiss_max} \\
            --output-tsv {j.excluded_tsv} \\
            --output-remove-list {j.remove_list}
        """,
    )

    b.write_output(j.excluded_tsv, output_excluded_path)
    b.write_output(j.remove_list, output_remove_list_path)

    return j
