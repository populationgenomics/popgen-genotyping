"""
Per-cohort PLINK2 QC: --missing / --freq / --hardy on the per-cohort PLINK 1.9 fileset.

Sibling of `Plink2Qc` (which runs on the merged PGEN). Outputs land in default
storage so per-batch metrics persist beyond the cohort's tmp intermediates.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch import Batch
    from hailtop.batch.job import BashJob


def run_cohort_plink_qc(
    bed_path: str,
    output_prefix: str,
    job_name: str = 'cohort_plink_qc',
) -> BashJob:
    """
    Run PLINK2 site-level QC on the per-cohort PLINK 1.9 binary fileset.

    Args:
        bed_path: Cloud path to the cohort `.bed`; `.bim` / `.fam` are assumed
            to sit alongside with matching basename.
        output_prefix: Cloud path prefix (without extension) for the QC outputs.
            Files `{prefix}.vmiss`, `.afreq`, `.hardy`, `.log` are written.
        job_name: Hail Batch job display name.

    Returns:
        BashJob with a declared resource group that includes all four outputs.
    """
    b: Batch = get_batch()
    j: BashJob = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'cohort_plink_qc'],
        image=config_retrieve(['workflow', 'plink_image']),
        default_cpu=2,
        default_storage='20G',
    )

    bim_path = bed_path.replace('.bed', '.bim')
    fam_path = bed_path.replace('.bed', '.fam')

    plink_files = b.read_input_group(
        bed=bed_path,
        bim=bim_path,
        fam=fam_path,
    )

    j.declare_resource_group(
        plink_qc_outputs={
            'vmiss': '{root}.vmiss',
            'afreq': '{root}.afreq',
            'hardy': '{root}.hardy',
            'log': '{root}.log',
        },
    )

    j.command(
        f"""
        set -ex
        plink2 --bfile {plink_files} \\
            --missing variant-only \\
            --freq \\
            --hardy \\
            --out {j.plink_qc_outputs}
        """,
    )

    b.write_output(j.plink_qc_outputs, output_prefix)

    return j
