"""
Plink2 outputs from pgen/pvar
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch import Batch
    from hailtop.batch.job import BashJob


def run_plink2_qc(
    pgen_path: str,
    outputs_path: str,
    job_name: str = 'plink2_qc',
) -> BashJob:
    """
    Generate plink2 QC outputs from pgen/pvar files.

    Args:
        pgen_path (str): Cloud path to the input plink pgen, pvar file path implied at same loc.
        outputs_path (str): Cloud path for the QC output files.
        job_name (str): Name for the Hail Batch job.

    Returns:
        BashJob: The queued Hail Batch job.
    """
    b: Batch = get_batch()
    j: BashJob = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'plink2_qc'],
        image=config_retrieve(['workflow', 'plink_image']),
        default_cpu=2,
        default_storage='50G',
    )

    # Read inputs into local env. Plink2 requires 3 files in the same location
    pvar_path = pgen_path.replace('.pgen', '.pvar')
    psam_path = pgen_path.replace('.pgen', '.psam')

    # Read inputs to local env
    pgen_files = b.read_input_group(
        pgen=pgen_path,
        pvar=pvar_path,
        psam=psam_path,
    )

    # Define output files with datestamp prefix
    datestamp: str = datetime.now(tz=timezone.utc).strftime('%Y%m%d')
    output_base = f'{datestamp}_qc'

    j.declare_resource_group(
        plink_qc_outputs={
            'smiss': f'{{root}}/{output_base}.smiss',
            'vmiss': f'{{root}}/{output_base}.vmiss',
            'afreq': f'{{root}}/{output_base}.afreq',
            'hwe': f'{{root}}/{output_base}.hwe',
            'het': f'{{root}}/{output_base}.het',
            'sexcheck': f'{{root}}/{output_base}.sexcheck',
            'kin': f'{{root}}/{output_base}.kin0',
        },
    )

    check_sex_args: str = (
        'max-female-xf=0.25 min-male-xf=0.75 max-female-ycount=100 '
        "min-male-ycount=1000 cols='status,pedsex,xf,ycount,yrate'"
    )

    j.command(
        f"""
        set -ex
        plink2 --pfile {pgen_files} \\
            --missing \\
            --freq \\
            --hardy \\
            --het \\
            --check-sex {check_sex_args} \\
            --make-king-table \\
            --out {j.plink_qc_outputs}/{output_base}
        """
    )

    b.write_output(j.plink_qc_outputs, outputs_path)

    return j
