"""
Stage: ExportCohortDatasets reproduction.
"""

from pathlib import Path

from scripts.test_utils import (
    BCFTOOLS_IMAGE,
    DATA_DIR,
    PLINK_IMAGE,
    run_docker,
    to_container,
)


def run_export_cohort_datasets(plink1_prefix_host: Path) -> Path:
    """
    Run the ExportCohortDatasets stage in Docker.

    Args:
        plink1_prefix_host (Path): Host path prefix for the PLINK 1.9 files.

    Returns:
        Path: Host path to the final BCF.
    """
    print('\n>>> Stage: ExportCohortDatasets (PLINK2/BCF) <<<')
    plink2_prefix: Path = DATA_DIR / 'cohort_plink2'
    plink2_int: str = to_container(plink2_prefix)
    plink1_int: str = to_container(plink1_prefix_host)
    final_bcf: Path = DATA_DIR / 'cohort.final.bcf'
    final_int: str = to_container(final_bcf)

    export_cmd: str = (
        f"bash -c 'plink2 --bfile {plink1_int} --allow-extra-chr --make-pgen --out {plink2_int} && "
        f'plink2 --pfile {plink2_int} --allow-extra-chr --export bcf --out {plink2_int} && '
        f"mv {plink2_int}.bcf {final_int}'"
    )
    run_docker(PLINK_IMAGE, export_cmd)
    run_docker(BCFTOOLS_IMAGE, f'bcftools index {final_int}')

    return final_bcf
