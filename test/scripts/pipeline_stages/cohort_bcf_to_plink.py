"""
Stage: CohortBcfToPlink reproduction.
"""

from pathlib import Path

from scripts.testing_utils import (
    DATA_DIR,
    PLINK_IMAGE,
    run_docker,
    to_container,
)


def run_cohort_bcf_to_plink(samples: list[str], light_bcf_host: Path) -> Path:
    """
    Run the CohortBcfToPlink stage in Docker.

    Args:
        samples (list[str]): List of sample names.
        light_bcf_host (Path): Host path to the light BCF.

    Returns:
        Path: Host path prefix for the generated PLINK 1.9 files.
    """
    print('\n>>> Stage: CohortBcfToPlink (PLINK 1.9) <<<')
    plink1_prefix: Path = DATA_DIR / 'cohort_plink1'
    plink1_int: str = to_container(plink1_prefix)
    light_int: str = to_container(light_bcf_host)

    sex_mapping_file: Path = DATA_DIR / 'sex_mapping.txt'
    with open(sex_mapping_file, 'w') as f:
        for s in samples:
            f.write(f'0 {s} 1\n')
    sex_int: str = to_container(sex_mapping_file)

    plink1_cmd: str = (
        f'bash -c \'plink2 --bcf {light_int} --max-alleles 2 --split-par hg38 --set-all-var-ids "@:#:\\$r:\\$a" '
        f"--update-sex {sex_int} --allow-extra-chr --make-bed --out {plink1_int}'"
    )
    run_docker(PLINK_IMAGE, plink1_cmd)

    return plink1_prefix
