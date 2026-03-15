"""
Stage: BafRegress reproduction.
"""

from pathlib import Path

from scripts.test_utils import (
    BCFTOOLS_IMAGE,
    DATA_DIR,
    run_docker,
    to_container,
)


def run_bafregress(heavy_bcf_host: Path, af_ref_int: str) -> Path:
    """
    Run the BafRegress stage in Docker.

    Args:
        heavy_bcf_host (Path): Host path to the heavy BCF.
        af_ref_int (str): Container path to the population AF reference.

    Returns:
        Path: Host path to the BafRegress output.
    """
    print('\n>>> Stage: BafRegress (Cohort Level) <<<')
    baf_out_host: Path = DATA_DIR / 'cohort.BAFRegress.txt'
    baf_out_int: str = to_container(baf_out_host)
    heavy_int: str = to_container(heavy_bcf_host)

    baf_cmd: str = f"bash -c 'bcftools +BAFregress -a {af_ref_int} --tag AF {heavy_int} > {baf_out_int}'"
    run_docker(BCFTOOLS_IMAGE, baf_cmd)

    return baf_out_host
