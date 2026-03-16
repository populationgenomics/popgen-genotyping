"""
Stage: BafRegress reproduction.
"""

from pathlib import Path

from scripts.testing_utils import (
    BCFTOOLS_IMAGE,
    DATA_DIR,
    run_docker,
    to_container,
)


def run_bafregress(heavy_bcf_host: Path, af_ref_int: str | None = None) -> Path:
    """
    Run the BafRegress stage in Docker.

    Args:
        heavy_bcf_host (Path): Host path to the heavy BCF.
        af_ref_int (str, optional): Container path to the population AF reference.
            If None, BAFRegress is run without an external AF reference, and AF
            is estimated from the input BCF using fill-tags.

    Returns:
        Path: Host path to the BafRegress output.
    """
    print('\n>>> Stage: BafRegress (Cohort Level) <<<')
    baf_out_host: Path = DATA_DIR / 'cohort.BAFRegress.txt'
    baf_out_int: str = to_container(baf_out_host)
    heavy_int: str = to_container(heavy_bcf_host)

    if af_ref_int:
        af_arg: str = f'-a {af_ref_int}'
        input_bcf: str = heavy_int
        baf_cmd: str = f"bash -c 'bcftools +BAFregress {af_arg} --tag AF {input_bcf} > {baf_out_int}'"
    else:
        # Estimate AF from the cohort itself if no reference is provided
        print('No AF reference provided. Estimating AF from the cohort using fill-tags...')
        baf_cmd: str = (
            f"bash -c 'bcftools +fill-tags {heavy_int} -- -t AF | bcftools +BAFregress --tag AF > {baf_out_int}'"
        )

    run_docker(BCFTOOLS_IMAGE, baf_cmd)

    return baf_out_host
