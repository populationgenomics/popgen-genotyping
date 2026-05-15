"""
Data preparation for genotyping pipeline reproduction.
"""

import os
import subprocess
from typing import TYPE_CHECKING

from scripts.testing_utils import (
    BCFTOOLS_IMAGE,
    DATA_DIR,
    SITES_FILE,
    run_docker,
    to_container,
)

if TYPE_CHECKING:
    from pathlib import Path


def assign_sexes(samples: list[str]) -> dict[str, str]:
    """
    Assign deterministic alternating sexes (M/F) to a sample list.

    Args:
        samples (list[str]): Sample names in cohort order.

    Returns:
        dict[str, str]: Mapping of sample name to 'M' or 'F'. Even-indexed
        samples (0, 2, 4, ...) get 'M'; odd-indexed get 'F'. With an even
        sample count this gives a 50/50 split.
    """
    return {s: ('M' if i % 2 == 0 else 'F') for i, s in enumerate(samples)}


def generate_gtcs(
    samples: list[str],
    num_snps: int,
    base_seed: int = 42,
) -> tuple[list[str], dict[str, str]]:
    """
    Generate synthetic GTC files for a list of samples.

    Args:
        samples (list[str]): List of sample names.
        num_snps (int): Number of SNPs per GTC.
        base_seed (int): Base random seed. Each sample i gets base_seed + i.

    Returns:
        tuple[list[str], dict[str, str]]: (container paths to GTCs in sample
            order, sample-name → 'M'/'F' mapping).
    """
    print(f'>>> Generating synthetic GTC files ({len(samples)} samples, base seed: {base_seed})...')
    sex_mapping: dict[str, str] = assign_sexes(samples)

    gtc_paths: list[str] = []
    for i, sample in enumerate(samples):
        gtc_path: Path = DATA_DIR / f'{sample}.gtc'
        # Always re-generate if we want to ensure unique seeds for this run
        if gtc_path.exists():
            os.remove(gtc_path)

        cmd: list[str] = [
            'python3',
            'test/scripts/generate_synthetic_gtc.py',
            str(gtc_path),
            '--num',
            str(num_snps),
            '--seed',
            str(base_seed + i),
            '--contam',
            '0.05',
            '--sex',
            sex_mapping[sample],
        ]
        subprocess.run(cmd, check=True, capture_output=True)  # noqa: S603
        gtc_paths.append(to_container(gtc_path))
    return gtc_paths, sex_mapping


def prepare_af_reference(num_snps: int) -> str:
    """
    Generate and index the population AF reference.

    Args:
        num_snps (int): Number of SNPs to include in the reference.

    Returns:
        str: Internal container path to the gzipped AF reference.
    """
    af_vcf_gz: Path = DATA_DIR / 'pop_ref.vcf.gz'
    if not af_vcf_gz.exists():
        print('\n>>> Generating population AF reference...')
        af_vcf: Path = DATA_DIR / 'pop_ref.vcf'
        gen_cmd: list[str] = ['python3', 'test/scripts/generate_af_reference.py', str(SITES_FILE), str(af_vcf)]
        subprocess.run(gen_cmd, check=True, capture_output=True)  # noqa: S603

        # Filter pop_ref to match num_snps (header + sites)
        with open(af_vcf) as f_in, open(f'{af_vcf}.tmp', 'w') as f_out:
            for i, line in enumerate(f_in):
                if i >= (num_snps + 4):  # 4 header lines
                    break
                f_out.write(line)
        os.rename(f'{af_vcf}.tmp', af_vcf)

        af_vcf_int: str = to_container(af_vcf)
        af_vcf_gz_int: str = to_container(af_vcf_gz)
        index_cmd: str = (
            f"bash -c 'bcftools view {af_vcf_int} -O z -o {af_vcf_gz_int} && bcftools index {af_vcf_gz_int}'"
        )
        run_docker(BCFTOOLS_IMAGE, index_cmd)
        if af_vcf.exists():
            os.remove(af_vcf)
    return to_container(af_vcf_gz)
