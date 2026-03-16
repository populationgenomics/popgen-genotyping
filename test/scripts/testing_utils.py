"""
Shared utilities for genotyping pipeline reproduction and testing.
"""

import subprocess
from pathlib import Path

# --- Configuration ---
BCFTOOLS_IMAGE: str = 'australia-southeast1-docker.pkg.dev/cpg-common/images/bcftools:1.23-1'
PLINK_IMAGE: str = 'australia-southeast1-docker.pkg.dev/cpg-common/images/plink:1.9-20250819-PLINK-2.0-20260228-1'

ROOT: Path = Path(__file__).parent.parent.parent
LOCAL_DIR: Path = ROOT / 'test' / 'local'
DATA_DIR: Path = LOCAL_DIR / 'reproduce_cohort'
SITES_FILE: Path = ROOT / 'test' / 'data' / 'sites.txt.gz'

# Reference files on host
BPM_HOST: Path = LOCAL_DIR / 'reference' / 'GDA-8v1-0_D2.bpm'
EGT_HOST: Path = LOCAL_DIR / 'reference' / 'GDA-8v1-0_D1_ClusterFile.egt'
FASTA_HOST: Path = LOCAL_DIR / 'reference' / 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'


def to_container(path: Path) -> str:
    """
    Convert a host Path to an internal container path string.

    Args:
        path (Path): The host path.

    Returns:
        str: The internal container path (prefixed with /data).
    """
    return f'/data/{path.relative_to(ROOT)}'


def run_docker(image: str, command: str, entrypoint: str | None = None, shm_size: str = '8g') -> None:
    """
    Run a command inside a Docker container with standard mounting.

    Args:
        image (str): Docker image to run.
        command (str): Command to execute inside the container.
        entrypoint (str, optional): Custom entrypoint for the container.
        shm_size (str): Shared memory size. Defaults to '8g'.
    """
    entrypoint_arg: str = f'--entrypoint {entrypoint}' if entrypoint else ''
    full_cmd: str = f'docker run --rm --shm-size={shm_size} -v {ROOT}:/data -w /data {entrypoint_arg} {image} {command}'
    print(f'Running: {full_cmd}')
    subprocess.run(full_cmd, shell=True, check=True)  # noqa: S602
