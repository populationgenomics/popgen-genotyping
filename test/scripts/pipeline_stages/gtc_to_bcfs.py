"""
Stage: GtcToBcfs reproduction.
"""

from pathlib import Path

from scripts.testing_utils import (
    BCFTOOLS_IMAGE,
    BPM_HOST,
    DATA_DIR,
    EGT_HOST,
    FASTA_HOST,
    run_docker,
    to_container,
)


def run_gtc_to_bcfs(samples: list[str], gtc_paths: list[str]) -> tuple[Path, Path]:
    """
    Run the GtcToBcfs stage in Docker.

    Args:
        samples (list[str]): List of sample names.
        gtc_paths (list[str]): List of container paths to GTC files.

    Returns:
        tuple[Path, Path]: (heavy_bcf_host_path, light_bcf_host_path)
    """
    print('\n>>> Stage: GtcToBcfs (Cohort Level) <<<')
    heavy_bcf: Path = DATA_DIR / 'cohort.heavy.bcf'
    light_bcf: Path = DATA_DIR / 'cohort.light.bcf'
    mapping_file: Path = DATA_DIR / 'reheader_map.txt'

    with open(mapping_file, 'w') as f:
        for s in samples:
            f.write(f'{s} {s}\n')

    gtc_args: str = ' '.join(gtc_paths)
    bpm_int: str = to_container(BPM_HOST)
    egt_int: str = to_container(EGT_HOST)
    fasta_int: str = to_container(FASTA_HOST)
    map_int: str = to_container(mapping_file)
    heavy_int: str = to_container(heavy_bcf)
    light_int: str = to_container(light_bcf)

    gtc_cmd: str = (
        f"bash -c 'mkdir -p ./bcftools-tmp && "
        f'bcftools +gtc2vcf --do-not-check-bpm -b {bpm_int} -e {egt_int} -f {fasta_int} {gtc_args} | '
        f'bcftools norm -m -both --no-version -c x -f {fasta_int} | '
        f'bcftools sort -T ./bcftools-tmp | '
        f'bcftools reheader -s {map_int} | '
        f'bcftools view -O b -o {heavy_int} --write-index && '
        f"bcftools annotate --no-version -x ^FORMAT/GT,FORMAT/GQ {heavy_int} -O b -o {light_int} --write-index'"
    )
    run_docker(BCFTOOLS_IMAGE, gtc_cmd)

    return heavy_bcf, light_bcf
