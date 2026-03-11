#!/usr/bin/env python3  # noqa: EXE001

"""
Reproduction script for the full genotyping pipeline:
1. Generate synthetic GTC files.
2. Convert GTC to BCF using bcftools:1.23-1 (matching gtc_to_bcfs_job.py).
3. Convert BCF to PLINK 1.9 and merge using official PLINK image.
4. Export to PLINK2 and BCF formats.
"""

import json
import subprocess
from pathlib import Path
from typing import Any

# Paths
ROOT: Path = Path(__file__).parent.parent.parent
LOCAL_DIR: Path = ROOT / 'test' / 'local'
DATA_DIR: Path = LOCAL_DIR / 'data'
REF_DIR: Path = LOCAL_DIR / 'reference'
OUTPUT_DIR: Path = LOCAL_DIR / 'output' / 'reproduce'
GTC_DIR: Path = DATA_DIR / 'gtc'

# Docker Images
BCFTOOLS_IMAGE: str = 'australia-southeast1-docker.pkg.dev/cpg-common/images/bcftools:1.23-1'
PLINK_IMAGE: str = 'australia-southeast1-docker.pkg.dev/cpg-common/images/plink:1.9-20250819-PLINK-2.0-20260228-1'

# Reference Files
BPM_MANIFEST: str = '/data/test/local/reference/GDA-8v1-0_D2.bpm'
EGT_CLUSTER: str = '/data/test/local/reference/GDA-8v1-0_D1_ClusterFile.egt'
FASTA_REF: str = '/data/test/local/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'


def run_cmd(cmd: list[str] | str, check: bool = True) -> None:
    """
    Execute a shell command using subprocess.

    Args:
        cmd (list[str] | str): The command to run.
        check (bool): If True, raise CalledProcessError on non-zero exit. Defaults to True.
    """
    print(f"Running: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    subprocess.run(cmd, shell=isinstance(cmd, str), check=check)  # noqa: S603


def run_docker(image: str, command: str, entrypoint: str | None = None) -> None:
    """
    Run a command inside a Docker container.

    Args:
        image (str): Docker image name.
        command (str): Command to execute inside the container.
        entrypoint (str, optional): Custom entrypoint for the container.
    """
    docker_cmd: list[str] = [
        'docker', 'run', '--rm',
        '-v', f'{ROOT}:/data',
        '-w', '/data'
    ]
    if entrypoint:
        docker_cmd.extend(['--entrypoint', entrypoint])

    docker_cmd.append(image)

    if entrypoint == 'bash':
        docker_cmd.extend(['-c', f'set -ex; {command}'])
    else:
        docker_cmd.append(command)

    run_cmd(docker_cmd)


def main() -> None:
    """
    Execute the full genotyping pipeline reproduction.
    """
    # 1. Prepare Directories
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / 'bcfs').mkdir(exist_ok=True)
    (OUTPUT_DIR / 'plink').mkdir(exist_ok=True)

    # 2. Generate Synthetic GTCs
    samples: list[str] = ['CPGSYN001', 'CPGSYN002']
    num_snps: int = 100000  # Reduced for speed in reproduction

    for sg_id in samples:
        gtc_path: Path = GTC_DIR / f'{sg_id}.gtc'
        if not gtc_path.exists():
            print(f'\n>>> Generating GTC for {sg_id} <<<')
            run_cmd([
                'python3', str(ROOT / 'test' / 'scripts' / 'generate_synthetic_gtc.py'),
                str(gtc_path),
                '--num', str(num_snps),
                '--name', 'GDA-8v1-0_D2.bpm'
            ])

    # 3. GTC to BCF (matching gtc_to_bcfs_job.py)
    for sg_id in samples:
        print(f'\n>>> GTC to BCF for {sg_id} <<<')
        gtc_rel: str = f'/data/test/local/data/gtc/{sg_id}.gtc'
        heavy_bcf: str = f'/data/test/local/output/reproduce/bcfs/{sg_id}.heavy.bcf'
        light_bcf: str = f'/data/test/local/output/reproduce/bcfs/{sg_id}.light.bcf'
        metadata_tsv: str = f'/data/test/local/output/reproduce/bcfs/{sg_id}_metadata.tsv'

        cmd: str = f"""
        mkdir -p /tmp/bcftools-tmp && \\
        bcftools +gtc2vcf \\
            --no-version \\
            --bpm {BPM_MANIFEST} \\
            --egt {EGT_CLUSTER} \\
            --fasta-ref {FASTA_REF} \\
            --extra {metadata_tsv} \\
            {gtc_rel} | \\
        bcftools norm -m -both --no-version -c x -f {FASTA_REF} | \\
        bcftools sort -T /tmp/bcftools-tmp | \\
        bcftools reheader -n {sg_id} | \\
        bcftools view -O b -o {heavy_bcf} --write-index

        bcftools annotate --no-version -x ^FORMAT/GT,FORMAT/GQ {heavy_bcf} \\
        -O b -o {light_bcf} --write-index
        """
        run_docker(BCFTOOLS_IMAGE, cmd, entrypoint='bash')

    # 4. Cohort Merge (PLINK 1.9)
    print('\n>>> Merging into Cohort (PLINK 1.9) <<<')
    manifest: dict[str, Any] = {
        'manifest': {
            sg_id: f'/data/test/local/output/reproduce/bcfs/{sg_id}.light.bcf'
            for sg_id in samples
        }
    }
    manifest_path: Path = OUTPUT_DIR / 'manifest.json'
    with open(manifest_path, 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=4)

    sex_tsv: Path = OUTPUT_DIR / 'sex.tsv'
    with open(sex_tsv, 'w', encoding='utf-8') as f:
        for sg_id in samples:
            f.write(f'0\t{sg_id}\t1\n')

    # Run vcf_to_plink.py via official PLINK image
    rel_manifest: str = '/data/test/local/output/reproduce/manifest.json'
    rel_sex: str = '/data/test/local/output/reproduce/sex.tsv'
    rel_out_prefix: str = '/data/test/local/output/reproduce/plink/cohort_merged'

    convert_cmd: str = (
        f'python3 /data/src/popgen_genotyping/scripts/vcf_to_plink.py '
        f'--manifest {rel_manifest} --sex-tsv {rel_sex} --out-prefix {rel_out_prefix} --threads 4'
    )
    run_docker(PLINK_IMAGE, convert_cmd, entrypoint='bash')

    # 5. Export Datasets (PLINK2 & BCF)
    print('\n>>> Exporting Datasets (PLINK2 & BCF) <<<')
    export_dir: Path = OUTPUT_DIR / 'export'
    export_dir.mkdir(exist_ok=True)

    rel_merged_bfile: str = '/data/test/local/output/reproduce/plink/cohort_merged'
    rel_export_prefix: str = '/data/test/local/output/reproduce/export/cohort'

    export_cmd: str = (
        f'plink2 --bfile {rel_merged_bfile} --allow-extra-chr '
        f'--make-pgen --export bcf --out {rel_export_prefix}'
    )
    run_docker(PLINK_IMAGE, export_cmd, entrypoint='bash')

    print('\nFull pipeline reproduction complete!')


if __name__ == '__main__':
    print('Ensure Docker Desktop is running (open -a Docker)')
    main()
