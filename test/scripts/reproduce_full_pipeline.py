#!/usr/bin/env python3  # noqa: EXE001

"""
Modular reproduction script for the refactored cohort-level genotyping pipeline.
"""

import argparse
import os
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


def generate_gtcs(samples: list[str], num_snps: int) -> list[str]:
    """
    Generate synthetic GTC files for a list of samples.

    Args:
        samples (list[str]): List of sample names.
        num_snps (int): Number of SNPs per GTC.

    Returns:
        list[str]: Internal container paths to the generated GTCs.
    """
    print(f'>>> Generating synthetic GTC files ({len(samples)} samples)...')
    gtc_paths: list[str] = []
    for sample in samples:
        gtc_path: Path = DATA_DIR / f'{sample}.gtc'
        if not gtc_path.exists():
            cmd: list[str] = [
                'python3',
                'test/scripts/generate_synthetic_gtc.py',
                str(gtc_path),
                '--num',
                str(num_snps),
                '--contam',
                '0.05',
            ]
            subprocess.run(cmd, check=True, capture_output=True)  # noqa: S603
        gtc_paths.append(f'/data/{gtc_path.relative_to(ROOT)}')
    return gtc_paths


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

        index_cmd: str = (
            f"bash -c 'bcftools view /data/{af_vcf.relative_to(ROOT)} -O z -o /data/{af_vcf_gz.relative_to(ROOT)} && "
            f"bcftools index /data/{af_vcf_gz.relative_to(ROOT)}'"
        )
        run_docker(BCFTOOLS_IMAGE, index_cmd)
        if af_vcf.exists():
            os.remove(af_vcf)
    return f'/data/{af_vcf_gz.relative_to(ROOT)}'


def main() -> None:
    """
    Execute the full genotyping pipeline reproduction.
    """
    parser = argparse.ArgumentParser(description='Reproduce full genotyping pipeline locally.')
    parser.add_argument('--samples', type=int, default=10, help='Number of samples to simulate')
    parser.add_argument('--snps', type=int, default=100000, help='Number of SNPs to simulate')
    args: argparse.Namespace = parser.parse_args()

    os.makedirs(DATA_DIR, exist_ok=True)
    samples: list[str] = [f'CPGSYN{i:03d}' for i in range(1, args.samples + 1)]

    # 1 & 2. Data Preparation
    gtc_paths: list[str] = generate_gtcs(samples, args.snps)
    af_ref_int: str = prepare_af_reference(args.snps)

    # 3. GtcToBcfs
    print('\n>>> Stage: GtcToBcfs (Cohort Level) <<<')
    heavy_bcf: Path = DATA_DIR / 'cohort.heavy.bcf'
    light_bcf: Path = DATA_DIR / 'cohort.light.bcf'
    mapping_file: Path = DATA_DIR / 'reheader_map.txt'
    with open(mapping_file, 'w') as f:
        for s in samples:
            f.write(f'{s} {s}\n')

    gtc_args: str = ' '.join(gtc_paths)
    bpm_int: str = f'/data/{BPM_HOST.relative_to(ROOT)}'
    egt_int: str = f'/data/{EGT_HOST.relative_to(ROOT)}'
    fasta_int: str = f'/data/{FASTA_HOST.relative_to(ROOT)}'
    map_int: str = f'/data/{mapping_file.relative_to(ROOT)}'
    heavy_int: str = f'/data/{heavy_bcf.relative_to(ROOT)}'
    light_int: str = f'/data/{light_bcf.relative_to(ROOT)}'

    gtc_cmd: str = (
        f"bash -c 'mkdir -p /dev/shm/bcftools-tmp && "
        f'bcftools +gtc2vcf --do-not-check-bpm -b {bpm_int} -e {egt_int} -f {fasta_int} {gtc_args} | '
        f'bcftools norm -m -both --no-version -c x -f {fasta_int} | '
        f'bcftools sort -T /dev/shm/bcftools-tmp | '
        f'bcftools reheader -s {map_int} | '
        f'bcftools view -O b -o {heavy_int} --write-index && '
        f"bcftools annotate --no-version -x ^FORMAT/GT,FORMAT/GQ {heavy_int} -O b -o {light_int} --write-index'"
    )
    run_docker(BCFTOOLS_IMAGE, gtc_cmd)

    # 4. BAFRegress
    print('\n>>> Stage: BafRegress (Cohort Level) <<<')
    baf_out_host: Path = DATA_DIR / 'cohort.BAFRegress.txt'
    baf_out_int: str = f'/data/{baf_out_host.relative_to(ROOT)}'
    baf_cmd: str = f"bash -c 'bcftools +BAFregress -a {af_ref_int} --tag AF {heavy_int} > {baf_out_int}'"
    run_docker(BCFTOOLS_IMAGE, baf_cmd)

    # 5. CohortBcfToPlink
    print('\n>>> Stage: CohortBcfToPlink (PLINK 1.9) <<<')
    plink1_prefix: Path = DATA_DIR / 'cohort_plink1'
    plink1_int: str = f'/data/{plink1_prefix.relative_to(ROOT)}'
    sex_mapping_file: Path = DATA_DIR / 'sex_mapping.txt'
    with open(sex_mapping_file, 'w') as f:
        for s in samples:
            f.write(f'0 {s} 1\n')
    sex_int: str = f'/data/{sex_mapping_file.relative_to(ROOT)}'

    plink1_cmd: str = (
        f'bash -c \'plink2 --bcf {light_int} --max-alleles 2 --split-par hg38 --set-all-var-ids "@:#:\\$r:\\$a" '
        f"--update-sex {sex_int} --allow-extra-chr --make-bed --out {plink1_int}'"
    )
    run_docker(PLINK_IMAGE, plink1_cmd)

    # 6. ExportCohortDatasets
    print('\n>>> Stage: ExportCohortDatasets (PLINK2/BCF) <<<')
    plink2_prefix: Path = DATA_DIR / 'cohort_plink2'
    plink2_int: str = f'/data/{plink2_prefix.relative_to(ROOT)}'
    final_bcf: Path = DATA_DIR / 'cohort.final.bcf'
    final_int: str = f'/data/{final_bcf.relative_to(ROOT)}'

    export_cmd: str = (
        f"bash -c 'plink2 --bfile {plink1_int} --allow-extra-chr --make-pgen --out {plink2_int} && "
        f'plink2 --pfile {plink2_int} --allow-extra-chr --export bcf --out {plink2_int} && '
        f"mv {plink2_int}.bcf {final_int}'"
    )
    run_docker(PLINK_IMAGE, export_cmd)
    run_docker(BCFTOOLS_IMAGE, f'bcftools index {final_int}')

    print('\nRefactored pipeline reproduction complete!')


if __name__ == '__main__':
    main()
