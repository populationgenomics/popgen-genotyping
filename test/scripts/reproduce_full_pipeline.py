#!/usr/bin/env python3  # noqa: EXE001

"""
Updated reproduction script for the refactored cohort-level genotyping pipeline.
Optimized for speed using /dev/shm for sort and 10 samples.
Includes all stages: GTC -> multi-sample BCF -> BAFRegress -> cohort PLINK 1.9 -> PLINK2/BCF.
"""

import os
import subprocess
import time
from pathlib import Path

# --- Configuration ---
BCFTOOLS_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/bcftools:1.23-1'
PLINK_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/plink:1.9-20250819-PLINK-2.0-20260228-1'

ROOT = Path(__file__).parent.parent.parent
LOCAL_DIR = ROOT / 'test' / 'local'
DATA_DIR = LOCAL_DIR / 'reproduce_cohort'
SITES_FILE = ROOT / 'test' / 'data' / 'sites.txt.gz'

# Reference files on host
BPM_HOST = LOCAL_DIR / 'reference' / 'GDA-8v1-0_D2.bpm'
EGT_HOST = LOCAL_DIR / 'reference' / 'GDA-8v1-0_D1_ClusterFile.egt'
FASTA_HOST = LOCAL_DIR / 'reference' / 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'

# --- Utilities ---

def run_docker(image, command, entrypoint=None, shm_size='8g'):
    """Run a command inside a Docker container with large SHM size."""
    entrypoint_arg = f'--entrypoint {entrypoint}' if entrypoint else ''
    full_cmd = (
        f'docker run --rm '
        f'--shm-size={shm_size} '
        f'-v {ROOT}:/data '
        f'-w /data '
        f'{entrypoint_arg} '
        f'{image} '
        f'{command}'
    )
    print(f'Running: {full_cmd}')
    subprocess.run(full_cmd, shell=True, check=True)

# --- Execution ---

def main():
    os.makedirs(DATA_DIR, exist_ok=True)
    
    # 1. Generate 10 synthetic samples
    samples = [f'CPGSYN{i:03d}' for i in range(1, 11)]
    gtc_paths = []
    
    num_snps = 100000 
    
    print('>>> Generating synthetic GTC files (10 samples)...')
    for sample in samples:
        gtc_path = DATA_DIR / f'{sample}.gtc'
        if not gtc_path.exists():
            subprocess.run(
                f'python3 test/scripts/generate_synthetic_gtc.py {gtc_path} --num {num_snps} --contam 0.05',
                shell=True, check=True, capture_output=True
            )
        gtc_paths.append(f'/data/{gtc_path.relative_to(ROOT)}')

    # 2. Create pop_ref for BAFRegress (Only once)
    af_vcf_gz = DATA_DIR / 'pop_ref.vcf.gz'
    if not af_vcf_gz.exists():
        print('\n>>> Generating population AF reference (First time only)...')
        af_vcf = DATA_DIR / 'pop_ref.vcf'
        subprocess.run(
            f'python3 test/scripts/generate_af_reference.py {SITES_FILE} {af_vcf}',
            shell=True, check=True, capture_output=True
        )
        # Filter pop_ref to match num_snps
        subprocess.run(f'head -n {num_snps + 4} {af_vcf} > {af_vcf}.tmp && mv {af_vcf}.tmp {af_vcf}', shell=True, check=True)
        
        run_docker(BCFTOOLS_IMAGE, f"bash -c 'bcftools view /data/{af_vcf.relative_to(ROOT)} -O z -o /data/{af_vcf_gz.relative_to(ROOT)} && bcftools index /data/{af_vcf_gz.relative_to(ROOT)}'")
        if af_vcf.exists():
            os.remove(af_vcf)

    # 3. GTC to multi-sample cohort BCF
    print('\n>>> GtcToBcfs (Cohort Level) <<<')
    heavy_bcf = DATA_DIR / 'cohort.heavy.bcf'
    light_bcf = DATA_DIR / 'cohort.light.bcf'
    
    mapping_file = DATA_DIR / 'reheader_map.txt'
    with open(mapping_file, 'w') as f:
        for s in samples:
            f.write(f'{s} {s}\n')

    gtc_args = ' '.join(gtc_paths)
    
    # Internal paths
    bpm_int = f'/data/{BPM_HOST.relative_to(ROOT)}'
    egt_int = f'/data/{EGT_HOST.relative_to(ROOT)}'
    fasta_int = f'/data/{FASTA_HOST.relative_to(ROOT)}'
    map_int = f'/data/{mapping_file.relative_to(ROOT)}'
    heavy_int = f'/data/{heavy_bcf.relative_to(ROOT)}'
    light_int = f'/data/{light_bcf.relative_to(ROOT)}'
    
    gtc_cmd = (
        f"bash -c 'mkdir -p /dev/shm/bcftools-tmp && "
        f"bcftools +gtc2vcf --do-not-check-bpm -b {bpm_int} -e {egt_int} -f {fasta_int} {gtc_args} | "
        f"bcftools norm -m -both --no-version -c x -f {fasta_int} | "
        f"bcftools sort -T /dev/shm/bcftools-tmp | "
        f"bcftools reheader -s {map_int} | "
        f"bcftools view -O b -o {heavy_int} --write-index && "
        f"bcftools annotate --no-version -x ^FORMAT/GT,FORMAT/GQ {heavy_int} -O b -o {light_int} --write-index'"
    )
    run_docker(BCFTOOLS_IMAGE, gtc_cmd)

    # 4. BAFRegress (Cohort Level)
    print('\n>>> BafRegress (Cohort Level) <<<')
    baf_out_host = DATA_DIR / 'cohort.BAFRegress.txt'
    af_ref_int = f'/data/{af_vcf_gz.relative_to(ROOT)}'
    baf_out_int = f'/data/{baf_out_host.relative_to(ROOT)}'
    
    baf_cmd = (
        f"bash -c 'bcftools +BAFregress -a {af_ref_int} --tag AF {heavy_int} > {baf_out_int}'"
    )
    run_docker(BCFTOOLS_IMAGE, baf_cmd)
    
    print('\nBAFRegress Results (First 5 lines):')
    with open(baf_out_host) as f:
        for _ in range(5):
            print(f.readline(), end='')

    # 5. Cohort BCF to PLINK 1.9
    print('\n>>> CohortBcfToPlink (Simplified) <<<')
    plink1_prefix = DATA_DIR / 'cohort_plink1'
    plink1_prefix_int = f'/data/{plink1_prefix.relative_to(ROOT)}'
    
    sex_mapping_file = DATA_DIR / 'sex_mapping.txt'
    with open(sex_mapping_file, 'w') as f:
        for s in samples:
            f.write(f'0 {s} 1\n')
    sex_mapping_int = f'/data/{sex_mapping_file.relative_to(ROOT)}'

    plink1_cmd = (
        f"bash -c 'plink2 --bcf {light_int} "
        f"--max-alleles 2 --split-par hg38 --set-all-var-ids \"@:#:\\$r:\\$a\" "
        f"--update-sex {sex_mapping_int} "
        f"--allow-extra-chr --make-bed --out {plink1_prefix_int}'"
    )
    run_docker(PLINK_IMAGE, plink1_cmd)

    # 6. ExportCohortDatasets (PLINK 1.9 -> PLINK2/BCF)
    print('\n>>> ExportCohortDatasets (PLINK 1.9 -> PLINK2/BCF) <<<')
    plink2_prefix = DATA_DIR / 'cohort_plink2'
    plink2_prefix_int = f'/data/{plink2_prefix.relative_to(ROOT)}'
    final_bcf = DATA_DIR / 'cohort.final.bcf'
    final_bcf_int = f'/data/{final_bcf.relative_to(ROOT)}'

    export_cmd = (
        f"bash -c 'plink2 --bfile {plink1_prefix_int} "
        f"--allow-extra-chr --make-pgen --out {plink2_prefix_int} && "
        f"plink2 --pfile {plink2_prefix_int} --allow-extra-chr --export bcf --out {plink2_prefix_int} && "
        f"mv {plink2_prefix_int}.bcf {final_bcf_int}'"
    )
    run_docker(PLINK_IMAGE, export_cmd)
    
    # Index with bcftools
    run_docker(BCFTOOLS_IMAGE, f"bcftools index {final_bcf_int}")

    print('\nRefactored pipeline reproduction (Full DAG) complete!')

if __name__ == '__main__':
    main()
