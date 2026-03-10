#!/usr/bin/env python3

"""
Reproduction script for the full genotyping pipeline:
1. Generate synthetic GTC files.
2. Convert GTC to BCF using bcftools:1.23-1 (matching gtc_to_bcfs_job.py).
3. Convert BCF to PLINK 1.9 and merge using plink-multi:local.
"""

import subprocess
import os
import sys
import json
from pathlib import Path

# Paths
ROOT = Path(__file__).parent.parent.parent
LOCAL_DIR = ROOT / 'test' / 'local'
DATA_DIR = LOCAL_DIR / 'data'
REF_DIR = LOCAL_DIR / 'reference'
OUTPUT_DIR = LOCAL_DIR / 'output' / 'reproduce'
GTC_DIR = DATA_DIR / 'gtc'

# Docker Images
BCFTOOLS_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images-dev/bcftools:1.23-1'
PLINK_IMAGE = 'plink-multi:local'

# Reference Files
BPM_MANIFEST = '/data/test/local/reference/GDA-8v1-0_D2.bpm'
EGT_CLUSTER = '/data/test/local/reference/GDA-8v1-0_D1_ClusterFile.egt'
FASTA_REF = '/data/test/local/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'

def run_cmd(cmd, check=True):
    print(f"Running: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    subprocess.run(cmd, shell=isinstance(cmd, str), check=check)

def run_docker(image, command, entrypoint=None):
    docker_cmd = [
        'docker', 'run', '--rm',
        '-v', f"{ROOT}:/data",
        '-w', '/data'
    ]
    if entrypoint:
        docker_cmd.extend(['--entrypoint', entrypoint])
    
    docker_cmd.append(image)
    
    if entrypoint == 'bash':
        docker_cmd.extend(['-c', f"set -ex; {command}"])
    else:
        docker_cmd.append(command)
        
    run_cmd(docker_cmd)

def main():
    # 1. Prepare Directories
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / 'bcfs').mkdir(exist_ok=True)
    (OUTPUT_DIR / 'plink').mkdir(exist_ok=True)
    
    # 2. Generate Synthetic GTCs
    samples = ['CPGSYN001', 'CPGSYN002']
    num_snps = 100000 # Reduced for speed in reproduction
    
    for sg_id in samples:
        gtc_path = GTC_DIR / f"{sg_id}.gtc"
        if not gtc_path.exists():
            print(f"\n>>> Generating GTC for {sg_id} <<<")
            run_cmd([
                'python3', str(ROOT / 'test' / 'scripts' / 'generate_synthetic_gtc.py'),
                str(gtc_path),
                '--num', str(num_snps),
                '--name', 'GDA-8v1-0_D2.bpm'
            ])
            
    # 3. GTC to BCF (matching gtc_to_bcfs_job.py)
    for sg_id in samples:
        print(f"\n>>> GTC to BCF for {sg_id} <<<")
        gtc_rel = f"/data/test/local/data/gtc/{sg_id}.gtc"
        heavy_bcf = f"/data/test/local/output/reproduce/bcfs/{sg_id}.heavy.bcf"
        light_bcf = f"/data/test/local/output/reproduce/bcfs/{sg_id}.light.bcf"
        metadata_tsv = f"/data/test/local/output/reproduce/bcfs/{sg_id}_metadata.tsv"
        
        # Exact command from gtc_to_bcfs_job.py
        # Note: Added --allow-extra-chr to bcftools reheader if needed, but the original didn't have it.
        # Actually, let's keep it exact.
        cmd = f"""
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
    print("\n>>> Merging into Cohort (PLINK 1.9) <<<")
    manifest = {
        "manifest": {
            sg_id: f"/data/test/local/output/reproduce/bcfs/{sg_id}.light.bcf"
            for sg_id in samples
        }
    }
    manifest_path = OUTPUT_DIR / 'manifest.json'
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=4)
        
    sex_tsv = OUTPUT_DIR / 'sex.tsv'
    with open(sex_tsv, 'w') as f:
        for sg_id in samples:
            f.write(f"0\t{sg_id}\t1\n")
            
    # Run vcf_to_plink.py via plink-multi:local
    rel_manifest = "/data/test/local/output/reproduce/manifest.json"
    rel_sex = "/data/test/local/output/reproduce/sex.tsv"
    rel_out_prefix = "/data/test/local/output/reproduce/plink/cohort_merged"
    
    cmd = f"python3 /data/src/popgen_genotyping/scripts/vcf_to_plink.py --manifest {rel_manifest} --sex-tsv {rel_sex} --out-prefix {rel_out_prefix} --threads 4"
    run_docker(PLINK_IMAGE, cmd, entrypoint='bash')

    # 5. Export Datasets (PLINK2 & BCF)
    print("\n>>> Exporting Datasets (PLINK2 & BCF) <<<")
    export_dir = OUTPUT_DIR / 'export'
    export_dir.mkdir(exist_ok=True)
    
    rel_merged_bfile = "/data/test/local/output/reproduce/plink/cohort_merged"
    rel_export_prefix = "/data/test/local/output/reproduce/export/cohort"
    
    # Combined command matching ExportCohortDatasets job
    export_cmd = f"plink2 --bfile {rel_merged_bfile} --allow-extra-chr --make-pgen --export bcf --out {rel_export_prefix}"
    run_docker(PLINK_IMAGE, export_cmd, entrypoint='bash')

    print("\nFull pipeline reproduction complete!")

if __name__ == '__main__':
    # Ensure Docker Desktop is mentioned
    print("Ensure Docker Desktop is running (open -a Docker)")
    main()
