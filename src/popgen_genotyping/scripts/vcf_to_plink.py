#!/usr/bin/env python3

"""
Parallel VCF/BCF to PLINK2 conversion and merging.
Only uses the Python standard library.
"""

import argparse
import json
import logging
import os
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def run_command(command: str) -> None:
    """Run a shell command and log its output."""
    logger.info(f"Running command: {command}")
    cmd_args = shlex.split(command)
    result = subprocess.run(cmd_args, check=True, text=True, capture_output=True)
    if result.stdout:
        logger.info(result.stdout)
    if result.stderr:
        logger.warning(result.stderr)


def convert_bcf_to_pgen(sg_id: str, bcf_path: str, output_dir: str) -> str:
    """Convert a single BCF to PLINK2 format."""
    out_prefix = os.path.join(output_dir, sg_id)
    # Localise the BCF first if it's on GCS
    local_bcf = os.path.join(output_dir, f"{sg_id}.bcf")
    if bcf_path.startswith("gs://"):
        run_command(f"gsutil cp {bcf_path} {local_bcf}")
        # Also try to get index
        run_command(f"gsutil cp {bcf_path}.csi {local_bcf}.csi || true")
    else:
        local_bcf = bcf_path

    # Convert to PGEN
    run_command(f"plink2 --bcf {local_bcf} --make-pgen --out {out_prefix}")

    # Remove local BCF to save space
    if bcf_path.startswith("gs://") and os.path.exists(local_bcf):
        os.remove(local_bcf)
        if os.path.exists(f"{local_bcf}.csi"):
            os.remove(f"{local_bcf}.csi")

    return out_prefix


def main():
    parser = argparse.ArgumentParser(description="Convert and merge BCFs to PLINK2")
    parser.add_argument("--manifest", required=True, help="Path to JSON manifest")
    parser.add_argument("--out-prefix", required=True, help="Cloud output prefix (gs://...)")
    parser.add_argument("--threads", type=int, default=8, help="Number of parallel threads")
    args = parser.parse_args()

    # Read the manifest
    local_manifest = "manifest.json"
    if args.manifest.startswith("gs://"):
        run_command(f"gsutil cp {args.manifest} {local_manifest}")
    else:
        local_manifest = args.manifest

    with open(local_manifest) as f:
        manifest = json.load(f)

    # Dictionary of {sg_id: bcf_path}
    samples: dict[str, str] = manifest.get("manifest", {})

    output_dir = "plink_outputs"
    os.makedirs(output_dir, exist_ok=True)

    # 1. Parallel conversion
    logger.info(f"Converting {len(samples)} samples to PLINK2 format...")
    pgen_prefixes = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_sg = {
            executor.submit(convert_bcf_to_pgen, sg_id, bcf_path, output_dir): sg_id
            for sg_id, bcf_path in samples.items()
        }
        for future in as_completed(future_to_sg):
            sg_id = future_to_sg[future]
            try:
                pgen_prefix = future.result()
                pgen_prefixes.append(pgen_prefix)
            except Exception as e:
                logger.error(f"Sample {sg_id} failed: {e}")
                raise

    # 2. Final merge
    mergelist_path = "mergelist.txt"
    with open(mergelist_path, "w") as f:
        for prefix in pgen_prefixes:
            f.write(f"{prefix}\n")

    logger.info("Merging all samples into a single cohort dataset...")
    cohort_out = "cohort_merged"
    run_command(f"plink2 --pmerge-list {mergelist_path} --make-pgen --out {cohort_out}")

    # 3. Upload results
    logger.info(f"Uploading merged files to {args.out_prefix}...")
    for ext in ["pgen", "pvar", "psam"]:
        local_file = f"{cohort_out}.{ext}"
        if os.path.exists(local_file):
            run_command(f"gsutil cp {local_file} {args.out_prefix}.{ext}")

    logger.info("Done!")


if __name__ == "__main__":
    main()
