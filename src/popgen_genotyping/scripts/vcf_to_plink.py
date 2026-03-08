#!/usr/bin/env python3

"""
Parallel VCF/BCF to PLINK2 conversion and merging.
Only uses the Python standard library and expects local file paths.
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


def convert_bcf_to_pgen(sg_id: str, local_bcf: str, output_dir: str) -> str:
    """Convert a single BCF to PLINK2 format."""
    out_prefix = os.path.join(output_dir, sg_id)

    # Convert to PGEN
    run_command(f"plink2 --bcf {local_bcf} --make-pgen --out {out_prefix}")

    return out_prefix


def main():
    parser = argparse.ArgumentParser(description="Convert and merge BCFs to PLINK2")
    parser.add_argument("--manifest", required=True, help="Path to local JSON manifest")
    parser.add_argument("--out-prefix", required=True, help="Local output prefix for merged files")
    parser.add_argument("--threads", type=int, default=8, help="Number of parallel threads")
    args = parser.parse_args()

    # Read the manifest
    with open(args.manifest) as f:
        manifest = json.load(f)

    # Dictionary of {sg_id: local_bcf_path}
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
    run_command(f"plink2 --pmerge-list {mergelist_path} --make-pgen --out {args.out_prefix}")

    logger.info("Done!")


if __name__ == "__main__":
    main()
