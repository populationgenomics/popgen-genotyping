#!/usr/bin/env python3

"""
Parallel VCF/BCF to PLINK 1.9 conversion and merging orchestration.
Only uses the Python standard library and expects local file paths.
"""

import argparse
import json
import logging
import os
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def run_command(command: str) -> None:
    """
    Execute a shell command using subprocess.

    Args:
        command (str): The full shell command string to execute.

    Raises:
        subprocess.CalledProcessError: If the command returns a non-zero exit code.
    """
    logger.info(f'Running command: {command}')
    cmd_args = shlex.split(command)
    # Don't capture output so it streams to the console for debugging
    subprocess.run(cmd_args, check=True)  # noqa: S603


def convert_bcf_to_plink19(sg_id: str, local_bcf: str, output_dir: str, sex_tsv: str | None = None) -> str:
    """
    Convert a single BCF to a sorted PLINK 1.9 binary fileset using PLINK2.

    Args:
        sg_id (str): Sequencing Group ID.
        local_bcf (str): Local path to the input BCF file.
        output_dir (str): Local directory to write converted files.
        sex_tsv (str, optional): Local path to 3-column sex metadata TSV.

    Returns:
        str: The local prefix of the generated PLINK binary files.
    """
    out_prefix: str = os.path.join(output_dir, sg_id)

    # Use PLINK2 for initial conversion to leverage --max-alleles 2
    # This filters out multiallelic variants before they hit the PLINK 1.9 merge
    command: str = (
        f'plink2 --bcf {local_bcf} '
        f'--max-alleles 2 '
        f'--split-par hg38 '
        f"--set-all-var-ids '@:#:$r:$a' "
        f'--allow-extra-chr '
        f'--make-bed '
        f'--memory 2000 '
        f'--out {out_prefix}'
    )

    # Add sex update if provided
    if sex_tsv:
        command += f' --update-sex {sex_tsv}'

    run_command(command)

    return out_prefix


def main() -> None:
    """
    Main entry point for the orchestration script.
    """
    parser = argparse.ArgumentParser(description='Convert and merge BCFs to PLINK 1.9')
    parser.add_argument('--manifest', required=True, help='Path to local JSON manifest')
    parser.add_argument('--out-prefix', required=True, help='Local output prefix for merged files')
    parser.add_argument('--sex-tsv', help='Path to local 3-column sex metadata TSV')
    parser.add_argument('--threads', type=int, default=8, help='Number of parallel threads')
    args: argparse.Namespace = parser.parse_args()

    # Read the manifest
    with open(args.manifest, encoding='utf-8') as f:
        manifest: dict[str, Any] = json.load(f)

    # Dictionary of {sg_id: local_bcf_path}
    samples: dict[str, str] = manifest.get('manifest', {})

    if not samples:
        logger.error('No samples found in manifest')
        return

    output_dir: str = 'plink_outputs'
    os.makedirs(output_dir, exist_ok=True)

    # 1. Parallel conversion to sorted bed/bim/fam
    logger.info(f'Converting {len(samples)} samples to PLINK 1.9 format...')
    prefixes: list[str] = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_sg = {
            executor.submit(convert_bcf_to_plink19, sg_id, bcf_path, output_dir, args.sex_tsv): sg_id
            for sg_id, bcf_path in samples.items()
        }
        for future in as_completed(future_to_sg):
            sg_id: str = future_to_sg[future]
            try:
                prefix: str = future.result()
                prefixes.append(prefix)
            except Exception as e:
                logger.error(f'Sample {sg_id} failed: {e}')
                raise

    # 2. Final merge using --merge-list
    # Note: PLINK 1.9 --merge-list expects prefixes of all but the first dataset
    first_prefix: str = prefixes[0]
    rest_prefixes: list[str] = prefixes[1:]

    if not rest_prefixes:
        logger.info('Single sample detected, copying to final output...')
        run_command(f'plink1.9 --bfile {first_prefix} --allow-extra-chr --make-bed --out {args.out_prefix}')
    else:
        mergelist_path: str = 'mergelist_rest.txt'
        with open(mergelist_path, 'w', encoding='utf-8') as f:
            for prefix in rest_prefixes:
                f.write(f'{prefix}\n')

        logger.info(f'Merging {len(rest_prefixes)} samples using --merge-list...')
        run_command(
            f'plink1.9 --bfile {first_prefix} '
            f'--merge-list {mergelist_path} '
            f'--allow-extra-chr '
            f'--make-bed --out {args.out_prefix}'
        )

    logger.info('Done!')


if __name__ == '__main__':
    main()
