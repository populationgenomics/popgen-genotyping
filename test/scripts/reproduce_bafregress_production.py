#!/usr/bin/env python3  # noqa: EXE001

"""
Test BAFRegress for a production-representative cohort size without AF reference.
"""

import argparse
import os

from scripts.pipeline_stages.bafregress import run_bafregress
from scripts.pipeline_stages.gtc_to_bcfs import run_gtc_to_bcfs
from scripts.prepare_data import generate_gtcs
from scripts.testing_utils import DATA_DIR


def main() -> None:
    """
    Execute the BAFRegress production test.
    """
    parser = argparse.ArgumentParser(description='Test BAFRegress with production cohort size.')
    parser.add_argument('--samples', type=int, default=94, help='Number of samples to simulate')
    parser.add_argument('--snps', type=int, default=100000, help='Number of SNPs to simulate')
    args: argparse.Namespace = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(DATA_DIR, exist_ok=True)

    # Define samples
    samples: list[str] = [f'CPGSYN{i:03d}' for i in range(1, args.samples + 1)]

    # 1. Data Preparation
    # We generate synthetic GTCs for all samples.
    gtc_paths_int: list[str] = generate_gtcs(samples, args.snps)

    # 2. Stage: GtcToBcfs
    # We need the heavy BCF for BAFRegress.
    heavy_bcf_host, _light_bcf_host = run_gtc_to_bcfs(samples, gtc_paths_int)

    # 3. Stage: BafRegress (Without AF reference)
    print('\n>>> Running BAFRegress WITHOUT external AF reference...')
    baf_out_host = run_bafregress(heavy_bcf_host, af_ref_int=None)

    print(f'\nTest complete. BAFRegress output: {baf_out_host}')
    if baf_out_host.exists():
        print('Output file exists. Content preview:')
        with open(baf_out_host) as f:
            for _ in range(5):
                line = f.readline()
                if not line:
                    break
                print(f'  {line.strip()}')


if __name__ == '__main__':
    main()
