#!/usr/bin/env python3  # noqa: EXE001

"""
Modular reproduction script for the refactored cohort-level genotyping pipeline.
"""

import argparse
import os

from scripts.pipeline_stages.bafregress import run_bafregress
from scripts.pipeline_stages.cohort_bcf_to_plink import run_cohort_bcf_to_plink
from scripts.pipeline_stages.export_cohort_datasets import run_export_cohort_datasets
from scripts.pipeline_stages.gtc_to_bcfs import run_gtc_to_bcfs
from scripts.prepare_data import generate_gtcs, prepare_af_reference
from scripts.testing_utils import DATA_DIR


def main() -> None:
    """
    Execute the full genotyping pipeline reproduction.
    """
    parser = argparse.ArgumentParser(description='Reproduce full genotyping pipeline locally.')
    parser.add_argument('--samples', type=int, default=10, help='Number of samples to simulate')
    parser.add_argument('--snps', type=int, default=100000, help='Number of SNPs to simulate')
    args: argparse.Namespace = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(DATA_DIR, exist_ok=True)

    # Define samples
    samples: list[str] = [f'CPGSYN{i:03d}' for i in range(1, args.samples + 1)]

    # 1. Data Preparation
    gtc_paths_int: list[str] = generate_gtcs(samples, args.snps)
    af_ref_int: str = prepare_af_reference(args.snps)

    # 2. Stage: GtcToBcfs
    heavy_bcf_host, light_bcf_host = run_gtc_to_bcfs(samples, gtc_paths_int)

    # 3. Stage: BafRegress
    _baf_out_host = run_bafregress(heavy_bcf_host, af_ref_int)

    # 4. Stage: CohortBcfToPlink
    plink1_prefix_host = run_cohort_bcf_to_plink(samples, light_bcf_host)

    # 5. Stage: ExportCohortDatasets
    _final_bcf_host = run_export_cohort_datasets(plink1_prefix_host)

    print('\nModular pipeline reproduction complete!')


if __name__ == '__main__':
    main()
