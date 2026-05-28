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
from scripts.pipeline_stages.king_ibdseg import run_king_ibdseg
from scripts.prepare_data import generate_gtcs, prepare_af_reference
from scripts.testing_utils import DATA_DIR


def main() -> None:
    """
    Execute the full genotyping pipeline reproduction.
    """
    parser = argparse.ArgumentParser(description='Reproduce full genotyping pipeline locally.')
    parser.add_argument('--samples', type=int, default=50, help='Number of samples to simulate')
    parser.add_argument(
        '--snps',
        type=int,
        default=500000,
        help='Number of SNPs to simulate. >=119,078 is needed to include any chrX loci '
        '(default 500,000 covers ~5,600 chrX loci against the GDA-8v1-0_D2 BPM).',
    )
    args: argparse.Namespace = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(DATA_DIR, exist_ok=True)

    # Define samples
    samples: list[str] = [f'CPGSYN{i:03d}' for i in range(1, args.samples + 1)]

    # 1. Data Preparation
    gtc_paths_int, sex_mapping = generate_gtcs(samples, args.snps)
    af_ref_int: str = prepare_af_reference(args.snps)

    # 2. Stage: GtcToBcfs
    heavy_bcf_host, light_bcf_host = run_gtc_to_bcfs(samples, gtc_paths_int)

    # 3. Stage: BafRegress
    baf_out_host = run_bafregress(heavy_bcf_host, af_ref_int)

    # 4. Stage: CohortBcfToPlink
    plink1_prefix_host = run_cohort_bcf_to_plink(samples, light_bcf_host, sex_mapping)

    # 5. Stage: ExportCohortDatasets
    _final_bcf_host = run_export_cohort_datasets(plink1_prefix_host)

    # 6. Stage: KingIbdseg (runs against the merged PLINK 1.9 dataset; uses the
    #    BAFRegress output to exclude contaminated samples before IBD inference).
    king_outputs = run_king_ibdseg(plink1_prefix_host, [baf_out_host])
    print('\nKING outputs:')
    for key, path in king_outputs.items():
        exists = 'OK' if path.exists() else 'MISSING'
        print(f'  {key:<12} {path}  [{exists}]')

    print('\nModular pipeline reproduction complete!')


if __name__ == '__main__':
    main()
