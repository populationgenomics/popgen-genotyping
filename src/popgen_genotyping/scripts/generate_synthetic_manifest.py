"""
Script to generate a deterministic synthetic genotyping manifest CSV for testing.
"""

import csv
import argparse
from pathlib import Path


def generate_manifest(output_path: str | Path, num_samples: int = 10, prefix: str = 'CPGSYN'):
    """
    Generate a synthetic manifest CSV with all required columns.

    Args:
        output_path (str | Path): Path to write the CSV file.
        num_samples (int): Number of sample entries to generate.
        prefix (str): Prefix for the sequencing group IDs (e.g., 'CPGSYN').
    """
    headers = [
        'file_name',
        'sample_sheet_id',
        'sample_id',
        'md5sum',
        'bbv_barcode',
        'sentrix_barcode_a',
        'sentrix_position_a',
        'sample_plate',
        'sample_well',
        'sample_plate_position',
        'cpg_sample_id_internal',
        'cpg_gcp_filepath',
        'cpg_sequencing_group_id',
        'cpg_cohort_id'
    ]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()

        for i in range(1, num_samples + 1):
            sg_id = f'{prefix}{i:03d}'
            sample_id = f'AGDB{i:05d}'
            sentrix_barcode = '207794440004'
            sentrix_pos = f'R{i:02d}C01'

            writer.writerow({
                'file_name': f'{sentrix_barcode}_{sentrix_pos}.gtc',
                'sample_sheet_id': str(835729021 + i),
                'sample_id': sample_id,
                'md5sum': 'deterministic_md5_placeholder',
                'bbv_barcode': str(835729021 + i),
                'sentrix_barcode_a': sentrix_barcode,
                'sentrix_position_a': sentrix_pos,
                'sample_plate': 'AGDB_GDA_001',
                'sample_well': f'{sentrix_barcode}_{sentrix_pos}',
                'sample_plate_position': 'A01',  # simplified
                'cpg_sample_id_internal': f'XPGLCL{50000 + i}',
                'cpg_gcp_filepath': f'gs://cpg-test-main/gtc/{sg_id}.gtc',
                'cpg_sequencing_group_id': sg_id,
                'cpg_cohort_id': 'COH1003'
            })


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate synthetic genotyping manifest CSV.')
    parser.add_argument('output', type=str, help='Output CSV path')
    parser.add_argument('--num', type=int, default=10, help='Number of samples')
    parser.add_argument('--prefix', type=str, default='CPGSYN', help='SG ID prefix')

    args = parser.parse_args()
    generate_manifest(args.output, args.num, args.prefix)
