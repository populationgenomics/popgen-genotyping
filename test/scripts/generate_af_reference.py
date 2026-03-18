"""
Generate a population allele frequency reference VCF from a sites file.
"""

import argparse
import gzip
import random
from pathlib import Path


def generate_aligned_af_vcf(sites_file: str | Path, output_vcf: str | Path, seed: int = 42) -> None:
    """
    Generate an aligned AF reference VCF from a sites file (CHROM POS REF ALT).

    Args:
        sites_file (str | Path): Input sites file, optionally gzipped.
        output_vcf (str | Path): Path to the output AF reference VCF.
        seed (int): Random seed for reproducible AF generation. Defaults to 42.
    """
    rng: random.Random = random.Random(seed)  # noqa: S311

    print(f'Generating aligned AF reference from {sites_file} to {output_vcf}...')

    # Open sites_file with gzip if it ends in .gz, otherwise open normally
    is_gz: bool = str(sites_file).endswith('.gz')
    open_func = gzip.open if is_gz else open
    mode: str = 'rt' if is_gz else 'r'

    with open_func(sites_file, mode) as f_sites, open(output_vcf, 'w') as f_out:
        f_out.write('##fileformat=VCFv4.2\n')
        f_out.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        f_out.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        f_out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

        # Use a constant for the magic column count
        min_columns: int = 4
        for i, line in enumerate(f_sites):
            parts: list[str] = line.strip().split('\t')
            if len(parts) < min_columns:
                continue
            chrom, pos, ref, alt = parts[:4]
            af: float = rng.betavariate(0.5, 0.5)
            f_out.write(f'{chrom}\t{pos}\tsnp{i}\t{ref}\t{alt}\t.\tPASS\tAF={af:.4f}\n')


def main() -> None:
    """
    Main entry point for AF reference generation.
    """
    parser = argparse.ArgumentParser(description='Generate population AF reference VCF from sites file.')
    parser.add_argument('sites_file', help='Input sites.txt file (CHROM POS REF ALT)')
    parser.add_argument('output_vcf', help='Output AF reference VCF file')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for AF generation')

    args: argparse.Namespace = parser.parse_args()
    generate_aligned_af_vcf(args.sites_file, args.output_vcf, args.seed)


if __name__ == '__main__':
    main()
