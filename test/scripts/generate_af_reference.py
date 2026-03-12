import argparse
import random
import os
import gzip

def generate_aligned_af_vcf(sites_file, output_vcf, seed=42):
    rng = random.Random(seed)
    
    print(f"Generating aligned AF reference from {sites_file} to {output_vcf}...")
    
    # Open sites_file with gzip if it ends in .gz, otherwise open normally
    open_func = gzip.open if sites_file.endswith('.gz') else open
    mode = 'rt' if sites_file.endswith('.gz') else 'r'
    
    with open_func(sites_file, mode) as f_sites, open(output_vcf, "w") as f_out:
        f_out.write("##fileformat=VCFv4.2\n")
        f_out.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        f_out.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
        f_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        for i, line in enumerate(f_sites):
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            chrom, pos, ref, alt = parts[:4]
            af = rng.betavariate(0.5, 0.5)
            f_out.write(f"{chrom}\t{pos}\tsnp{i}\t{ref}\t{alt}\t.\tPASS\tAF={af:.4f}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate population AF reference VCF from sites file.")
    parser.add_argument("sites_file", help="Input sites.txt file (CHROM POS REF ALT)")
    parser.add_argument("output_vcf", help="Output AF reference VCF file")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for AF generation")
    
    args = parser.parse_args()
    generate_aligned_af_vcf(args.sites_file, args.output_vcf, args.seed)
