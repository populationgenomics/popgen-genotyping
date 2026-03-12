import argparse
import random
import struct
from pathlib import Path


# Tag Definitions
TAGS = {
    'num_snps': 1, 'sample_name': 101, 'norm_xforms': 400,
    'raw_x': 1000, 'raw_y': 1001, 'genotypes': 1002,
    'base_calls': 1003, 'scores': 1004, 'baf': 1012, 'lrr': 1013
}


def leb128_encode(n: int) -> bytes:
    out: bytearray = bytearray()
    while True:
        byte: int = n & 0x7F
        n >>= 7
        if n == 0:
            out.append(byte)
            return bytes(out)
        out.append(byte | 0x80)

def generate_gtc(output_path: str | Path, num_snps: int, sample_name: str, 
                 contamination: float = 0.05, seed: int = 42) -> list[float]:
    """
    Generates a GTC where BAF is a function of Allele Frequency and Contamination.
    Returns the list of population frequencies used.
    """
    print(f'Generating GTC: {sample_name} | Contamination: {contamination*100}% | Seed: {seed}')
    rng = random.Random(seed)



    # 1. Biological Model Generation
    pop_afs = []
    gts = []
    baf_list = []

    # Standard array noise (standard deviation)
    sigma = 0.035

    for _ in range(num_snps):
        # Generate a population Allele Frequency (AF) for this site
        af = rng.betavariate(0.5, 0.5) 
        pop_afs.append(af)

        # Hardy-Weinberg Equilibrium for Genotypes
        # p^2 (AA), 2pq (AB), q^2 (BB)
        p_aa, p_ab, p_bb = (1-af)**2, 2*af*(1-af), af**2
        gt = rng.choices([1, 2, 3], weights=[p_aa, p_ab, p_bb], k=1)[0]
        gts.append(gt)

        # Contamination Model: BAF = (1-e)*Expected + e*Population_AF
        # For AA (Expected BAF 0): BAF ~ epsilon * af
        # For BB (Expected BAF 1): BAF ~ 1 - [epsilon * (1-af)]
        drift = rng.uniform(-0.02, 0.02)
        if gt == 1: # AA
            mu = (contamination * af) + drift
            baf_val = rng.gauss(mu, sigma)
        elif gt == 3: # BB
            mu = 1.0 - (contamination * (1.0 - af)) + drift
            baf_val = rng.gauss(mu, sigma)
        else: # AB
            baf_val = rng.gauss(0.5, 0.025)
        if rng.random() < 0.02: # 2% outliers
            baf_val = rng.uniform(0,1)
        baf_list.append(max(0.0, min(1.0, baf_val)))

    # 2. Binary Encoding
    def make_block(data: bytes) -> bytes:
        return struct.pack('<i', num_snps) + data

    genotypes_data = bytes(gts)
    scores_data = struct.pack(f'<{num_snps}f', *([0.99 for _ in range(num_snps)]))
    baf_data = struct.pack(f'<{num_snps}f', *baf_list)
    lrr_data = struct.pack(f'<{num_snps}f', *([rng.gauss(0, 0.05) for _ in range(num_snps)]))

    # Raw intensities (simplified)
    raw_x_data = struct.pack(f'<{num_snps}H', *([1000 for _ in range(num_snps)]))
    raw_y_data = struct.pack(f'<{num_snps}H', *([1000 for _ in range(num_snps)]))

    # Base Calls (Dummy placeholders)
    base_calls_data = b''.join([b'AA' if g==1 else b'AB' if g==2 else b'BB' for g in gts])

    # Normalization (Required structure)
    xforms_data = struct.pack('<i', 1) + struct.pack('<i', 1) + struct.pack('<12f', *([0.0]*12))

    name_bytes = sample_name.encode('ascii')
    sample_name_data = leb128_encode(len(name_bytes)) + name_bytes

    tags = [
        (TAGS['num_snps'], struct.pack('<i', num_snps)),
        (TAGS['sample_name'], sample_name_data),
        (TAGS['norm_xforms'], xforms_data),
        (TAGS['raw_x'], make_block(raw_x_data)),
        (TAGS['raw_y'], make_block(raw_y_data)),
        (TAGS['genotypes'], make_block(genotypes_data)),
        (TAGS['base_calls'], make_block(base_calls_data)),
        (TAGS['scores'], make_block(scores_data)),
        (TAGS['baf'], make_block(baf_data)),
        (TAGS['lrr'], make_block(lrr_data)),
    ]

    # Write File
    current_offset = 8 + (len(tags) * 6)
    toc, blocks = [], []
    for ident, data in tags:
        toc.append(struct.pack('<Hi', ident, current_offset))
        blocks.append(data)
        current_offset += len(data)

    with open(output_path, 'wb') as f:
        f.write(b'gtc')
        f.write(struct.pack('<bi', 5, len(tags)))
        for entry in toc: f.write(entry)
        for block in blocks: f.write(block)

    return pop_afs

def write_af_vcf(output_path: str, afs: list[float]):
    """Creates a sites-only VCF to provide BAFregress with the 'true' AF."""
    with open(output_path, 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for i, af in enumerate(afs):
            f.write(f"chr1\t{i+1}\tsnp{i}\tA\tG\t.\tPASS\tAF={af:.4f}\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=str, help='Output GTC path')
    parser.add_argument('--num', type=int, default=5000, help='Number of SNPs')
    parser.add_argument('--contam', type=float, default=0.05, help='Contamination (0.0 to 1.0)')
    parser.add_argument('--af_out', type=str, default='pop_ref.vcf', help='Path for sidecar AF VCF')

    args = parser.parse_args()

    # Generate the GTC
    pop_afs = generate_gtc(args.output, args.num, Path(args.output).stem, contamination=args.contam)

    # Generate the AF reference file
    write_af_vcf(args.af_out, pop_afs)
    print(f"Generated sidecar AF reference: {args.af_out}")

if __name__ == '__main__':
    main()
