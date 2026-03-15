"""
Generate a deterministic synthetic GTC file for testing.
"""

import argparse
import random
import struct
from pathlib import Path

# --- Constants ---
GTC_VERSION: int = 5
GTC_HEADER_ID: bytes = b'gtc'
DEFAULT_SNP_COUNT: int = 1000
DEFAULT_CONTAMINATION: float = 0.05
OUTLIER_PERCENTAGE: float = 0.02
DEFAULT_SEED: int = 42
AF_SEED: int = 1337  # Fixed seed for population AFs to be consistent across a cohort

# Official Illumina Tags
TAG_NUM_SNPS: int = 1
TAG_SAMPLE_NAME: int = 101
TAG_NORM_XFORMS: int = 400
TAG_RAW_X: int = 1000
TAG_RAW_Y: int = 1001
TAG_GENOTYPES: int = 1002
TAG_BASE_CALLS: int = 1003
TAG_GENCALL_SCORES: int = 1004
TAG_B_ALLELE_FREQS: int = 1012
TAG_NUM_SNPS_INT: int = 1013


def leb128_encode(n: int) -> bytes:
    """
    Encode an integer using Little Endian Base 128.

    Args:
        n (int): The integer to encode.

    Returns:
        bytes: The LEB128 encoded bytes.
    """
    out: bytearray = bytearray()
    while True:
        byte: int = n & 0x7F
        n >>= 7
        if n == 0:
            out.append(byte)
            return bytes(out)
        out.append(byte | 0x80)


def generate_gtc(  # noqa: PLR0915
    output_path: str | Path,
    num_snps: int,
    sample_name: str,
    contamination: float = DEFAULT_CONTAMINATION,
    sample_seed: int = DEFAULT_SEED,
) -> list[float]:
    """
    Generate a realistic deterministic synthetic GTC for a given number of loci.

    Args:
        output_path (str | Path): Path to save the GTC file.
        num_snps (int): Number of SNPs to generate.
        sample_name (str): Name of the sample to embed in the GTC.
        contamination (float): Simulated contamination level. Defaults to 0.05.
        sample_seed (int): Random seed for sample-specific genotypes. Defaults to 42.

    Returns:
        list[float]: The generated population allele frequencies for the sites.
    """
    print(f'Generating GTC: {sample_name} | Contamination: {contamination * 100}% | Seed: {sample_seed}')

    # 1. Simulate population AFs (Consistent across cohort)
    af_rng: random.Random = random.Random(AF_SEED)  # noqa: S311
    pop_afs: list[float] = [af_rng.betavariate(0.5, 0.5) for _ in range(num_snps)]

    # 2. Simulate sample genotypes (Unique per sample)
    sample_rng: random.Random = random.Random(sample_seed)  # noqa: S311

    # Genotypes: 1=AA, 2=AB, 3=BB, 0=NC
    gts: list[int] = []
    for af in pop_afs:
        p_aa: float = (1 - af) ** 2
        p_ab: float = 2 * af * (1 - af)
        p_bb: float = af**2
        gts.append(sample_rng.choices([1, 2, 3], weights=[p_aa, p_ab, p_bb], k=1)[0])

    # 3. BAF and LRR simulation
    baf_list: list[float] = []
    sigma: float = 0.03
    drift: float = sample_rng.gauss(0, 0.01)

    for i, gt in enumerate(gts):
        af: float = pop_afs[i]
        if gt == 1:  # AA
            mu = (contamination * af) + drift
            baf_val = sample_rng.gauss(mu, sigma)
        elif gt == 3:  # BB  # noqa: PLR2004
            mu = 1.0 - (contamination * (1.0 - af)) + drift
            baf_val = sample_rng.gauss(mu, sigma)
        else:  # AB
            baf_val = sample_rng.gauss(0.5, 0.025)

        if sample_rng.random() < OUTLIER_PERCENTAGE:
            baf_val = sample_rng.uniform(0, 1)
        baf_list.append(max(0.0, min(1.0, baf_val)))

    # 4. Pack data blocks
    def make_array_block(data: bytes) -> bytes:
        return struct.pack('<i', num_snps) + data

    genotypes_data: bytes = make_array_block(bytes(gts))
    baf_data: bytes = make_array_block(struct.pack(f'<{num_snps}f', *baf_list))
    lrr_data: bytes = make_array_block(
        struct.pack(f'<{num_snps}f', *([sample_rng.gauss(0, 0.05) for _ in range(num_snps)]))
    )
    scores_data: bytes = make_array_block(
        struct.pack(f'<{num_snps}f', *([0.98 + sample_rng.random() * 0.02 for _ in range(num_snps)]))
    )
    raw_x_data: bytes = make_array_block(
        struct.pack(f'<{num_snps}H', *([1000 + int(sample_rng.gauss(0, 100)) for _ in range(num_snps)]))
    )
    raw_y_data: bytes = make_array_block(
        struct.pack(f'<{num_snps}H', *([1000 + int(sample_rng.gauss(0, 100)) for _ in range(num_snps)]))
    )

    # Base Calls (Dummy placeholders)
    base_calls_list: list[bytes] = []
    for g in gts:
        if g == 1:
            base_calls_list.append(b'AA')
        elif g == 2:  # noqa: PLR2004
            base_calls_list.append(b'AB')
        else:
            base_calls_list.append(b'BB')
    base_calls_data: bytes = make_array_block(b''.join(base_calls_list))

    # 5. Normalization and Metadata
    xforms_data: bytes = struct.pack('<i', 1) + struct.pack('<i', 1) + struct.pack('<12f', *([0.0] * 12))
    name_bytes: bytes = sample_name.encode('ascii')
    sample_name_data: bytes = leb128_encode(len(name_bytes)) + name_bytes

    tags: list[tuple[int, bytes]] = [
        (TAG_NUM_SNPS, struct.pack('<i', num_snps)),
        (TAG_SAMPLE_NAME, sample_name_data),
        (TAG_NORM_XFORMS, xforms_data),
        (TAG_RAW_X, raw_x_data),
        (TAG_RAW_Y, raw_y_data),
        (TAG_GENOTYPES, genotypes_data),
        (TAG_BASE_CALLS, base_calls_data),
        (TAG_GENCALL_SCORES, scores_data),
        (TAG_B_ALLELE_FREQS, baf_data),
        (TAG_NUM_SNPS_INT, lrr_data),
    ]

    # 6. Write file
    current_offset: int = 8 + (len(tags) * 6)
    toc: list[bytes] = []
    blocks: list[bytes] = []
    for identifier, data in tags:
        toc.append(struct.pack('<Hi', identifier, current_offset))
        blocks.append(data)
        current_offset += len(data)

    with open(output_path, 'wb') as f:
        f.write(GTC_HEADER_ID)
        f.write(struct.pack('<bi', GTC_VERSION, len(tags)))
        for entry in toc:
            f.write(entry)
        for block in blocks:
            f.write(block)

    return pop_afs


def write_af_vcf(output_path: str | Path, afs: list[float]) -> None:
    """
    Write a sidecar AF reference VCF for the synthetic GTC.

    Args:
        output_path (str | Path): Path to save the VCF.
        afs (list[float]): Population allele frequencies.
    """
    with open(output_path, 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Population AF">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for i, af in enumerate(afs):
            f.write(f'chr1\t{i + 1}\t.\tA\tG\t.\t.\tAF={af:.4f}\n')


def main() -> None:
    """
    Main entry point for GTC generation.
    """
    parser = argparse.ArgumentParser(description='Generate synthetic GTC')
    parser.add_argument('output', type=str, help='Output GTC path')
    parser.add_argument('--num', type=int, default=DEFAULT_SNP_COUNT, help='Number of SNPs')
    parser.add_argument('--name', type=str, default='GDA-8v1-0_D2.bpm', help='Sample/Manifest name')
    parser.add_argument('--seed', type=int, default=DEFAULT_SEED, help='Random seed')
    parser.add_argument('--contam', type=float, default=DEFAULT_CONTAMINATION, help='Contamination level')
    parser.add_argument('--af-out', type=str, default='pop_ref.vcf', help='Output AF reference path')

    args: argparse.Namespace = parser.parse_args()
    pop_afs: list[float] = generate_gtc(
        args.output,
        args.num,
        Path(args.output).stem,
        contamination=args.contam,
        sample_seed=args.seed,
    )

    write_af_vcf(args.af_out, pop_afs)
    print(f'Generated sidecar AF reference: {args.af_out}')


if __name__ == '__main__':
    main()
