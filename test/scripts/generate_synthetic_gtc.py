"""
Generate a deterministic synthetic GTC file for testing.
"""

import argparse
import gzip
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

# BPM chromosome strings that imply hemizygous calls for a given sex.
# - 'X' is the non-PAR chrX block (Illumina marks PAR as 'XY' separately).
# - 'Y' is the non-PAR chrY block.
# - 'MT' is mitochondrial (effectively haploid for both sexes).
HEMIZYGOUS_MALE_CHROMS: frozenset[str] = frozenset({'X', 'Y', 'MT'})
HEMIZYGOUS_FEMALE_CHROMS: frozenset[str] = frozenset({'MT'})
FEMALE_NULL_CHROMS: frozenset[str] = frozenset({'Y'})  # females have no chrY signal


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


def load_bpm_chroms(cache_path: Path, limit: int | None = None) -> list[str]:
    """
    Load the BPM locus-index to chromosome mapping from the gzipped cache.

    Args:
        cache_path (Path): Path to test/data/bpm_chrom.txt.gz.
        limit (int | None): If set, return only the first `limit` chromosomes.

    Returns:
        list[str]: Chromosome strings indexed by (BPM_index - 1).
    """
    chroms: list[str] = []
    with gzip.open(cache_path, 'rt', encoding='utf-8') as f:
        for i, line in enumerate(f):
            if limit is not None and i >= limit:
                break
            chroms.append(line.strip())
    return chroms


def generate_gtc(  # noqa: PLR0915
    output_path: str | Path,
    num_snps: int,
    sample_name: str,
    contamination: float = DEFAULT_CONTAMINATION,
    sample_seed: int = DEFAULT_SEED,
    sex: str = 'M',
    chroms: list[str] | None = None,
) -> list[float]:
    """
    Generate a realistic deterministic synthetic GTC for a given number of loci.

    Args:
        output_path (str | Path): Path to save the GTC file.
        num_snps (int): Number of SNPs to generate.
        sample_name (str): Name of the sample to embed in the GTC.
        contamination (float): Simulated contamination level. Defaults to 0.05.
        sample_seed (int): Random seed for sample-specific genotypes. Defaults to 42.
        sex (str): 'M' or 'F'. Drives ploidy of chrX/chrY/MT calls. Defaults to 'M'.
        chroms (list[str] | None): BPM chromosomes for each locus (length >= num_snps).
            If None, all loci are treated as autosomes (legacy behaviour).

    Returns:
        list[float]: The generated population allele frequencies for the sites.
    """
    print(f'Generating GTC: {sample_name} | Sex: {sex} | Contamination: {contamination * 100}% | Seed: {sample_seed}')

    if chroms is not None and len(chroms) < num_snps:
        raise ValueError(
            f'chroms length {len(chroms)} is shorter than num_snps {num_snps}; '
            f'regenerate the BPM chrom cache or reduce --num',
        )

    sex = sex.upper()
    if sex not in {'M', 'F'}:
        raise ValueError(f"sex must be 'M' or 'F', got {sex!r}")

    # 1. Simulate population AFs (Consistent across cohort)
    af_rng: random.Random = random.Random(AF_SEED)  # noqa: S311
    pop_afs: list[float] = [af_rng.betavariate(0.5, 0.5) for _ in range(num_snps)]

    # 2. Simulate sample genotypes (Unique per sample)
    sample_rng: random.Random = random.Random(sample_seed)  # noqa: S311

    hemi_chroms: frozenset[str] = HEMIZYGOUS_MALE_CHROMS if sex == 'M' else HEMIZYGOUS_FEMALE_CHROMS

    # Genotypes: 1=AA, 2=AB, 3=BB, 0=NC
    gts: list[int] = []
    for i, pop_af in enumerate(pop_afs):
        chrom: str | None = chroms[i] if chroms is not None else None

        if chrom in FEMALE_NULL_CHROMS and sex == 'F':
            # Females have no chrY signal — emit no-call.
            gts.append(0)
            continue

        if chrom in hemi_chroms:
            # Hemizygous: draw one allele only, never AB.
            gts.append(sample_rng.choices([1, 3], weights=[1 - pop_af, pop_af], k=1)[0])
            continue

        # Diploid HWE draw (autosomes, PAR, female chrX).
        p_aa: float = (1 - pop_af) ** 2
        p_ab: float = 2 * pop_af * (1 - pop_af)
        p_bb: float = pop_af**2
        gts.append(sample_rng.choices([1, 2, 3], weights=[p_aa, p_ab, p_bb], k=1)[0])

    # 3. BAF and LRR simulation
    baf_list: list[float] = []
    sigma: float = 0.03
    drift: float = sample_rng.gauss(0, 0.01)

    for i, gt in enumerate(gts):
        af: float = pop_afs[i]
        if gt == 0:
            # No-call → BAF undefined; use a neutral 0.5 sentinel so the
            # downstream BCF still has a value.
            baf_list.append(0.5)
            continue
        if gt == 1:  # AA (or hemizygous A)
            mu = (contamination * af) + drift
            baf_val = sample_rng.gauss(mu, sigma)
        elif gt == 3:  # BB (or hemizygous B)
            mu = 1.0 - (contamination * (1.0 - af)) + drift
            baf_val = sample_rng.gauss(mu, sigma)
        else:  # AB (diploid het only — hemizygous draws never produce this)
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
        elif g == 2:
            base_calls_list.append(b'AB')
        elif g == 3:
            base_calls_list.append(b'BB')
        else:
            base_calls_list.append(b'NN')
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
    repo_root: Path = Path(__file__).resolve().parents[2]
    default_cache: Path = repo_root / 'test' / 'data' / 'bpm_chrom.txt.gz'

    parser = argparse.ArgumentParser(description='Generate synthetic GTC')
    parser.add_argument('output', type=str, help='Output GTC path')
    parser.add_argument('--num', type=int, default=DEFAULT_SNP_COUNT, help='Number of SNPs')
    parser.add_argument('--name', type=str, default='GDA-8v1-0_D2.bpm', help='Sample/Manifest name')
    parser.add_argument('--seed', type=int, default=DEFAULT_SEED, help='Random seed')
    parser.add_argument('--contam', type=float, default=DEFAULT_CONTAMINATION, help='Contamination level')
    parser.add_argument('--af-out', type=str, default='pop_ref.vcf', help='Output AF reference path')
    parser.add_argument('--sex', type=str, default='M', choices=['M', 'F'], help="Sample sex ('M' or 'F')")
    parser.add_argument(
        '--bpm-chrom-cache',
        type=Path,
        default=default_cache,
        help='Path to the gzipped BPM chrom cache (line N = chrom of BPM index N)',
    )

    args: argparse.Namespace = parser.parse_args()

    chroms: list[str] | None = None
    if args.bpm_chrom_cache.exists():
        chroms = load_bpm_chroms(args.bpm_chrom_cache, limit=args.num)
    else:
        print(f'BPM chrom cache not found at {args.bpm_chrom_cache}; treating all loci as autosomes')

    pop_afs: list[float] = generate_gtc(
        args.output,
        args.num,
        Path(args.output).stem,
        contamination=args.contam,
        sample_seed=args.seed,
        sex=args.sex,
        chroms=chroms,
    )

    write_af_vcf(args.af_out, pop_afs)
    print(f'Generated sidecar AF reference: {args.af_out}')


if __name__ == '__main__':
    main()
