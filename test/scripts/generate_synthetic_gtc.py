"""
Generate a deterministic synthetic GTC file for testing.
"""

import argparse
import random
import struct
from pathlib import Path


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


def generate_gtc(output_path: str | Path, num_snps: int, sample_name: str, seed: int = 42) -> None:  # noqa: PLR0915
    """
    Generate a realistic deterministic synthetic GTC for a given number of loci.

    Args:
        output_path (str | Path): Path to save the GTC file.
        num_snps (int): Number of SNPs to generate.
        sample_name (str): Name of the sample to embed in the GTC.
        seed (int): Random seed for determinism. Defaults to 42.
    """
    print(f'Generating realistic deterministic synthetic GTC for {num_snps} loci (Seed: {seed})...')
    rng: random.Random = random.Random(seed)  # noqa: S311

    # Official Illumina Tags
    tag_num_snps: int = 1
    tag_sample_name: int = 101
    tag_normalization_transforms: int = 400
    tag_raw_x: int = 1000
    tag_raw_y: int = 1001
    tag_genotypes: int = 1002
    tag_base_calls: int = 1003
    tag_gencall_scores: int = 1004
    tag_b_allele_freqs: int = 1012
    tag_logr_ratios: int = 1013

    def make_array_block(data: bytes) -> bytes:
        """
        Helper to prefix array data with the element count.
        """
        return struct.pack('<i', num_snps) + data

    # 1. Genotypes: 1 byte per locus. 1=AA, 2=AB, 3=BB, 0=NC
    gts: list[int] = rng.choices([1, 2, 3], weights=[0.3, 0.4, 0.3], k=num_snps)
    genotypes_data: bytes = bytes(gts)

    # Base Calls: 2 bytes per locus. We'll just use 'AA', 'AB', 'BB' based on GT
    base_calls_list: list[bytes] = []
    for gt in gts:
        if gt == 1:
            base_calls_list.append(b'AA')
        elif gt == 2:  # noqa: PLR2004
            base_calls_list.append(b'AB')
        else:
            base_calls_list.append(b'BB')
    base_calls_data: bytes = b''.join(base_calls_list)

    # Scores: float32. High confidence
    scores_data: bytes = struct.pack(f'<{num_snps}f', *([0.98 + rng.random() * 0.02 for _ in range(num_snps)]))

    # BAF and LRR: float32
    baf_list: list[float] = []
    lrr_list: list[float] = []
    for gt in gts:
        # Add some noise to make it look real
        noise: float = rng.gauss(0, 0.02)
        if gt == 1:
            baf_list.append(max(0.0, min(1.0, 0.0 + abs(noise))))
        elif gt == 2:  # noqa: PLR2004
            baf_list.append(max(0.0, min(1.0, 0.5 + noise)))
        else:
            baf_list.append(max(0.0, min(1.0, 1.0 - abs(noise))))
        lrr_list.append(rng.gauss(0, 0.05))

    baf_data: bytes = struct.pack(f'<{num_snps}f', *baf_list)
    lrr_data: bytes = struct.pack(f'<{num_snps}f', *lrr_list)

    # Intensities: uint16  # noqa: ERA001
    raw_x_data: bytes = struct.pack(f'<{num_snps}H', *([1000 + int(rng.gauss(0, 100)) for _ in range(num_snps)]))
    raw_y_data: bytes = struct.pack(f'<{num_snps}H', *([1000 + int(rng.gauss(0, 100)) for _ in range(num_snps)]))

    # 2. Normalization Transforms (dummy but required structure)
    num_xforms: int = 1
    # version, offset_x, offset_y, scale_x, scale_y, shear, theta, cvx, cvy, nn12, rr12, taa, tbb
    # 1 int32 + 12 float32
    identity_xform: bytes = struct.pack('<i', 1) + struct.pack(
        '<12f', 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )
    xforms_data: bytes = struct.pack('<i', num_xforms) + identity_xform * num_xforms

    # 3. Sample Name
    name_bytes: bytes = sample_name.encode('ascii')
    sample_name_data: bytes = leb128_encode(len(name_bytes)) + name_bytes

    # Collect all tags
    tags: list[tuple[int, bytes, bool]] = [
        (tag_num_snps, struct.pack('<i', num_snps), False),
        (tag_sample_name, sample_name_data, False),
        (tag_normalization_transforms, xforms_data, False),
        (tag_raw_x, make_array_block(raw_x_data), False),
        (tag_raw_y, make_array_block(raw_y_data), False),
        (tag_genotypes, make_array_block(genotypes_data), False),
        (tag_base_calls, make_array_block(base_calls_data), False),
        (tag_gencall_scores, make_array_block(scores_data), False),
        (tag_b_allele_freqs, make_array_block(baf_data), False),
        (tag_logr_ratios, make_array_block(lrr_data), False),
    ]

    num_tags: int = len(tags)
    header_size: int = 8
    toc_size: int = num_tags * 6
    current_offset: int = header_size + toc_size

    toc_entries: list[bytes] = []
    data_blocks: list[bytes] = []

    for identifier, data, _ in tags:
        toc_entries.append(struct.pack('<Hi', identifier, current_offset))
        data_blocks.append(data)
        current_offset += len(data)

    with open(output_path, 'wb') as f:
        # GTC Header: 'GTC', version 3, num_tags
        f.write(b'GTC')
        f.write(struct.pack('<bi', 3, num_tags))
        for entry in toc_entries:
            f.write(entry)
        for block in data_blocks:
            f.write(block)

    print(f'Successfully generated realistic synthetic GTC: {output_path}')


def main() -> None:
    """
    CLI main for GTC generator.
    """
    parser = argparse.ArgumentParser(description='Generate synthetic GTC')
    parser.add_argument('output', type=str, help='Output GTC path')
    parser.add_argument('--num', type=int, default=1000, help='Number of SNPs')
    parser.add_argument('--name', type=str, default='GDA-8v1-0_D2.bpm', help='Sample/Manifest name')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')

    args: argparse.Namespace = parser.parse_args()
    generate_gtc(args.output, args.num, args.name, args.seed)


if __name__ == '__main__':
    main()
