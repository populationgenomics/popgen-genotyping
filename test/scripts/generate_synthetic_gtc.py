import struct
import random
import os

def leb128_encode(n):
    """LEB128 encode an integer."""
    out = bytearray()
    while True:
        byte = n & 0x7f
        n >>= 7
        if n == 0:
            out.append(byte)
            return bytes(out)
        out.append(byte | 0x80)

def generate_gtc(output_path, num_snps, sample_name, seed=42):
    print(f"Generating realistic deterministic synthetic GTC for {num_snps} loci (Seed: {seed})...")
    rng = random.Random(seed)
    
    # Official Illumina Tags
    TAG_NUM_SNPS = 1
    TAG_SAMPLE_NAME = 101
    TAG_NORMALIZATION_TRANSFORMS = 400
    TAG_RAW_X = 1000
    TAG_RAW_Y = 1001
    TAG_GENOTYPES = 1002
    TAG_BASE_CALLS = 1003
    TAG_GENCALL_SCORES = 1004
    TAG_B_ALLELE_FREQS = 1012
    TAG_LOGR_RATIOS = 1013
    
    def make_array_block(data):
        return struct.pack('<i', num_snps) + data

    # 1. Generate randomized data using seeded rng
    # Genotypes: 1=AA, 2=AB, 3=BB
    gts = rng.choices([1, 2, 3], weights=[0.3, 0.4, 0.3], k=num_snps)
    genotypes_data = bytes(gts)
    
    # Base Calls: 2 bytes per locus. We'll just use 'AA', 'AB', 'BB' based on GT
    base_calls_list = []
    for gt in gts:
        if gt == 1: base_calls_list.append(b'AA')
        elif gt == 2: base_calls_list.append(b'AB')
        else: base_calls_list.append(b'BB')
    base_calls_data = b''.join(base_calls_list)
    
    # Scores: float32. High confidence
    scores_data = struct.pack(f'<{num_snps}f', *([0.98 + rng.random()*0.02 for _ in range(num_snps)]))
    
    # BAF and LRR: float32
    baf_list = []
    lrr_list = []
    for gt in gts:
        # Add some noise to make it look real
        noise = rng.gauss(0, 0.02)
        if gt == 1: baf_list.append(max(0.0, min(1.0, 0.0 + abs(noise))))
        elif gt == 2: baf_list.append(max(0.0, min(1.0, 0.5 + noise)))
        else: baf_list.append(max(0.0, min(1.0, 1.0 - abs(noise))))
        lrr_list.append(rng.gauss(0, 0.05))
        
    baf_data = struct.pack(f'<{num_snps}f', *baf_list)
    lrr_data = struct.pack(f'<{num_snps}f', *lrr_list)
    
    # Intensities: uint16
    raw_x_data = struct.pack(f'<{num_snps}H', *([1000 + int(rng.gauss(0, 100)) for _ in range(num_snps)]))
    raw_y_data = struct.pack(f'<{num_snps}H', *([1000 + int(rng.gauss(0, 100)) for _ in range(num_snps)]))

    # 2. Normalization Transforms (Tag 400)
    # Fix: Identity transform must have scale = 1.0
    num_xforms = 100
    # version, offset_x, offset_y, scale_x, scale_y, shear, theta, cvx, cvy, nn12, rr12, taa, tbb
    # 1 int32 + 12 float32
    identity_xform = struct.pack('<i', 1) + struct.pack('<12f', 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    xforms_data = struct.pack('<i', num_xforms) + identity_xform * num_xforms
    
    # 3. Sample Name
    name_bytes = sample_name.encode('ascii')
    sample_name_block = leb128_encode(len(name_bytes)) + name_bytes

    # 4. Assembly
    tags = [
        (TAG_NUM_SNPS, num_snps, True),
        (TAG_SAMPLE_NAME, sample_name_block, False),
        (TAG_NORMALIZATION_TRANSFORMS, xforms_data, False),
        (TAG_RAW_X, make_array_block(raw_x_data), False),
        (TAG_RAW_Y, make_array_block(raw_y_data), False),
        (TAG_GENOTYPES, make_array_block(genotypes_data), False),
        (TAG_BASE_CALLS, make_array_block(base_calls_data), False),
        (TAG_GENCALL_SCORES, make_array_block(scores_data), False),
        (TAG_B_ALLELE_FREQS, make_array_block(baf_data), False),
        (TAG_LOGR_RATIOS, make_array_block(lrr_data), False),
    ]
    
    num_tags = len(tags)
    header_size = 8
    toc_size = num_tags * 6
    current_offset = header_size + toc_size
    
    toc_entries = []
    data_blocks = []
    for tag_id, content, is_direct in tags:
        if is_direct:
            toc_entries.append(struct.pack('<Hi', tag_id, content))
        else:
            toc_entries.append(struct.pack('<Hi', tag_id, current_offset))
            data_blocks.append(content)
            current_offset += len(content)

    with open(output_path, 'wb') as f:
        f.write(b'gtc\x05')
        f.write(struct.pack('<i', num_tags))
        for entry in toc_entries:
            f.write(entry)
        for block in data_blocks:
            f.write(block)

    print(f"Successfully generated realistic synthetic GTC: {output_path}")

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Generate a synthetic Illumina GTC file.')
    parser.add_argument('output', type=str, help='Output path')
    parser.add_argument('--num', type=int, default=1904599, help='Number of SNPs')
    parser.add_argument('--name', type=str, default='GDA-8v1-0_D2.bpm', help='Sample/Manifest name')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    
    args = parser.parse_args()
    generate_gtc(args.output, args.num, args.name, args.seed)
