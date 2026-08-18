#!/usr/bin/env python3
"""Generate a synthetic chrX dataset shaped like the cohort that hangs plink2 AVX2 --hwe.

Cohort shape: 2660 samples = 1760 female, 897 male, 3 unknown sex.
Emits a VCF sweeping the female (hom-ref/het/hom-alt) and male (ref/alt) genotype-count
space, including degenerate corners (0 hets, monomorphic, high missingness), plus a
matching .fam-style sex file. Convert and run with plink2:

    python3 make_synthetic_chrx.py
    plink2 --vcf synthetic_chrx.vcf --update-sex synthetic_sex.txt --split-par hg38 \
        --make-pgen --out synthetic_chrx
    plink2 --pfile synthetic_chrx --output-chr chrM --hwe 1e-05 0.001 --write-snplist \
        --out synthetic_hwe
"""

N_FEMALE = 1760
N_MALE = 897
N_UNKNOWN = 3

samples = (
    [f'F{i:04d}' for i in range(N_FEMALE)] + [f'M{i:04d}' for i in range(N_MALE)] + [f'U{i}' for i in range(N_UNKNOWN)]
)

with open('synthetic_sex.txt', 'w') as f:
    f.write('#IID\tSEX\n')
    for s in samples:
        sex = '2' if s.startswith('F') else '1' if s.startswith('M') else 'NA'
        f.write(f'{s}\t{sex}\n')


def genotypes(f_het, f_homalt, f_miss, m_alt, m_miss, u_gt):
    """Genotype columns for one variant: counts for females, males, then unknowns."""
    cols = []
    f_homref = N_FEMALE - f_het - f_homalt - f_miss
    assert f_homref >= 0
    cols += ['0/0'] * f_homref + ['0/1'] * f_het + ['1/1'] * f_homalt + ['./.'] * f_miss
    m_ref = N_MALE - m_alt - m_miss
    assert m_ref >= 0
    cols += ['0'] * m_ref + ['1'] * m_alt + ['.'] * m_miss
    cols += [u_gt] * N_UNKNOWN
    return cols


variants = []
pos = 3_000_000  # inside hg38 chrX non-PAR
# Sweep: het counts including 0 and near-boundary values, hom-alt from monomorphic-ref
# to monomorphic-alt, male alt counts across the range, sprinkled missingness, and
# unknown-sex genotypes as hom, het, haploid and missing.
f_het_vals = [0, 1, 2, 3, 5, 10, 50, 100, 400, 879, 880, 1000, 1500, 1755, 1760]
f_homalt_fracs = [0.0, 0.01, 0.1, 0.5, 0.9]
m_alt_fracs = [0.0, 0.01, 0.5, 0.99, 1.0]
u_gts = ['0/0', '0/1', '1', './.']
miss_pairs = [(0, 0), (35, 18), (176, 90)]  # (female missing, male missing) ~0/2/10%

for f_het in f_het_vals:
    for haf in f_homalt_fracs:
        for maf in m_alt_fracs:
            for f_miss, m_miss in miss_pairs:
                for u_gt in u_gts:
                    f_homalt = int((N_FEMALE - f_het - f_miss) * haf)
                    if f_het + f_homalt + f_miss > N_FEMALE:
                        continue
                    m_alt = int((N_MALE - m_miss) * maf)
                    variants.append((f_het, f_homalt, f_miss, m_alt, m_miss, u_gt))

with open('synthetic_chrx.vcf', 'w') as f:
    f.write('##fileformat=VCFv4.2\n')
    f.write('##contig=<ID=chrX,length=156040895>\n')
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples) + '\n')
    for i, (f_het, f_homalt, f_miss, m_alt, m_miss, u_gt) in enumerate(variants):
        cols = genotypes(f_het, f_homalt, f_miss, m_alt, m_miss, u_gt)
        f.write(
            f'chrX\t{pos + i * 100}\tvar{i}_fh{f_het}_fa{f_homalt}_fm{f_miss}'
            f'_ma{m_alt}_mm{m_miss}_u{u_gt.replace("/", "")}\tA\tG\t.\t.\t.\tGT\t' + '\t'.join(cols) + '\n'
        )

print(f'{len(variants)} variants, {len(samples)} samples written to synthetic_chrx.vcf')
