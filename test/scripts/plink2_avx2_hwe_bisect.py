#!/usr/bin/env python3
"""Bisect the plink2 AVX2 chrX --hwe hang to the exact trigger variant(s).

Queues one Hail Batch job (AVX2 worker, plink image) that generates the synthetic
dataset from test/scripts/make_synthetic_chrx.py, scans it in 100-variant position
windows with the AVX2 binary under a short timeout, then tests every variant of each
hanging window individually. Variant IDs encode their genotype counts
(var<i>_fh<female het>_fa<female homalt>_fm<female missing>_ma<male alt>_mm<male
missing>_u<unknown-sex GT>), so the HANG lines in the job log name the trigger count
combinations directly.

Run via analysis-runner (synthetic data only, test access):

    analysis-runner \\
        --dataset ourdna \\
        --access-level test \\
        --image australia-southeast1-docker.pkg.dev/cpg-common/images/popgen_genotyping:latest \\
        --output-dir popgen_genotyping \\
        --description 'plink2 AVX2 chrX hwe bisect' \\
        test/scripts/plink2_avx2_hwe_bisect.py
"""

from pathlib import Path

from cpg_utils.hail_batch import get_batch

PLINK_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/plink:1.9-20250819-PLINK-2.0-20260808-2'
CHUNK_TIMEOUT = 15
SINGLE_TIMEOUT = 5
BASE_BP = 3_000_000
STEP_BP = 100
N_VARIANTS = 4100
MAX_HANGING_CHUNKS = 6  # per-variant scan cap; remaining hanging chunks are only listed

generator_source = (Path(__file__).parent / 'make_synthetic_chrx.py').read_text()

b = get_batch('plink2 AVX2 chrX --hwe bisect')
j = b.new_bash_job('plink2_avx2_hwe_bisect')
j.image(PLINK_IMAGE)
j.cpu(4)
j.storage('10G')

j.command(
    f"""
set -u
cd $BATCH_TMPDIR

cat <<'GENERATOR_EOF' > make_synthetic_chrx.py
{generator_source}
GENERATOR_EOF

python3 make_synthetic_chrx.py
plink2 --vcf synthetic_chrx.vcf --update-sex synthetic_sex.txt --split-par hg38 \\
    --make-pgen --out synthetic_chrx > /dev/null 2>&1

AVX2=/usr/local/bin/variants/plink2_intel_avx2
$AVX2 --version

scan() {{  # scan <timeout> <from_bp> <to_bp>; echoes exit status
    timeout "$1" $AVX2 --pfile synthetic_chrx --chr X --from-bp "$2" --to-bp "$3" \\
        --hwe 1e-05 0.001 midp keep-fewhet --write-snplist --out scan_out > /dev/null 2>&1
    echo $?
}}

echo "=== chunk scan ({CHUNK_TIMEOUT}s timeout per 100-variant window) ==="
hanging_chunks=""
i=0
while [ $i -lt {N_VARIANTS} ]; do
    end=$((i + 99)); [ $end -ge {N_VARIANTS} ] && end=$(({N_VARIANTS} - 1))
    bp1=$(({BASE_BP} + i * {STEP_BP})); bp2=$(({BASE_BP} + end * {STEP_BP}))
    status=$(scan {CHUNK_TIMEOUT} $bp1 $bp2)
    # Exit 13 is plink2's "no variants remaining after main filters" — a normal
    # outcome when every variant in the window fails the HWE filter, not a fault.
    if [ "$status" = "124" ]; then
        echo "CHUNK HANG: variants $i-$end"
        hanging_chunks="$hanging_chunks $i"
    elif [ "$status" != "0" ] && [ "$status" != "13" ]; then
        echo "CHUNK CRASH (exit $status): variants $i-$end"
        hanging_chunks="$hanging_chunks $i"
    fi
    i=$((i + 100))
done
echo "hanging chunks:$hanging_chunks"

echo "=== per-variant scan ({SINGLE_TIMEOUT}s timeout) ==="
scanned=0
for chunk in $hanging_chunks; do
    if [ $scanned -ge {MAX_HANGING_CHUNKS} ]; then
        echo "chunk $chunk skipped (MAX_HANGING_CHUNKS reached)"
        continue
    fi
    scanned=$((scanned + 1))
    end=$((chunk + 99)); [ $end -ge {N_VARIANTS} ] && end=$(({N_VARIANTS} - 1))
    i=$chunk
    while [ $i -le $end ]; do
        bp=$(({BASE_BP} + i * {STEP_BP}))
        status=$(scan {SINGLE_TIMEOUT} $bp $bp)
        if [ "$status" != "0" ] && [ "$status" != "13" ]; then
            vid=$(awk -v p=$bp '$2 == p {{print $3}}' synthetic_chrx.pvar)
            [ "$status" = "124" ] && echo "HANG: $vid" || echo "CRASH (exit $status): $vid"
        fi
        i=$((i + 1))
    done
done
echo "=== bisect done ==="
"""
)

b.run(wait=False)
