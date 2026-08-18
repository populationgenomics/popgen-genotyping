#!/usr/bin/env python3
"""Queue one Hail Batch job that tests the plink2 AVX2 chrX --hwe hang on real AVX2 hardware.

The job runs in the plink image on a batch worker (which has AVX2), generates the
synthetic chrX dataset from test/scripts/make_synthetic_chrx.py entirely inside the job,
and runs the --hwe command under a timeout with both the auto-selected (AVX2) and the
generic plink2 binaries. Exit markers in the job log show HANG/OK per binary.

Run via analysis-runner (test access is fine; the job touches no cohort data):

    analysis-runner \\
        --dataset ourdna \\
        --access-level test \\
        --image australia-southeast1-docker.pkg.dev/cpg-common/images/popgen_genotyping:latest \\
        --output-dir popgen_genotyping \\
        --description 'plink2 AVX2 chrX hwe reproducer' \\
        test/scripts/plink2_avx2_hwe_repro.py
"""

from pathlib import Path

from cpg_utils.hail_batch import get_batch

PLINK_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/plink:1.9-20250819-PLINK-2.0-20260808-1'
TIMEOUT_SECONDS = 300

generator_source = (Path(__file__).parent / 'make_synthetic_chrx.py').read_text()

b = get_batch('plink2 AVX2 chrX --hwe reproducer')
j = b.new_bash_job('plink2_avx2_hwe_repro')
j.image(PLINK_IMAGE)
j.cpu(4)
j.storage('10G')

j.command(
    f"""
set -uxo pipefail
cd $BATCH_TMPDIR

cat <<'GENERATOR_EOF' > make_synthetic_chrx.py
{generator_source}
GENERATOR_EOF

python3 make_synthetic_chrx.py
plink2 --vcf synthetic_chrx.vcf --update-sex synthetic_sex.txt --split-par hg38 \\
    --make-pgen --out synthetic_chrx

grep -o "avx2" /proc/cpuinfo | head -1

for variant in plink2_intel_avx2 plink2_generic; do
    binary=/usr/local/bin/variants/$variant
    $binary --version
    timeout {TIMEOUT_SECONDS} $binary --pfile synthetic_chrx --output-chr chrM \\
        --hwe 1e-05 0.001 midp keep-fewhet --write-snplist --out hwe_$variant
    status=$?
    if [ $status -eq 124 ]; then
        echo "RESULT $variant: HANG (killed after {TIMEOUT_SECONDS}s)"
    elif [ $status -ne 0 ]; then
        echo "RESULT $variant: CRASH (exit $status)"
    else
        echo "RESULT $variant: OK"
    fi
done
"""
)

b.run(wait=False)
