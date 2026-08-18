#!/usr/bin/env python3
"""Verify the plink2 AVX2 chrX --hwe fix against the production COH15003 merged pgen.

Queues one Hail Batch job that reads the merged fileset produced by the certified
bootstrap aggregate (ExportCohortDatasets, version 1) and runs the exact SnpQcReport
--hwe command with both the AVX2 and generic plink2 builds under a timeout, then
compares the resulting snplists. This is the command that hung in the pre-20260818
AVX2 builds; RESULT markers in the job log show OK/HANG per binary and MATCH/MISMATCH
for the snplist comparison. Nothing is written to the dataset bucket or metamist.

Run via analysis-runner (full access: the job reads from the main bucket):

    analysis-runner \\
        --dataset ourdna \\
        --access-level full \\
        --image australia-southeast1-docker.pkg.dev/cpg-common/images/popgen_genotyping:latest \\
        --output-dir popgen_genotyping \\
        --description 'plink2 AVX2 chrX hwe fix verification on prod COH15003' \\
        test/scripts/plink2_avx2_hwe_prod_verify.py
"""

from cpg_utils.hail_batch import get_batch

PLINK_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/plink:1.9-20250819-PLINK-2.0-20260818-1'
MERGED_PREFIX = 'gs://cpg-ourdna-main/genotypingarray/popgen_genotyping/ExportCohortDatasets/1/COH15003'
HWE_ARGS = '1e-05 0.001 midp keep-fewhet'
TIMEOUT_SECONDS = 900

b = get_batch('plink2 AVX2 chrX --hwe fix verification on prod COH15003')
j = b.new_bash_job('plink2_avx2_hwe_prod_verify')
j.image(PLINK_IMAGE)
j.cpu(4)
j.storage('20G')

merged_pgen = b.read_input_group(
    pgen=f'{MERGED_PREFIX}.pgen',
    pvar=f'{MERGED_PREFIX}.pvar',
    psam=f'{MERGED_PREFIX}.psam',
)

j.command(
    f"""
set -uxo pipefail
cd $BATCH_TMPDIR

for variant in plink2_intel_avx2 plink2_generic; do
    binary=/usr/local/bin/variants/$variant
    $binary --version
    timeout {TIMEOUT_SECONDS} $binary --pfile {merged_pgen} \\
        --output-chr chrM \\
        --hwe {HWE_ARGS} \\
        --write-snplist \\
        --out hwe_$variant
    status=$?
    if [ $status -eq 124 ]; then
        echo "RESULT $variant: HANG (killed after {TIMEOUT_SECONDS}s)"
    elif [ $status -ne 0 ]; then
        echo "RESULT $variant: CRASH (exit $status)"
    else
        echo "RESULT $variant: OK"
    fi
done

if cmp -s hwe_plink2_intel_avx2.snplist hwe_plink2_generic.snplist; then
    echo "RESULT snplist: MATCH ($(wc -l < hwe_plink2_generic.snplist) variants)"
else
    echo "RESULT snplist: MISMATCH"
    diff hwe_plink2_intel_avx2.snplist hwe_plink2_generic.snplist | head -20
fi
"""
)

b.run(wait=False)
