"""
Job logic for inferring identity-by-descent segments with KING `--ibdseg`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


# Header for the autosomal pairwise IBD summary (KING --ibdseg).
_AUTOSOME_SEG_HEADER: str = 'FID1\tID1\tFID2\tID2\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType'

# Header for the X-chromosome pairwise IBD summary. KING records per-sample
# sex and per-individual maximum IBD-segment proportions on chrX, but does not
# emit an InfType column (relationship inference on X depends on the sex
# combination of the pair, not on a single categorical).
_X_SEG_HEADER: str = 'FID1\tID1\tFID2\tID2\tSex1\tSex2\tMaxIBD1\tMaxIBD2\tIBD1Seg\tIBD2Seg\tPropIBD'

# Maximum BAFRegress contamination estimate accepted before a sample is excluded
# from the IBD calculation. The 3% cutoff follows the BAFRegress convention
# (Jun et al. 2012, doi:10.1016/j.ajhg.2012.09.004); samples whose estimate is
# above this, non-numeric, or absent from any BAFRegress output are dropped.
BAFREGRESS_THRESHOLD: float = 0.03


def run_king_ibdseg(
    bed_path: str,
    bim_path: str,
    fam_path: str,
    bafregress_paths: list[str],
    output_seg_path: str,
    output_segments_path: str,
    output_seg_x_path: str,
    output_segments_x_path: str,
    output_log_path: str,
    job_name: str = 'king_ibdseg',
) -> BashJob:
    """
    Run KING `--ibdseg` on a merged PLINK 1.9 dataset.

    KING 2.3.2 writes at the chosen prefix: ``{prefix}.seg`` (autosomal pairwise
    IBD summary), ``{prefix}.segments.gz`` (autosomal per-segment detail) and a
    ``{prefix}allsegs.txt`` informational file (analysis metadata, not
    relatedness). When the input dataset includes chrX SNPs it additionally
    writes ``{prefix}X.seg`` and ``{prefix}X.segments.gz`` (no dot before the
    trailing ``X``). KING streams its log to stdout, so this job tees the run
    into a ``.log`` file captured alongside the relatedness outputs. Refs:
    Manichaikul et al. (2010) doi:10.1093/bioinformatics/btq559; Chen et al.
    (2024).

    Args:
        bed_path: Cloud path to the PLINK 1.9 .bed file.
        bim_path: Cloud path to the PLINK 1.9 .bim file.
        fam_path: Cloud path to the PLINK 1.9 .fam file.
        bafregress_paths: Per-cohort BAFRegress output files. Samples whose
            ``estimate`` column exceeds ``BAFREGRESS_THRESHOLD``, is
            non-numeric, or is absent from every file are dropped from the
            KING run via ``--remove``.
        output_seg_path: Cloud path for the autosomal .seg pairwise summary.
        output_segments_path: Cloud path for the autosomal .segments.gz detail.
        output_seg_x_path: Cloud path for the X-chr .seg pairwise summary.
        output_segments_x_path: Cloud path for the X-chr .segments.gz detail.
        output_log_path: Cloud path for the captured KING stdout log.
        job_name: Name for the Hail Batch job.

    Returns:
        The queued Hail Batch BashJob.
    """
    b = get_batch()
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'king_ibdseg'],
        image=config_retrieve(['workflow', 'king_image']),
        default_cpu=8,
        default_memory='highmem',
        default_storage='50G',
    )

    plink_input = b.read_input_group(bed=bed_path, bim=bim_path, fam=fam_path)
    bafregress_files = [b.read_input(p) for p in bafregress_paths]

    j.declare_resource_group(
        king_outputs={
            'seg': '{root}.seg',
            'segments_gz': '{root}.segments.gz',
            'seg_x': '{root}X.seg',
            'segments_gz_x': '{root}X.segments.gz',
            'log': '{root}.log',
        },
    )

    # Build a KING --remove list of samples that fail the BAFRegress
    # contamination filter. "Fail" means: estimate > BAFREGRESS_THRESHOLD,
    # estimate is non-numeric, or the sample is absent from every BAFRegress
    # file. Built as (all .fam IIDs) \ (IIDs with a numeric estimate <= threshold)
    # so that the absent-from-BAFRegress case is handled by the same set diff.
    #
    # KING omits the .seg / .segments.gz files entirely when no pairs cross the
    # IBD-segment threshold, and omits the X variants entirely when the input
    # has no chrX SNPs. Backfill empty placeholders so Hail Batch's output check
    # succeeds and downstream stages see a well-formed (header-only) table.
    bafregress_inputs_block = ' '.join(bafregress_files) if bafregress_files else ''
    j.command(
        f"""
        set -euo pipefail

        # All FID/IID pairs in the merged PLINK fileset, sorted by IID.
        awk '{{print $1"\\t"$2}}' {plink_input.fam} | sort -k2,2 > all_samples.tsv

        # IIDs with a numeric BAFRegress estimate <= {BAFREGRESS_THRESHOLD}.
        : > good_iids.txt
        for baf in {bafregress_inputs_block}; do
            awk -v thresh={BAFREGRESS_THRESHOLD} '
                NR==1 {{
                    sid_col = 0; est_col = 0
                    for (i=1; i<=NF; i++) {{
                        if ($i == "sample_id") sid_col = i
                        else if ($i == "estimate") est_col = i
                    }}
                    if (sid_col == 0 || est_col == 0) {{
                        print "ERROR: missing sample_id or estimate column in " FILENAME > "/dev/stderr"
                        exit 1
                    }}
                    next
                }}
                $est_col ~ /^-?[0-9]+\\.?[0-9]*([eE][-+]?[0-9]+)?$/ && ($est_col + 0) <= thresh {{
                    print $sid_col
                }}
            ' "$baf" >> good_iids.txt
        done
        sort -u good_iids.txt -o good_iids.txt

        # remove = samples in .fam whose IID is not in good_iids.
        awk 'NR==FNR{{good[$1]=1; next}} !($2 in good) {{print $1"\\t"$2}}' \\
            good_iids.txt all_samples.tsv > remove_samples.tsv

        n_total=$(wc -l < all_samples.tsv)
        n_remove=$(wc -l < remove_samples.tsv)
        echo "BAFRegress contamination filter: excluding $n_remove / $n_total samples" \\
             "(estimate > {BAFREGRESS_THRESHOLD}, non-numeric, or absent)"

        if [ "$n_remove" -gt 0 ]; then
            king -b {plink_input.bed} --ibdseg --degree 3 --cpus $(nproc) \\
                --remove remove_samples.tsv \\
                --prefix {j.king_outputs} 2>&1 | tee {j.king_outputs.log}
        else
            king -b {plink_input.bed} --ibdseg --degree 3 --cpus $(nproc) \\
                --prefix {j.king_outputs} 2>&1 | tee {j.king_outputs.log}
        fi

        if [ ! -s {j.king_outputs.seg} ]; then
            printf '%s\\n' '{_AUTOSOME_SEG_HEADER}' > {j.king_outputs.seg}
        fi
        if [ ! -s {j.king_outputs.segments_gz} ]; then
            : | gzip -c > {j.king_outputs.segments_gz}
        fi
        if [ ! -s {j.king_outputs.seg_x} ]; then
            printf '%s\\n' '{_X_SEG_HEADER}' > {j.king_outputs.seg_x}
        fi
        if [ ! -s {j.king_outputs.segments_gz_x} ]; then
            : | gzip -c > {j.king_outputs.segments_gz_x}
        fi
        """,
    )

    b.write_output(j.king_outputs.seg, output_seg_path)
    b.write_output(j.king_outputs.segments_gz, output_segments_path)
    b.write_output(j.king_outputs.seg_x, output_seg_x_path)
    b.write_output(j.king_outputs.segments_gz_x, output_segments_x_path)
    b.write_output(j.king_outputs.log, output_log_path)

    return j
