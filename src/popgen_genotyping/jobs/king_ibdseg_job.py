"""
Job logic for inferring identity-by-descent segments with KING `--ibdseg`.
"""

from __future__ import annotations

from importlib.resources import files
from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve, get_driver_image
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


def run_king_ibdseg(
    bed_path: str,
    bim_path: str,
    fam_path: str,
    output_seg_path: str,
    output_segments_path: str,
    output_seg_x_path: str,
    output_segments_x_path: str,
    output_log_path: str,
    bafregress_paths: list[str],
    smiss_path: str,
    contamination_max: float,
    fmiss_max: float,
    output_excluded_path: str,
    job_name: str = 'king_ibdseg',
) -> list[BashJob]:
    """
    Run KING `--ibdseg` on a merged PLINK 1.9 dataset.

    Three jobs are queued, in dependency order:

    1. ``exclusions``: runs ``scripts/relatedness_exclusions.py`` against the BAFRegress
       output(s) and the merged-cohort ``.smiss`` to select samples whose BAFRegress
       estimate exceeds ``contamination_max`` or whose per-sample ``F_MISS`` exceeds
       ``fmiss_max``. Either condition generates thousands of spurious pairwise
       relationships under KING, which would otherwise contaminate every clean sample's
       relatedness list. Writes an audit TSV (``IID``, ``REASON``, ``BAF_REGRESS``,
       ``F_MISS``; header-only when nothing is excluded) to ``output_excluded_path``, and
       an in-batch plink2 ``--remove`` sample-ID list (``#IID`` header, one IID per line;
       header-only when nothing is excluded).
    2. ``recode``: the merged ``.bim`` carries hg38 ``chr``-prefixed contig names
       (``MergeCohortPlink`` emits ``--output-chr chrM``); KING only parses numeric
       chromosome codes and discards everything else as "other", so a ``chr``-prefixed
       input leaves it with zero autosomal SNPs. This job runs ``plink2 --output-chr 26``
       to recode the fileset to numeric codes (``23``=X, ``24``=Y, ``25``=XY, ``26``=MT)
       that KING recognises, always passing ``--remove`` with the list from step 1 (a
       header-only list removes nothing, so no conditional flag is needed).
    3. ``king``: KING `--ibdseg` against the recoded, exclusion-filtered fileset.

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
        output_seg_path: Cloud path for the autosomal .seg pairwise summary.
        output_segments_path: Cloud path for the autosomal .segments.gz detail.
        output_seg_x_path: Cloud path for the X-chr .seg pairwise summary.
        output_segments_x_path: Cloud path for the X-chr .segments.gz detail.
        output_log_path: Cloud path for the captured KING stdout log.
        bafregress_paths: Cloud paths to the cohort's BAFRegress output file(s).
        smiss_path: Cloud path to the merged-cohort PLINK2 ``.smiss``.
        contamination_max: BAFRegress exclusion threshold (strictly greater than).
        fmiss_max: Per-sample ``F_MISS`` exclusion threshold (strictly greater than).
        output_excluded_path: Cloud path for the excluded-samples audit TSV.
        job_name: Name for the Hail Batch job.

    Returns:
        The queued jobs in dependency order: ``[exclusions, recode, king]``.
    """
    b = get_batch()

    # Job 1: select samples to exclude from KING for contamination or missingness.
    exclusions = register_job(
        batch=b,
        job_name=f'{job_name}_exclusions',
        config_path=['popgen_genotyping', 'king_ibdseg'],
        image=get_driver_image(),
        default_cpu=1,
        default_storage='10G',
    )
    script_path = str(files('popgen_genotyping.scripts').joinpath('relatedness_exclusions.py'))
    script = b.read_input(script_path)
    bafregress_files = [b.read_input(path) for path in bafregress_paths]
    bafregress_bash_args = ' '.join(str(f) for f in bafregress_files)
    smiss_file = b.read_input(smiss_path)

    exclusions.command(
        f"""
        set -ex
        python3 {script} \\
            --bafregress {bafregress_bash_args} \\
            --smiss {smiss_file} \\
            --contamination-max {contamination_max} \\
            --fmiss-max {fmiss_max} \\
            --output-tsv {exclusions.excluded_tsv} \\
            --output-remove-list {exclusions.remove_list}
        """,
    )
    b.write_output(exclusions.excluded_tsv, output_excluded_path)

    # Job 2: recode chromosomes to the numeric codes KING parses, and drop the excluded
    # samples in the same pass. plink2 reads the chr-prefixed input natively (plink 1.9
    # would treat each chr-prefixed code as an unmappable extra contig); --output-chr 26
    # emits 1-22, 23=X, 24=Y, 25=XY, 26=MT. --remove is always passed: a header-only list
    # (nothing excluded) removes no samples, so no conditional flag is needed. The recoded
    # fileset is consumed in-batch by KING and is not written to cloud.
    recode = register_job(
        batch=b,
        job_name=f'{job_name}_recode',
        config_path=['popgen_genotyping', 'king_ibdseg'],
        image=config_retrieve(['workflow', 'plink_image']),
        default_cpu=2,
        default_storage='50G',
    )
    recode.depends_on(exclusions)
    plink_input = b.read_input_group(bed=bed_path, bim=bim_path, fam=fam_path)
    recode.declare_resource_group(
        recoded={
            'bed': '{root}.bed',
            'bim': '{root}.bim',
            'fam': '{root}.fam',
        },
    )
    recode.command(
        f"""
        set -ex
        plink2 --bfile {plink_input} \\
            --allow-extra-chr \\
            --output-chr 26 \\
            --remove {exclusions.remove_list} \\
            --make-bed \\
            --out {recode.recoded}

        input_count=$(wc -l < {plink_input.fam})
        removed_count=$(($(wc -l < {exclusions.remove_list}) - 1))
        recoded_count=$(wc -l < {recode.recoded.fam})
        expected_count=$((input_count - removed_count))
        if [ "$recoded_count" -ne "$expected_count" ]; then
            echo "expected $expected_count samples after --remove ($removed_count excluded), got $recoded_count" >&2
            exit 1
        fi
        """,
    )

    # Job 3: KING --ibdseg against the recoded, exclusion-filtered fileset.
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'king_ibdseg'],
        image=config_retrieve(['workflow', 'king_image']),
        default_cpu=8,
        default_memory='highmem',
        default_storage='50G',
    )
    j.depends_on(recode)

    j.declare_resource_group(
        king_outputs={
            'seg': '{root}.seg',
            'segments_gz': '{root}.segments.gz',
            'seg_x': '{root}X.seg',
            'segments_gz_x': '{root}X.segments.gz',
            'log': '{root}.log',
        },
    )

    # KING omits the .seg / .segments.gz files entirely when no pairs cross the
    # IBD-segment threshold, and omits the X variants entirely when the input
    # has no chrX SNPs. Backfill empty placeholders so Hail Batch's output check
    # succeeds and downstream stages see a well-formed (header-only) table.
    # pipefail: KING exits 1 on a fatal error (e.g. its A1-major-allele check), but the
    # exit status of `king | tee` is tee's, so without it `set -e` never fires and the
    # backfill below would turn a KING failure into a "successful" header-only result.
    j.command(
        f"""
        set -exo pipefail
        king -b {recode.recoded.bed} --ibdseg --degree 3 --cpus $(nproc) \\
            --prefix {j.king_outputs} 2>&1 | tee {j.king_outputs.log}
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

    return [exclusions, recode, j]
