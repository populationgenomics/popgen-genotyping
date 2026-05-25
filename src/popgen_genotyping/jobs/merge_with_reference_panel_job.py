"""
Job logic for merging a cohort PLINK2 aggregate with an external reference panel.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def run_merge_with_reference_panel(
    cohort_pgen_paths: dict[str, str],
    reference_panel_paths: dict[str, str],
    fasta_ref_path: str,
    expected_contig_style: str,
    expected_variant_id_pattern: str,
    output_pgen_prefix: str,
    output_log_path: str,
    output_stats_path: str,
    job_name: str = 'merge_with_reference_panel',
) -> BashJob:
    """
    Merge a cohort PLINK2 aggregate with an external reference panel.

    The cohort PGEN/PVAR/PSAM is first round-tripped to PLINK 1.9 (the merge
    pipeline operates in PLINK 1.9 throughout), then normalized against the
    supplied FASTA reference: REF is set authoritatively from the FASTA (so
    A1=ALT/A2=REF regardless of any upstream allele-frequency swap), variant
    IDs are re-derived as `chr:pos:ref:alt`, and contigs are re-prefixed with
    `chr` to match the reference panel's convention. Conventions are then
    asserted against config expectations. The per-locus intersect is computed
    inline as an awk ID hash on the two BIMs; both sides are pre-filtered to
    that intersect via `plink --extract` before a single `plink --bmerge
    --keep-allele-order`. Because both BIMs use the FASTA-anchored
    `chr:pos:ref:alt` ID scheme, an ID match implies an allele match — the
    intersect eliminates the `.missnp` class of conflicts. Final output is
    PLINK2.

    Args:
        cohort_pgen_paths (dict[str, str]): Cloud paths for the merged-cohort
            PLINK2 fileset (`pgen`, `pvar`, `psam`). Typically the output of a
            prior `ExportCohortDatasets` run, pointed at by config.
        reference_panel_paths (dict[str, str]): Cloud paths for the reference
            panel PLINK 1.9 fileset (`bed`, `bim`, `fam`).
        fasta_ref_path (str): Cloud path to the GRCh38 FASTA (its `.fai` is
            expected at `<fasta_ref_path>.fai`).
        expected_contig_style (str): `'with_chr'` or `'no_chr'`. Both reference
            and post-normalization cohort BIMs are asserted to match.
        expected_variant_id_pattern (str): Extended regex (egrep) matching
            permitted variant IDs. Asserted on the first 100 rows of each side.
        output_pgen_prefix (str): Cloud prefix for the merged PLINK2 fileset
            (writes `.pgen`, `.pvar`, `.psam`).
        output_log_path (str): Cloud path for the captured merge log.
        output_stats_path (str): Cloud path for the variant-intersection stats TSV.
        job_name (str): Name for the Hail Batch job.

    Returns:
        BashJob: The queued Hail Batch job.
    """
    b = get_batch()
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'merge_with_reference_panel'],
        image=config_retrieve(['workflow', 'plink_image']),
        default_cpu=4,
        default_storage='100G',
    )

    # 1. Stage inputs
    cohort_input = b.read_input_group(
        pgen=cohort_pgen_paths['pgen'],
        pvar=cohort_pgen_paths['pvar'],
        psam=cohort_pgen_paths['psam'],
    )
    reference_input = b.read_input_group(
        bed=reference_panel_paths['bed'],
        bim=reference_panel_paths['bim'],
        fam=reference_panel_paths['fam'],
    )
    fasta_file = b.read_input_group(
        base=fasta_ref_path,
        fai=fasta_ref_path + '.fai',
    )

    # 2. Declare output resource groups
    j.declare_resource_group(
        merged_pgen={
            'pgen': '{root}.pgen',
            'pvar': '{root}.pvar',
            'psam': '{root}.psam',
        },
        merge_log={'log': '{root}.log'},
        stats={'tsv': '{root}.tsv'},
    )

    # 3. Build the combined command
    # ConvertCohortToPlink1, NormalizeCohort, ValidateAgainstExpectations,
    # ComputeIntersect, Merge, ConvertToPlink2, Stats - all in one BashJob
    # with `set -e`.
    j.command(
        f"""
        set -ex

        # ---- ConvertCohortToPlink1 -------------------------------------------
        # Round-trip the cohort PGEN/PVAR/PSAM to PLINK 1.9 BED/BIM/FAM. The
        # rest of the merge pipeline operates in PLINK 1.9 (the reference panel
        # is distributed as PLINK 1.9), and the next step re-anchors REF/ALT
        # against the FASTA anyway — so we don't need to preserve PGEN's native
        # orientation through the round-trip.
        plink2 --pfile {cohort_input} \\
            --allow-extra-chr \\
            --make-bed --out cohort_plink1

        # ---- NormalizeCohort -------------------------------------------------
        # Pass 1: anchor REF against the FASTA, restrict to bi-allelic ACGT SNPs,
        # re-derive variant IDs as chr:pos:ref:alt, re-prefix contigs.
        #
        # `--ref-from-fa force` makes the FASTA authoritative for REF at every
        # site, neutralizing any upstream A1/A2 frequency swaps (e.g. the
        # PLINK 1.9 `--make-bed` minor-allele rotation) before the merge.
        # `--set-all-var-ids` re-derives every ID from the FASTA-anchored
        # orientation, so stale IDs from a frequency-swapped upstream are
        # overwritten with the correct `chr:pos:ref:alt`.
        # `--snps-only just-acgt` drops indels and non-ACGT alleles; the
        # reference panel is bi-allelic SNPs only, so anything else is noise.
        plink2 --bfile cohort_plink1 \\
            --fa {fasta_file.base} --ref-from-fa force \\
            --set-all-var-ids '@:#:$r:$a' \\
            --new-id-max-allele-len 100 \\
            --snps-only just-acgt \\
            --max-alleles 2 \\
            --output-chr chrM \\
            --allow-extra-chr \\
            --make-bed --out normalized_cohort_pre_dedup

        # Pass 2: drop every variant at any position with more than one record.
        # `--rm-dup` keys off variant IDs and would not catch e.g. (A,C) and
        # (A,T) at the same coordinate, which still confuses the merge step.
        # awk identifies positions with multiple records and emits the IDs of
        # *all* of them so the next `--exclude` removes the whole group
        # (keeping none, per the merge contract).
        awk 'NR==FNR {{count[$1"\\t"$4]++; next}} count[$1"\\t"$4] > 1 {{print $2}}' \\
            normalized_cohort_pre_dedup.bim normalized_cohort_pre_dedup.bim \\
            > duplicate_position_var_ids.txt

        plink2 --bfile normalized_cohort_pre_dedup \\
            --exclude duplicate_position_var_ids.txt \\
            --allow-extra-chr \\
            --make-bed --out normalized_cohort

        # ---- ValidateAgainstExpectations -------------------------------------
        cohort_style=$(awk '{{print ($1 ~ /^chr/) ? "with_chr" : "no_chr"; exit}}' normalized_cohort.bim)
        ref_style=$(awk '{{print ($1 ~ /^chr/) ? "with_chr" : "no_chr"; exit}}' {reference_input.bim})
        if [ "$cohort_style" != "{expected_contig_style}" ]; then
            echo "normalized cohort contig style '$cohort_style' != expected '{expected_contig_style}'" >&2
            exit 1
        fi
        if [ "$ref_style" != "{expected_contig_style}" ]; then
            echo "reference contig style '$ref_style' != expected '{expected_contig_style}'" >&2
            exit 1
        fi
        if ! awk 'NR<=100 {{print $2}}' normalized_cohort.bim | grep -qE '{expected_variant_id_pattern}'; then
            echo "normalized cohort variant IDs do not match '{expected_variant_id_pattern}'" >&2
            exit 1
        fi
        if ! awk 'NR<=100 {{print $2}}' {reference_input.bim} | grep -qE '{expected_variant_id_pattern}'; then
            echo "reference variant IDs do not match '{expected_variant_id_pattern}'" >&2
            exit 1
        fi

        # ---- ComputeIntersect ------------------------------------------------
        # Build the per-locus extract list as the intersection of variant IDs.
        # Both BIMs use the FASTA-anchored chr:pos:ref:alt ID scheme (asserted
        # above), so an ID match implies an allele match — pre-filtering by
        # this list eliminates the .missnp class of merge conflicts and makes
        # a retry path unnecessary.
        awk 'NR==FNR {{ref[$2]=1; next}} ($2 in ref) {{print $2}}' \\
            {reference_input.bim} normalized_cohort.bim > common.ids

        # Fail fast on a fully-disjoint panel: bmerge on empty filesets would
        # abort downstream with a less informative error.
        if [ ! -s common.ids ]; then
            echo "no variants in common between normalized cohort and reference panel" >&2
            exit 1
        fi

        plink --bfile normalized_cohort --extract common.ids \\
            --allow-extra-chr --keep-allele-order \\
            --make-bed --out cohort_intersect
        plink --bfile {reference_input} --extract common.ids \\
            --allow-extra-chr --keep-allele-order \\
            --make-bed --out reference_intersect

        # ---- Merge -----------------------------------------------------------
        plink --bfile cohort_intersect \\
            --bmerge reference_intersect \\
            --allow-extra-chr \\
            --keep-allele-order \\
            --make-bed --out final_merge

        # ---- ConvertToPlink2 -------------------------------------------------
        plink2 --bfile final_merge --allow-extra-chr --make-pgen --out {j.merged_pgen}

        # ---- Stats -----------------------------------------------------------
        cohort_pre_dedup_n=$(wc -l < normalized_cohort_pre_dedup.bim)
        cohort_dup_position_n=$(wc -l < duplicate_position_var_ids.txt)
        cohort_n=$(wc -l < normalized_cohort.bim)
        ref_n=$(wc -l < {reference_input.bim})
        intersect_n=$(wc -l < common.ids)
        cohort_only_n=$((cohort_n - intersect_n))
        ref_only_n=$((ref_n - intersect_n))
        final_n=$(wc -l < final_merge.bim)
        {{
            printf 'metric\\tvalue\\n'
            printf 'cohort_variants_post_snp_filter\\t%s\\n' "$cohort_pre_dedup_n"
            printf 'cohort_variants_dropped_duplicate_position\\t%s\\n' "$cohort_dup_position_n"
            printf 'cohort_variants_pre_merge\\t%s\\n' "$cohort_n"
            printf 'reference_variants_pre_merge\\t%s\\n' "$ref_n"
            printf 'intersect_variants\\t%s\\n' "$intersect_n"
            printf 'cohort_only_variants\\t%s\\n' "$cohort_only_n"
            printf 'reference_only_variants\\t%s\\n' "$ref_only_n"
            printf 'final_merged_variants\\t%s\\n' "$final_n"
        }} > {j.stats.tsv}

        # ---- Capture log -----------------------------------------------------
        cp final_merge.log {j.merge_log.log}
        """
    )

    # 4. Write outputs back to cloud
    b.write_output(j.merged_pgen, output_pgen_prefix)
    b.write_output(j.merge_log.log, output_log_path)
    b.write_output(j.stats.tsv, output_stats_path)

    return j
