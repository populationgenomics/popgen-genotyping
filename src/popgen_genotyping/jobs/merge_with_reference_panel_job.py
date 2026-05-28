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
    drop_strand_ambiguous: bool = True,
    job_name: str = 'merge_with_reference_panel',
) -> BashJob:
    """
    Merge a cohort PLINK2 aggregate with an external reference panel.

    Cohort PGEN is round-tripped to PLINK 1.9, then FASTA-anchored
    (`--ref-from-fa force`), restricted to bi-allelic ACGT SNPs, given
    `chr:pos:ref:alt` variant IDs, and stripped of any position with more
    than one record. When `drop_strand_ambiguous` is true (default),
    strand-ambiguous ({A,T} / {C,G}) sites are also dropped; disabling
    that emits a runtime warning to stderr. Contig style and variant-ID
    pattern are asserted against config on both sides. The intersect is
    computed by ID hash across the two BIMs; both sides are pre-filtered
    with `plink --extract` and joined with a single `plink --bmerge`.
    Final output is PLINK2.

    Args:
        cohort_pgen_paths (dict[str, str]): Cohort PLINK2 fileset paths
            (`pgen`, `pvar`, `psam`).
        reference_panel_paths (dict[str, str]): Reference panel PLINK 1.9
            fileset paths (`bed`, `bim`, `fam`).
        fasta_ref_path (str): GRCh38 FASTA path; `.fai` expected at
            `<fasta_ref_path>.fai`.
        expected_contig_style (str): `'with_chr'` or `'no_chr'`. Asserted on
            both BIMs.
        expected_variant_id_pattern (str): Extended regex (egrep) matching
            permitted variant IDs. Asserted on the first 100 rows of each side.
        output_pgen_prefix (str): Cloud prefix for the merged PLINK2 fileset
            (writes `.pgen`, `.pvar`, `.psam`).
        output_log_path (str): Cloud path for the captured merge log.
        output_stats_path (str): Cloud path for the variant-intersection stats TSV.
        drop_strand_ambiguous (bool): When True (default), drop {A,T} / {C,G}
            sites before intersection. When False, retain them and emit a
            stderr warning — only safe when both sides share a strand
            convention.
        job_name (str): Name for the Hail Batch job.

    Returns:
        BashJob: The queued Hail Batch job.
    """
    strand_filter_enabled: str = 'true' if drop_strand_ambiguous else 'false'
    b = get_batch()
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'merge_with_reference_panel'],
        image=config_retrieve(['workflow', 'plink_image']),
        default_cpu=4,
        default_storage='100G',
    )

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

    j.declare_resource_group(
        merged_pgen={
            'pgen': '{root}.pgen',
            'pvar': '{root}.pvar',
            'psam': '{root}.psam',
        },
        merge_log={'log': '{root}.log'},
        stats={'tsv': '{root}.tsv'},
    )

    j.command(
        f"""
        set -ex

        # ---- ConvertCohortToPlink1 -------------------------------------------
        # PGEN → PLINK 1.9. The rest of the merge runs in PLINK 1.9, and the
        # next step re-anchors REF/ALT against the FASTA, so no need to
        # preserve PGEN's native orientation through the round-trip.
        plink2 --pfile {cohort_input} \\
            --allow-extra-chr \\
            --make-bed --out cohort_plink1

        # ---- NormalizeCohort -------------------------------------------------
        # Pass 1: anchor REF against the FASTA, restrict to bi-allelic ACGT
        # SNPs, re-derive variant IDs as chr:pos:ref:alt. `--ref-from-fa force`
        # makes the FASTA authoritative for REF regardless of the cohort BIM's
        # A1/A2 orientation; `--set-all-var-ids` rewrites IDs from that
        # FASTA-anchored orientation.
        plink2 --bfile cohort_plink1 \\
            --fa {fasta_file.base} --ref-from-fa force \\
            --set-all-var-ids '@:#:$r:$a' \\
            --new-id-max-allele-len 100 \\
            --snps-only just-acgt \\
            --max-alleles 2 \\
            --output-chr chrM \\
            --allow-extra-chr \\
            --make-bed --out normalized_cohort_pre_dedup

        # Pass 2: drop variants that would corrupt the merge.
        #
        # 1. Every variant at any position with more than one record.
        #    `--rm-dup` keys off variant IDs and would miss e.g. (A,C) and
        #    (A,T) at the same coordinate, both of which break the merge
        #    step. awk keys off chr+pos and emits every ID at any
        #    duplicated position.
        # 2. Strand-ambiguous sites: {{A,T}} and {{C,G}} have identical
        #    REF/ALT under strand-flip, so FASTA-anchoring cannot resolve
        #    whether the cohort and the reference panel called them on the
        #    same strand. Dropping them is the safe default; the toggle is
        #    for the case where both sides are known to share a strand
        #    convention.
        awk 'NR==FNR {{count[$1"\\t"$4]++; next}} count[$1"\\t"$4] > 1 {{print $2}}' \\
            normalized_cohort_pre_dedup.bim normalized_cohort_pre_dedup.bim \\
            > duplicate_position_var_ids.txt

        if [ "{strand_filter_enabled}" = "true" ]; then
            awk '
                ($5=="A" && $6=="T") || ($5=="T" && $6=="A") ||
                ($5=="C" && $6=="G") || ($5=="G" && $6=="C") {{print $2}}
            ' normalized_cohort_pre_dedup.bim > strand_ambiguous_var_ids.txt
        else
            echo "WARNING: drop_strand_ambiguous=false — {{A,T}} / {{C,G}} sites retained;" \\
                "merge may produce silent strand-flips if the panel uses a different convention" >&2
            : > strand_ambiguous_var_ids.txt
        fi

        cat duplicate_position_var_ids.txt strand_ambiguous_var_ids.txt \\
            | sort -u > normalize_exclude.txt

        plink2 --bfile normalized_cohort_pre_dedup \\
            --exclude normalize_exclude.txt \\
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
        # Intersect the two BIMs by variant ID. Both sides use the
        # FASTA-anchored chr:pos:ref:alt ID scheme (asserted above), so an
        # ID match implies an allele match; pre-filtering to the intersect
        # avoids allele-set conflicts at merge time.
        awk 'NR==FNR {{ref[$2]=1; next}} ($2 in ref) {{print $2}}' \\
            {reference_input.bim} normalized_cohort.bim > common.ids

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
        cohort_strand_ambig_n=$(wc -l < strand_ambiguous_var_ids.txt)
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
            printf 'cohort_variants_dropped_strand_ambiguous\\t%s\\n' "$cohort_strand_ambig_n"
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

    b.write_output(j.merged_pgen, output_pgen_prefix)
    b.write_output(j.merge_log.log, output_log_path)
    b.write_output(j.stats.tsv, output_stats_path)

    return j
