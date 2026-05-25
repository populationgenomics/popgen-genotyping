"""
Job logic for merging the cohort PLINK 1.9 dataset with an external reference panel.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def run_merge_with_reference_panel(
    cohort_plink_paths: dict[str, str],
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
    Merge the cohort PLINK 1.9 fileset with an external reference panel.

    The cohort is first normalized against the supplied FASTA reference: REF is
    set authoritatively from the FASTA (so A1=ALT/A2=REF regardless of any
    upstream allele-frequency swap), variant IDs are re-derived as
    `chr:pos:ref:alt`, and contigs are re-prefixed with `chr` to match the
    reference panel's convention. Conventions are then asserted against config
    expectations; the merge runs with `plink --bmerge --keep-allele-order` and
    retries once after dropping any sites in `.missnp`. Final output is PLINK2.

    Args:
        cohort_plink_paths (dict[str, str]): Cloud paths for the merged-cohort
            PLINK 1.9 fileset (`bed`, `bim`, `fam`).
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
        bed=cohort_plink_paths['bed'],
        bim=cohort_plink_paths['bim'],
        fam=cohort_plink_paths['fam'],
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
    # NormalizeCohort, ValidateAgainstExpectations, Merge (+ missnp retry),
    # ConvertToPlink2, Stats - all in one BashJob with `set -e`.
    j.command(
        f"""
        set -ex

        # ---- NormalizeCohort -------------------------------------------------
        # `--ref-from-fa force` makes the FASTA authoritative for REF at every
        # site, neutralizing any upstream A1/A2 frequency swaps (e.g. the
        # PLINK 1.9 `--make-bed` minor-allele rotation) before the merge.
        # `--set-all-var-ids` re-derives every ID from the FASTA-anchored
        # orientation, so stale IDs from a frequency-swapped upstream are
        # overwritten with the correct `chr:pos:ref:alt`.
        plink2 --bfile {cohort_input} \\
            --fa {fasta_file.base} --ref-from-fa force \\
            --set-all-var-ids '@:#:$r:$a' \\
            --new-id-max-allele-len 100 \\
            --max-alleles 2 \\
            --output-chr chrM \\
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

        # ---- Merge (first attempt) -------------------------------------------
        # Allow a single failure: PLINK 1.9 returns non-zero on .missnp.
        set +e
        plink --bfile normalized_cohort \\
            --bmerge {reference_input} \\
            --allow-extra-chr \\
            --keep-allele-order \\
            --make-bed --out first_merge
        first_rc=$?
        set -e

        if [ "$first_rc" -ne 0 ]; then
            if [ ! -s first_merge.missnp ]; then
                echo "merge step failed without producing .missnp; check first_merge.log" >&2
                exit "$first_rc"
            fi
            # ---- Drop .missnp on both sides, retry -----------------------------
            plink --bfile normalized_cohort --exclude first_merge.missnp \\
                --allow-extra-chr --keep-allele-order \\
                --make-bed --out cohort_clean
            plink --bfile {reference_input} --exclude first_merge.missnp \\
                --allow-extra-chr --keep-allele-order \\
                --make-bed --out reference_clean
            plink --bfile cohort_clean \\
                --bmerge reference_clean \\
                --allow-extra-chr --keep-allele-order \\
                --make-bed --out final_merge
        else
            mv first_merge.bed final_merge.bed
            mv first_merge.bim final_merge.bim
            mv first_merge.fam final_merge.fam
            mv first_merge.log final_merge.log
        fi

        # ---- ConvertToPlink2 -------------------------------------------------
        plink2 --bfile final_merge --allow-extra-chr --make-pgen --out {j.merged_pgen}

        # ---- Stats -----------------------------------------------------------
        cohort_n=$(wc -l < normalized_cohort.bim)
        ref_n=$(wc -l < {reference_input.bim})
        if [ -s first_merge.missnp ]; then
            missnp_n=$(wc -l < first_merge.missnp)
        else
            missnp_n=0
        fi
        final_n=$(wc -l < final_merge.bim)
        {{
            printf 'metric\\tvalue\\n'
            printf 'cohort_variants_pre_merge\\t%s\\n' "$cohort_n"
            printf 'reference_variants_pre_merge\\t%s\\n' "$ref_n"
            printf 'missnp_dropped\\t%s\\n' "$missnp_n"
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
