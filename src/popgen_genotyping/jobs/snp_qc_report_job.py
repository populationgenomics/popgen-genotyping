"""Hail Batch jobs that produce the per-SNP QC inclusion list.

Three jobs:

1. ``bcftools_image`` extract — pulls ``GenTrain_Score`` and ``Cluster_Sep`` (plus
   CHROM/POS/REF/ALT/ID) out of the references-repo EGT info BCF as a
   header-less TSV.
2. ``plink_image`` merged-set variant metrics — runs ``plink2 --missing`` and
   ``plink2 --hwe ... --write-snplist`` against the merged PGEN to write the
   per-variant F_MISS table and the Hardy-Weinberg pass snplist. Two plink2
   invocations are required because ``--hwe`` is a variant filter and would
   otherwise restrict ``--missing`` to the post-filter variant set.
3. ``driver_image`` runs ``scripts/snp_qc_report.py``, joining the EGT TSV,
   the merged-set ``.vmiss``, and the HWE pass snplist, and applying the
   thresholds declared under
   ``[popgen_genotyping.snp_qc_report.thresholds]``.
"""

from __future__ import annotations

from importlib.resources import files
from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve, get_driver_image
from cpg_utils.hail_batch import get_batch

from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch import Batch
    from hailtop.batch.job import BashJob


def run_snp_qc_report(
    *,
    egt_info_bcf_path: str,
    egt_info_bcf_index_path: str,
    merged_pgen_path: str,
    merged_pvar_path: str,
    merged_psam_path: str,
    gentrain_min: float,
    cluster_sep_min: float,
    fmiss_max: float,
    hwe_p: float,
    hwe_k: float,
    hwe_midp: bool,
    hwe_keep_fewhet: bool,
    output_audit_tsv_path: str,
    output_inclusion_list_path: str,
    output_summary_tsv_path: str,
    job_name: str = 'snp_qc_report',
) -> list[BashJob]:
    """Queue the jobs.

    Args:
        egt_info_bcf_path: Path to the references-repo EGT info BCF.
        egt_info_bcf_index_path: Path to its ``.csi`` index.
        merged_pgen_path: Path to the merged PLINK2 ``.pgen``.
        merged_pvar_path: Path to the merged PLINK2 ``.pvar``.
        merged_psam_path: Path to the merged PLINK2 ``.psam``.
        gentrain_min: Inclusive lower bound for ``GenTrain_Score``.
        cluster_sep_min: Inclusive lower bound for ``Cluster_Sep``.
        fmiss_max: Inclusive upper bound for the merged-set ``F_MISS``.
        hwe_p: ``p`` argument to ``plink2 --hwe`` (per-variant p-value floor;
            effective threshold is ``p · 10^(-n·k)``).
        hwe_k: ``k`` argument to ``plink2 --hwe`` (sample-size scaling exponent).
        hwe_midp: If True, pass the ``midp`` modifier to ``plink2 --hwe``.
        hwe_keep_fewhet: If True, pass the ``keep-fewhet`` modifier so only
            excess-heterozygosity variants are filtered.
        output_audit_tsv_path: Output path for the bgzipped audit TSV.
        output_inclusion_list_path: Output path for the inclusion ``.snplist``
            (passing variant IDs, for ``plink2 --extract``).
        output_summary_tsv_path: Output path for the summary TSV.
        job_name: Display name for the filter job (the extract jobs append
            ``_extract_egt_info`` / ``_hwe`` suffixes).

    Returns:
        Queued jobs in dependency order: ``[extract_cluster_file_info, hwe_snplist, generate_snp_inclusion_list]``.
    """
    b: Batch = get_batch()

    extract_cluster_file_info: BashJob = register_job(
        batch=b,
        job_name=f'{job_name}_extract_egt_info',
        config_path=['popgen_genotyping', 'snp_qc_report', 'extract'],
        image=config_retrieve(['workflow', 'bcftools_image']),
        default_cpu=1,
        default_storage='10G',
    )

    egt_input = b.read_input_group(bcf=egt_info_bcf_path, csi=egt_info_bcf_index_path)
    extract_cluster_file_info.declare_resource_group(egt_info={'tsv': '{root}.egt_info.tsv'})

    egt_info_fields: str = (
        r'%CHROM\t%POS\t%REF\t%ALT\t%ID\t'
        r'%INFO/GenTrain_Score\t%INFO/Cluster_Sep\n'
    )
    extract_cluster_file_info.command(
        f"""
        set -euxo pipefail
        bcftools query \\
            -f '{egt_info_fields}' \\
            {egt_input.bcf} > {extract_cluster_file_info.egt_info.tsv}
        """,
    )

    hwe_snplist: BashJob = register_job(
        batch=b,
        job_name=f'{job_name}_hwe',
        config_path=['popgen_genotyping', 'snp_qc_report', 'hwe'],
        image=config_retrieve(['workflow', 'plink_image']),
        default_cpu=2,
        default_storage='50G',
    )

    merged_pgen = b.read_input_group(
        pgen=merged_pgen_path,
        pvar=merged_pvar_path,
        psam=merged_psam_path,
    )
    hwe_snplist.declare_resource_group(
        missing={'vmiss': '{root}.vmiss'},
        hwe_pass={'snplist': '{root}.snplist'},
    )

    hwe_modifiers: str = ' '.join(m for m, on in (('midp', hwe_midp), ('keep-fewhet', hwe_keep_fewhet)) if on)
    hwe_args: str = f'{hwe_p} {hwe_k}' + (f' {hwe_modifiers}' if hwe_modifiers else '')

    # Two plink2 runs with distinct --out prefixes: --missing writes the F_MISS
    # table (missing.vmiss) over the full merged variant set, and --hwe
    # --write-snplist writes the variants surviving the HWE filter
    # (hwe_pass.snplist). Separate prefixes keep the F_MISS table over the full
    # variant set rather than the HWE-filtered subset.
    hwe_snplist.command(
        f"""
        set -euxo pipefail
        plink2 --pfile {merged_pgen} \\
            --output-chr chrM \\
            --missing \\
            --out {hwe_snplist.missing}
        plink2 --pfile {merged_pgen} \\
            --output-chr chrM \\
            --hwe {hwe_args} \\
            --write-snplist \\
            --out {hwe_snplist.hwe_pass}
        """,
    )

    generate_snp_inclusion_list: BashJob = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'snp_qc_report'],
        image=get_driver_image(),
        default_cpu=2,
        default_memory='standard',
        default_storage='20G',
    )
    generate_snp_inclusion_list.depends_on(extract_cluster_file_info, hwe_snplist)

    script_path: str = str(files('popgen_genotyping.scripts').joinpath('snp_qc_report.py'))
    script_resource = b.read_input(script_path)

    generate_snp_inclusion_list.declare_resource_group(
        outputs={
            'audit_tsv_gz': '{root}.audit.tsv.gz',
            'inclusion_list': '{root}.include.snplist',
            'summary_tsv': '{root}.summary.tsv',
        },
    )

    generate_snp_inclusion_list.command(
        f"""
        set -euxo pipefail
        python3 {script_resource} \\
            --egt-info-tsv {extract_cluster_file_info.egt_info.tsv} \\
            --merged-vmiss {hwe_snplist.missing.vmiss} \\
            --hwe-pass-snplist {hwe_snplist.hwe_pass.snplist} \\
            --gentrain-min {gentrain_min} \\
            --cluster-sep-min {cluster_sep_min} \\
            --fmiss-max {fmiss_max} \\
            --output-audit-tsv {generate_snp_inclusion_list.outputs.audit_tsv_gz} \\
            --output-inclusion-list {generate_snp_inclusion_list.outputs.inclusion_list} \\
            --output-summary-tsv {generate_snp_inclusion_list.outputs.summary_tsv}
        """,
    )

    b.write_output(generate_snp_inclusion_list.outputs.audit_tsv_gz, output_audit_tsv_path)
    b.write_output(generate_snp_inclusion_list.outputs.inclusion_list, output_inclusion_list_path)
    b.write_output(generate_snp_inclusion_list.outputs.summary_tsv, output_summary_tsv_path)

    return [extract_cluster_file_info, hwe_snplist, generate_snp_inclusion_list]
