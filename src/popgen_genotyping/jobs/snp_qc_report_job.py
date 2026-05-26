"""Hail Batch jobs to produces the per-SNP QC exclusion list.

Two jobs:

1. ``bcftools_image`` extract — pulls ``GenTrain_Score`` and ``Cluster_Sep`` (plus
   CHROM/POS/REF/ALT/ID) out of the references-repo EGT info BCF as a
   header-less TSV.
2. ``driver_image`` runs ``scripts/snp_qc_report.py``, joining the
   EGT TSV with the merged-set ``plink2 --missing`` output and applying the
   threshold filters declared under
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
    merged_vmiss_path: str,
    gentrain_min: float,
    cluster_sep_min: float,
    fmiss_max: float,
    exclude_strand_ambiguous: bool,
    output_audit_tsv_path: str,
    output_exclusion_list_path: str,
    output_summary_tsv_path: str,
    job_name: str = 'snp_qc_report',
) -> list[BashJob]:
    """Queue the jobs.

    Args:
        egt_info_bcf_path: Path to the references-repo EGT info BCF.
        egt_info_bcf_index_path: Path to its ``.csi`` index.
        merged_vmiss_path: Path to the merged-set ``.vmiss`` from ``Plink2Qc``.
        gentrain_min: Inclusive lower bound for ``GenTrain_Score``.
        cluster_sep_min: Inclusive lower bound for ``Cluster_Sep``.
        fmiss_max: Inclusive upper bound for the merged-set ``F_MISS``.
        exclude_strand_ambiguous: When True, also exclude ``{A,T}``/``{C,G}`` SNPs.
        output_audit_tsv_path: Output path for the bgzipped audit TSV.
        output_exclusion_list_path: Output path for the exclusion ``.snplist``.
        output_summary_tsv_path: Output path for the summary TSV.
        job_name: Display name for the synthesis job (the extract job appends
            an ``_extract_egt_info`` suffix).

    Returns:
        Both queued jobs in dependency order: ``[extract_cluster_file_info, generate_snp_exclusion_list]``.
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

    generate_snp_exclusion_list: BashJob = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'snp_qc_report'],
        image=get_driver_image(),
        default_cpu=2,
        default_memory='standard',
        default_storage='20G',
    )
    generate_snp_exclusion_list.depends_on(extract_cluster_file_info)

    script_path: str = str(files('popgen_genotyping.scripts').joinpath('snp_qc_report.py'))
    script_resource = b.read_input(script_path)
    merged_vmiss = b.read_input(merged_vmiss_path)

    generate_snp_exclusion_list.declare_resource_group(
        outputs={
            'audit_tsv_gz': '{root}.audit.tsv.gz',
            'exclusion_list': '{root}.exclude.snplist',
            'summary_tsv': '{root}.summary.tsv',
        },
    )

    strand_flag: str = '--exclude-strand-ambiguous' if exclude_strand_ambiguous else ''

    generate_snp_exclusion_list.command(
        f"""
        set -euxo pipefail
        python3 {script_resource} \\
            --egt-info-tsv {extract_cluster_file_info.egt_info.tsv} \\
            --merged-vmiss {merged_vmiss} \\
            --gentrain-min {gentrain_min} \\
            --cluster-sep-min {cluster_sep_min} \\
            --fmiss-max {fmiss_max} \\
            {strand_flag} \\
            --output-audit-tsv {generate_snp_exclusion_list.outputs.audit_tsv_gz} \\
            --output-exclusion-list {generate_snp_exclusion_list.outputs.exclusion_list} \\
            --output-summary-tsv {generate_snp_exclusion_list.outputs.summary_tsv}
        """,
    )

    b.write_output(generate_snp_exclusion_list.outputs.audit_tsv_gz, output_audit_tsv_path)
    b.write_output(generate_snp_exclusion_list.outputs.exclusion_list, output_exclusion_list_path)
    b.write_output(generate_snp_exclusion_list.outputs.summary_tsv, output_summary_tsv_path)

    return [extract_cluster_file_info, generate_snp_exclusion_list]
