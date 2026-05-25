"""
Two-job chain that materialises the multi-cohort SNP QC report.

Job 1 (`bcftools_image`) pulls the per-variant INFO payload (`GenTrain_Score`,
`Cluster_Sep`, training cluster counts) out of the references-repo EGT info
BCF and emits a header-less TSV. Job 2 (`driver_image`, depends on Job 1)
runs the pandas synthesis script, joining EGT INFO with merged Plink2Qc,
per-cohort CohortPlinkQc, and per-cohort CohortClusterStats. Outputs:
audit TSV (.gz), exclusion .snplist, summary TSV.
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
    merged_afreq_path: str,
    merged_hardy_path: str,
    cohort_plink_qc: dict[str, dict[str, str]],
    cohort_cluster_stats: dict[str, str],
    thresholds: dict[str, float],
    output_audit_tsv_path: str,
    output_exclusion_list_path: str,
    output_summary_tsv_path: str,
    job_name: str = 'snp_qc_report',
) -> list[BashJob]:
    """
    Queue the EGT-INFO-extract + synthesis job pair.

    Args:
        egt_info_bcf_path: Cloud path to the EGT info BCF (from references repo).
        egt_info_bcf_index_path: Cloud path to its `.csi` index.
        merged_vmiss_path: Cloud path to merged-set `.vmiss` (from `Plink2Qc`).
        merged_afreq_path: Cloud path to merged-set `.afreq`.
        merged_hardy_path: Cloud path to merged-set `.hardy`.
        cohort_plink_qc: ``{cohort_id: {'vmiss': str, 'afreq': str, 'hardy': str}}``
            from `CohortPlinkQc` for every cohort in the current multi-cohort.
        cohort_cluster_stats: ``{cohort_id: cluster_stats.tsv.gz path}`` from
            `CohortClusterStats`. May be a subset of `cohort_plink_qc` keys.
        thresholds: Dict with keys ``gentrain_min``, ``cluster_sep_min``,
            ``vmiss_max``, ``hwe_p_min``, ``call_rate_range_max``,
            ``af_chi2_p_min``. Passed verbatim to the synthesis CLI.
        output_audit_tsv_path: Cloud path for the bgzipped per-variant audit TSV.
        output_exclusion_list_path: Cloud path for the exclusion `.snplist`.
        output_summary_tsv_path: Cloud path for the per-filter summary TSV.
        job_name: Display name for the synthesis job (the extract job appends
            a suffix).

    Returns:
        Both queued jobs, in dependency order: ``[extract_job, synth_job]``.
    """
    b: Batch = get_batch()

    # --- Job 1: bcftools query → EGT INFO TSV ---
    extract_job: BashJob = register_job(
        batch=b,
        job_name=f'{job_name}_extract_egt_info',
        config_path=['popgen_genotyping', 'snp_qc_report', 'extract'],
        image=config_retrieve(['workflow', 'bcftools_image']),
        default_cpu=1,
        default_storage='10G',
    )

    egt_input = b.read_input_group(bcf=egt_info_bcf_path, csi=egt_info_bcf_index_path)
    extract_job.declare_resource_group(egt_info={'tsv': '{root}.egt_info.tsv'})

    egt_info_fields = (
        r'%CHROM\t%POS\t%REF\t%ALT\t%ID\t'
        r'%INFO/GenTrain_Score\t%INFO/Cluster_Sep\t'
        r'%INFO/N_AA\t%INFO/N_AB\t%INFO/N_BB\t'
        r'%INFO/meanTHETA_AA\t%INFO/meanTHETA_AB\t%INFO/meanTHETA_BB\t'
        r'%INFO/devTHETA_AA\t%INFO/devTHETA_AB\t%INFO/devTHETA_BB\n'
    )
    extract_job.command(
        f"""
        set -euxo pipefail
        bcftools query \\
            -f '{egt_info_fields}' \\
            {egt_input.bcf} > {extract_job.egt_info.tsv}
        """,
    )

    # --- Job 2: pandas synthesis ---
    synth_job: BashJob = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'snp_qc_report'],
        image=get_driver_image(),
        default_cpu=4,
        default_memory='highmem',
        default_storage='50G',
    )
    synth_job.depends_on(extract_job)

    script_path = str(files('popgen_genotyping.scripts').joinpath('snp_qc_report.py'))
    script_resource = b.read_input(script_path)

    merged_vmiss = b.read_input(merged_vmiss_path)
    merged_afreq = b.read_input(merged_afreq_path)
    merged_hardy = b.read_input(merged_hardy_path)

    cohort_vmiss_args: list[str] = []
    cohort_afreq_args: list[str] = []
    cohort_hardy_args: list[str] = []
    for cohort_id, paths in cohort_plink_qc.items():
        vmiss_res = b.read_input(paths['vmiss'])
        afreq_res = b.read_input(paths['afreq'])
        hardy_res = b.read_input(paths['hardy'])
        cohort_vmiss_args.append(f'--cohort-vmiss {cohort_id}:{vmiss_res}')
        cohort_afreq_args.append(f'--cohort-afreq {cohort_id}:{afreq_res}')
        cohort_hardy_args.append(f'--cohort-hardy {cohort_id}:{hardy_res}')

    cohort_stats_args: list[str] = []
    for cohort_id, stats_path in cohort_cluster_stats.items():
        stats_res = b.read_input(stats_path)
        cohort_stats_args.append(f'--cohort-cluster-stats {cohort_id}:{stats_res}')

    synth_job.declare_resource_group(
        outputs={
            'audit_tsv_gz': '{root}.audit.tsv.gz',
            'exclusion_list': '{root}.exclude.snplist',
            'summary_tsv': '{root}.summary.tsv',
        },
    )

    cohort_args_joined = ' '.join([*cohort_vmiss_args, *cohort_afreq_args, *cohort_hardy_args, *cohort_stats_args])

    synth_job.command(
        f"""
        set -euxo pipefail
        python3 {script_resource} \\
            --egt-info-tsv {extract_job.egt_info.tsv} \\
            --merged-vmiss {merged_vmiss} \\
            --merged-afreq {merged_afreq} \\
            --merged-hardy {merged_hardy} \\
            {cohort_args_joined} \\
            --gentrain-min {thresholds['gentrain_min']} \\
            --cluster-sep-min {thresholds['cluster_sep_min']} \\
            --vmiss-max {thresholds['vmiss_max']} \\
            --hwe-p-min {thresholds['hwe_p_min']} \\
            --call-rate-range-max {thresholds['call_rate_range_max']} \\
            --af-chi2-p-min {thresholds['af_chi2_p_min']} \\
            --egt-centroid-delta-max {thresholds['egt_centroid_delta_max']} \\
            --cluster-icc-max {thresholds['cluster_icc_max']} \\
            --cluster-spread-max {thresholds['cluster_spread_max']} \\
            --cluster-separation-min {thresholds['cluster_separation_min']} \\
            --output-audit-tsv {synth_job.outputs.audit_tsv_gz} \\
            --output-exclusion-list {synth_job.outputs.exclusion_list} \\
            --output-summary-tsv {synth_job.outputs.summary_tsv}
        """,
    )

    b.write_output(synth_job.outputs.audit_tsv_gz, output_audit_tsv_path)
    b.write_output(synth_job.outputs.exclusion_list, output_exclusion_list_path)
    b.write_output(synth_job.outputs.summary_tsv, output_summary_tsv_path)

    return [extract_job, synth_job]
