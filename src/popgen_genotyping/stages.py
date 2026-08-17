"""
This file exists to define all the Stages for the workflow.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from cpg_flow.stage import CohortStage, stage
from cpg_utils.config import config_retrieve, reference_path
from loguru import logger

from popgen_genotyping.jobs.baf_regress_job import run_bafregress
from popgen_genotyping.jobs.cohort_bcf_to_plink_job import run_cohort_bcf_to_plink
from popgen_genotyping.jobs.export_cohort_datasets_job import run_export_cohort_datasets
from popgen_genotyping.jobs.gtc_to_bcfs_job import run_gtc_to_bcfs
from popgen_genotyping.jobs.king_ibdseg_job import run_king_ibdseg
from popgen_genotyping.jobs.merge_cohort_plink_job import run_merge_plink
from popgen_genotyping.jobs.plink2_qc_job import run_plink2_qc
from popgen_genotyping.jobs.plink2_to_plink1_job import run_plink2_to_plink1
from popgen_genotyping.jobs.qc_report_job import run_qc_report
from popgen_genotyping.jobs.snp_qc_report_job import run_snp_qc_report
from popgen_genotyping.metamist_utils import (
    format_merge_plan,
    query_reported_sex,
    resolve_bafregress_map,
    resolve_cohort_gtc_mapping,
    resolve_merge_inputs,
)
from popgen_genotyping.utils import get_output_prefix

if TYPE_CHECKING:
    from cpg_flow.stage import StageInput, StageOutput
    from cpg_flow.targets import Cohort
    from cpg_utils import Path
    from hailtop.batch.job import BashJob
    from hailtop.batch.resource import ResourceGroup


@stage
class GtcToBcfs(CohortStage):
    """
    Convert all cohort GTC files to Heavy and Light multi-sample BCFs.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Define the expected cohort-level multi-sample BCF and metadata outputs.
        """
        prefix: Path = get_output_prefix(dataset=cohort.dataset, stage_name=self.name, tmp=True)
        return {
            'heavy_bcf': prefix / f'{cohort.id}.heavy.bcf',
            'heavy_bcf_index': prefix / f'{cohort.id}.heavy.bcf.csi',
            'light_bcf': prefix / f'{cohort.id}.light.bcf',
            'light_bcf_index': prefix / f'{cohort.id}.light.bcf.csi',
            'metadata_tsv': prefix / f'{cohort.id}_gtc_metadata.tsv',
        }

    def queue_jobs(self, cohort: Cohort, _inputs: StageInput) -> StageOutput:
        """
        Queue the cohort-level multi-sample GTC to BCF conversion job.
        """
        outputs: dict[str, Path] = self.expected_outputs(cohort)

        # Retrieve reference paths from config
        fasta_ref: str = config_retrieve(['popgen_genotyping', 'references', 'fasta_ref_path'])
        bpm_manifest: str = config_retrieve(['popgen_genotyping', 'references', 'bpm_manifest_path'])
        egt_cluster: str = config_retrieve(['popgen_genotyping', 'references', 'egt_cluster_path'])

        # Resolve GTC mapping (SG_ID -> {gtc, old_name})
        mapping_data: dict[str, dict[str, str]] = resolve_cohort_gtc_mapping(cohort=cohort)

        gtc_paths: list[str] = [d['gtc'] for d in mapping_data.values()]
        # sample_mapping: barcode_pos -> SG_ID
        sample_mapping: dict[str, str] = {d['old_name']: sg_id for sg_id, d in mapping_data.items()}

        j: BashJob = run_gtc_to_bcfs(
            gtc_paths=gtc_paths,
            sample_mapping=sample_mapping,
            output_heavy_bcf_path=str(outputs['heavy_bcf']),
            output_light_bcf_path=str(outputs['light_bcf']),
            output_metadata_path=str(outputs['metadata_tsv']),
            bpm_manifest_path=bpm_manifest,
            egt_cluster_path=egt_cluster,
            fasta_ref_path=fasta_ref,
            job_name=f'GtcToBcfs_{cohort.name}',
        )

        return self.make_outputs(cohort, data=outputs, jobs=[j])


@stage(required_stages=[GtcToBcfs], analysis_type='array_bafregress')
class BafRegress(CohortStage):
    """
    Estimate sample contamination for the entire cohort using BAFRegress.

    Output is written to durable, version-independent storage and registered against the
    cohort as an ``array_bafregress`` analysis — computed once per plate and reused across
    runs. The phase-2 QC report resolves these by querying all ``array_bafregress`` analyses
    (``resolve_bafregress_map``), so the final table covers every constituent cohort of the
    aggregate.
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        """
        Define the expected BAFRegress output text file in durable, unversioned storage.
        """
        # Durable, unversioned per-cohort path: processed once, reused across versions.
        prefix: Path = get_output_prefix(dataset=cohort.dataset, stage_name=self.name, versioned=False)
        return prefix / f'{cohort.id}.BAFRegress.txt'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Queue the cohort-level BAFRegress estimation job.
        """
        outputs: Path = self.expected_outputs(cohort)

        # Pull the cohort heavy BCF output from the previous stage
        heavy_bcf_path: Path = inputs.as_path(cohort, GtcToBcfs, 'heavy_bcf')

        # Population AF reference optional for BAFRegress
        af_ref_path: str | None = config_retrieve(['popgen_genotyping', 'references', 'af_ref_path'], default=None)

        j: BashJob = run_bafregress(
            bcf_path=str(heavy_bcf_path),
            output_path=str(outputs),
            af_ref_path=af_ref_path,
            job_name=f'BafRegress_{cohort.name}',
        )

        return self.make_outputs(cohort, data=outputs, jobs=[j])


@stage(required_stages=[GtcToBcfs], analysis_type='array_cohort_bed', analysis_keys=['bed'])
class CohortBcfToPlink(CohortStage):
    """
    Convert the cohort-level light BCF to PLINK 1.9 format.

    Output is written to durable, version-independent storage and registered against the
    cohort as an ``array_cohort_bed`` analysis. Phase 2 discovers these filesets by querying the
    ``array_cohort_bed`` analyses in Metamist (``resolve_merge_inputs``), not via stage-wiring;
    the fileset is processed once per plate and reused across runs.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Define the expected PLINK 1.9 binary fileset in durable, unversioned storage.
        """
        # Durable, unversioned per-cohort path: processed once, reused across versions.
        prefix: Path = get_output_prefix(dataset=cohort.dataset, stage_name=self.name, versioned=False)
        return {
            'bed': prefix / f'{cohort.id}.bed',
            'bim': prefix / f'{cohort.id}.bim',
            'fam': prefix / f'{cohort.id}.fam',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Queue the cohort-level BCF to PLINK 1.9 conversion job.
        """
        outputs: dict[str, Path] = self.expected_outputs(cohort=cohort)

        # Pull the cohort light BCF path from GtcToBcfs
        light_bcf_path: Path = inputs.as_path(target=cohort, stage=GtcToBcfs, key='light_bcf')

        # Fetch reported sex metadata for these SGs
        full_sex_mapping: dict[str, str] = query_reported_sex(project=cohort.dataset.name)
        cohort_sg_ids: set[str] = set(cohort.get_sequencing_group_ids())
        sex_mapping: dict[str, str] = {
            sg_id: sex_code for sg_id, sex_code in full_sex_mapping.items() if sg_id in cohort_sg_ids
        }

        # Define the job via the job utility
        j: BashJob = run_cohort_bcf_to_plink(
            bcf_path=str(light_bcf_path),
            output_prefix=str(outputs['bed']).replace('.bed', ''),
            sex_mapping=sex_mapping,
            job_name=f'CohortBcfToPlink_{cohort.name}',
        )

        return self.make_outputs(cohort, data=outputs, jobs=[j])


@stage
class MergeCohortPlink(CohortStage):
    """
    Merge the phase-1 per-plate PLINK 1.9 datasets into the super cohort, with rolling aggregate.

    Runs in phase 2 against the manually-created super cohort. Inputs are resolved from Metamist
    (``resolve_merge_inputs``), not stage-wiring — hence no ``required_stages``: the super cohort's
    membership is the source of truth, the previous aggregate is carried forward from the configured
    ``previous_aggregate_cohort_id``, and the contributing plates are derived as
    ``NEW = super - previous aggregate`` and mapped to their registered ``array_cohort_bed`` filesets.
    The resolved plan is logged for operator confirmation. Whole plate filesets are merged and a final
    ``--keep`` trims the result to super-cohort membership. Output is stored in tmp.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Define the expected merged PLINK 1.9 fileset in temporary storage.
        """
        # Store in tmp per requirement. Keyed by the super-cohort ID: every rolling aggregate
        # gets a new super cohort, so two runs at the same workflow.version land on distinct
        # paths — otherwise cpg-flow's skip-if-exists would silently reuse the previous super
        # cohort's merge (and the in-job --keep count assert would never fire).
        prefix: Path = get_output_prefix(dataset=cohort.dataset, stage_name=self.name, tmp=True)
        return {
            'bed': prefix / f'{cohort.id}_merged.bed',
            'bim': prefix / f'{cohort.id}_merged.bim',
            'fam': prefix / f'{cohort.id}_merged.fam',
        }

    def queue_jobs(self, cohort: Cohort, _inputs: StageInput) -> StageOutput:
        """
        Queue the PLINK 1.9 merge job for the super cohort, resolving inputs from Metamist.
        """
        outputs: dict[str, Path] = self.expected_outputs(cohort)

        # 1. Resolve the merge plan from Metamist. The super cohort is the source of truth;
        #    new plates are derived (NEW = super - previous aggregate), not listed in config.
        previous_aggregate_cohort_id: str | None = config_retrieve(
            ['popgen_genotyping', 'merge_cohort_plink', 'previous_aggregate_cohort_id'], default=None
        )
        super_cohort_sg_ids: list[str] = cohort.get_sequencing_group_ids()
        resolved: dict = resolve_merge_inputs(
            super_cohort_sg_ids=super_cohort_sg_ids,
            previous_aggregate_cohort_id=previous_aggregate_cohort_id,
        )

        # Log the plan so an operator can confirm the derived plates match the phase-1 runs.
        # loguru, not stdlib logging: cpg-flow configures loguru, while unconfigured stdlib
        # logging drops INFO — the plan would never reach the driver log.
        logger.info(format_merge_plan(resolved, previous_aggregate_cohort_id))

        cohort_plink_paths: list[dict[str, str]] = [
            {'bed': plate['bed'], 'bim': plate['bim'], 'fam': plate['fam']} for plate in resolved['plate_merge_list']
        ]

        # 2. Carry the previous aggregate forward, if any. It is stored as PLINK2, so convert
        #    it back to PLINK 1.9 before the merge.
        previous_aggregate_plink1_resource: ResourceGroup | None = None
        merge_job_dependencies: list[BashJob] = []

        if resolved['previous_aggregate_paths']:
            plink1_prefix = get_output_prefix(dataset=cohort.dataset, stage_name='Plink2ToPlink1', tmp=True)

            conversion_job, converted_plink1_resource = run_plink2_to_plink1(
                pfile_prefix=resolved['previous_aggregate_paths'],
                output_prefix=str(plink1_prefix / f'{cohort.id}_plink1'),
                job_name='Plink2ToPlink1',
            )
            merge_job_dependencies.append(conversion_job)
            previous_aggregate_plink1_resource = converted_plink1_resource

        # 3. Call merge job. keep_samples trims the whole-plate merge to super-cohort membership.
        j: BashJob = run_merge_plink(
            cohort_plink_paths=cohort_plink_paths,
            output_prefix=str(outputs['bed']).replace('.bed', ''),
            keep_samples=super_cohort_sg_ids,
            previous_aggregate_resource=previous_aggregate_plink1_resource,
            samples_to_remove=resolved['samples_to_remove'],
            job_name='MergeCohortPlink',
        )

        if merge_job_dependencies:
            j.depends_on(*merge_job_dependencies)

        return self.make_outputs(cohort, data=outputs, jobs=[j])


@stage(required_stages=[MergeCohortPlink], analysis_type='array_aggregate_pgen', analysis_keys=['pgen'])
class ExportCohortDatasets(CohortStage):
    """
    Export the merged cohort to PLINK2 format for long-term storage.
    BCF output goes to tmp for analysis, as it is too large for long term storage.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Define the expected PLINK2 outputs to long-term storage.
        BCF output goes to tmp for analysis, as it is too large for long term storage.
        """
        # No datestamp: the cohort ID (plus the versioned prefix) already distinguishes
        # aggregates, and paths must be stable across days for skip-if-exists and
        # single-stage reruns to find upstream outputs. Same for all phase-2 stages.
        prefix: Path = get_output_prefix(dataset=cohort.dataset, stage_name=self.name)
        tmp_bcf_prefix: Path = get_output_prefix(dataset=cohort.dataset, stage_name=self.name, tmp=True)
        return {
            'pgen': prefix / f'{cohort.id}.pgen',
            'pvar': prefix / f'{cohort.id}.pvar',
            'psam': prefix / f'{cohort.id}.psam',
            'bcf': tmp_bcf_prefix / f'{cohort.id}.bcf',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Queue the dataset export job using PLINK2.
        """
        outputs: dict[str, Path] = self.expected_outputs(cohort=cohort)

        # 1. Pull input from MergeCohortPlink
        input_plink: dict[str, Path] = inputs.as_dict(target=cohort, stage=MergeCohortPlink)

        # 2. Call export job
        j: BashJob = run_export_cohort_datasets(
            input_plink_prefix={
                'bed': str(input_plink['bed']),
                'bim': str(input_plink['bim']),
                'fam': str(input_plink['fam']),
            },
            output_prefix=str(outputs['pgen']).replace('.pgen', ''),
            bcf_output_path=str(outputs['bcf']),
            job_name='ExportCohortDatasets',
        )

        return self.make_outputs(cohort, data=outputs, jobs=[j])


@stage(required_stages=[ExportCohortDatasets], analysis_type='array_qc_raw', analysis_keys=['log'])
class Plink2Qc(CohortStage):
    """
    Per-sample PLINK2 QC on the merged pgen/pvar/psam: missingness, inbreeding,
    sex-check, plus per-variant allele frequencies. Per-variant missingness and
    Hardy-Weinberg are computed inside ``SnpQcReport``.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Define the expected PLINK2 QC output files for the cohort.
        """
        prefix: Path = get_output_prefix(dataset=cohort.dataset, stage_name=self.name)
        output_base_name = f'{cohort.id}_qc'
        return {
            'smiss': prefix / f'{output_base_name}.smiss',
            'afreq': prefix / f'{output_base_name}.afreq',
            'het': prefix / f'{output_base_name}.het',
            'sexcheck': prefix / f'{output_base_name}.sexcheck',
            'log': prefix / f'{output_base_name}.log',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Queue the PLINK2 QC job for the cohort.
        """
        outputs: dict[str, Path] = self.expected_outputs(cohort=cohort)

        # Get the input PGEN file path from the ExportCohortDatasets stage
        input_plink_pgen: Path = inputs.as_path(target=cohort, stage=ExportCohortDatasets, key='pgen')

        # The outputs_path for the run_plink2_qc job is the base prefix for all QC files.
        output_plink2_prefix = str(outputs['smiss']).removesuffix('.smiss')

        # Call the Hail Batch job function
        j: BashJob = run_plink2_qc(
            pgen_path=str(input_plink_pgen),
            outputs_path=output_plink2_prefix,
            job_name=f'Plink2Qc_{cohort.name}',
        )

        # Return the expected outputs of this stage, referencing the outputs generated by the job
        return self.make_outputs(cohort, data=outputs, jobs=[j])


@stage(
    required_stages=[MergeCohortPlink],
    analysis_type='array_relatedness_ibdseg',
    analysis_keys=['seg', 'seg_x'],
)
class KingIbdseg(CohortStage):
    """
    Infer pairwise IBD segments across the merged cohort with KING `--ibdseg`.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Define the expected KING `--ibdseg` outputs for the cohort.

        KING 2.3.2 produces an X-chromosome companion (`{prefix}X.seg` /
        `{prefix}X.segments.gz`, no dot before the `X`) whenever the merged
        PLINK fileset has chrX SNPs. The job backfills header-only placeholders
        when the input has no chrX, so these outputs are always materialised.
        """
        prefix: Path = get_output_prefix(dataset=cohort.dataset, stage_name=self.name)
        output_base_name = f'{cohort.id}_king'
        return {
            'seg': prefix / f'{output_base_name}.seg',
            'segments': prefix / f'{output_base_name}.segments.gz',
            'seg_x': prefix / f'{output_base_name}X.seg',
            'segments_x': prefix / f'{output_base_name}X.segments.gz',
            'log': prefix / f'{output_base_name}.log',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Queue the KING `--ibdseg` job against the merged PLINK 1.9 dataset.
        """
        outputs: dict[str, Path] = self.expected_outputs(cohort=cohort)

        merged_plink: dict[str, Path] = inputs.as_dict(target=cohort, stage=MergeCohortPlink)

        jobs: list[BashJob] = run_king_ibdseg(
            bed_path=str(merged_plink['bed']),
            bim_path=str(merged_plink['bim']),
            fam_path=str(merged_plink['fam']),
            output_seg_path=str(outputs['seg']),
            output_segments_path=str(outputs['segments']),
            output_seg_x_path=str(outputs['seg_x']),
            output_segments_x_path=str(outputs['segments_x']),
            output_log_path=str(outputs['log']),
            job_name=f'KingIbdseg_{cohort.name}',
        )

        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage(
    required_stages=[ExportCohortDatasets],
    analysis_type='array_snp_qc',
    analysis_keys=['inclusion_list'],
)
class SnpQcReport(CohortStage):
    """
    Per-SNP QC: vendor cluster scores + call rate + Hardy-Weinberg.

    Runs ``plink2 --missing`` and ``plink2 --hwe`` against the merged PGEN to
    derive merged-set F_MISS and a Hardy-Weinberg pass list, joins those with
    the references-repo EGT INFO BCF (``GenTrain_Score`` / ``Cluster_Sep``),
    and applies the thresholds declared under
    ``[popgen_genotyping.snp_qc_report.thresholds]``. Emits an audit TSV, an
    inclusion ``.snplist`` (passing variant IDs, for ``plink2 --extract``), and
    a per-filter summary TSV. The PGEN released by ``ExportCohortDatasets`` is
    not modified.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Define the audit TSV, inclusion list, and summary TSV outputs.
        """
        prefix: Path = get_output_prefix(dataset=cohort.dataset, stage_name=self.name)
        return {
            'audit_tsv': prefix / f'{cohort.id}_snp_qc.audit.tsv.gz',
            'inclusion_list': prefix / f'{cohort.id}_snp_qc.include.snplist',
            'summary_tsv': prefix / f'{cohort.id}_snp_qc.summary.tsv',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Queue the EGT-INFO extract, merged-set variant-metrics, and filter jobs.
        """
        outputs: dict[str, Path] = self.expected_outputs(cohort=cohort)
        merged_pgen_path: Path = inputs.as_path(target=cohort, stage=ExportCohortDatasets, key='pgen')
        merged_pvar_path: Path = inputs.as_path(target=cohort, stage=ExportCohortDatasets, key='pvar')
        merged_psam_path: Path = inputs.as_path(target=cohort, stage=ExportCohortDatasets, key='psam')

        stage_path: list[str] = ['popgen_genotyping', 'snp_qc_report']
        thresholds_path: list[str] = [*stage_path, 'thresholds']
        gentrain_min: float = config_retrieve([*thresholds_path, 'gentrain_min'], 0.7)
        cluster_sep_min: float = config_retrieve([*thresholds_path, 'cluster_sep_min'], 0.4)
        fmiss_max: float = config_retrieve([*thresholds_path, 'fmiss_max'], 0.02)
        hwe_p: float = config_retrieve([*thresholds_path, 'hwe_p'], 1e-5)
        hwe_k: float = config_retrieve([*thresholds_path, 'hwe_k'], 0.001)
        hwe_midp: bool = config_retrieve([*thresholds_path, 'hwe_midp'], True)
        hwe_keep_fewhet: bool = config_retrieve([*thresholds_path, 'hwe_keep_fewhet'], True)

        egt_info_bcf_key: str = config_retrieve(
            [*stage_path, 'egt_info_bcf_reference_key'],
            'illumina_microarray/GDA_8v1_0_D1_ClusterFile_egt_info_bcf',
        )
        egt_info_bcf_index_key: str = config_retrieve(
            [*stage_path, 'egt_info_bcf_index_reference_key'],
            'illumina_microarray/GDA_8v1_0_D1_ClusterFile_egt_info_bcf_index',
        )

        jobs: list[BashJob] = run_snp_qc_report(
            egt_info_bcf_path=reference_path(egt_info_bcf_key),
            egt_info_bcf_index_path=reference_path(egt_info_bcf_index_key),
            merged_pgen_path=str(merged_pgen_path),
            merged_pvar_path=str(merged_pvar_path),
            merged_psam_path=str(merged_psam_path),
            gentrain_min=gentrain_min,
            cluster_sep_min=cluster_sep_min,
            fmiss_max=fmiss_max,
            hwe_p=hwe_p,
            hwe_k=hwe_k,
            hwe_midp=hwe_midp,
            hwe_keep_fewhet=hwe_keep_fewhet,
            output_audit_tsv_path=str(outputs['audit_tsv']),
            output_inclusion_list_path=str(outputs['inclusion_list']),
            output_summary_tsv_path=str(outputs['summary_tsv']),
            job_name=f'SnpQcReport_{cohort.name}',
        )

        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage(required_stages=[Plink2Qc, KingIbdseg], analysis_type='array_qc_report')
class QcReport(CohortStage):
    """
    Create the QC report for the super cohort.

    ``BafRegress`` is a phase-1 (per-plate) stage that does not run in phase 2, so it is not a
    ``required_stages`` dependency. Its per-plate contamination outputs are resolved by querying
    all registered ``array_bafregress`` analyses in Metamist (``resolve_bafregress_map``) and
    selecting the super cohort's full membership — covering every constituent plate, not just the
    new ones.
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        """
        Define the expected QC report output file for the cohort.
        """
        prefix: Path = get_output_prefix(dataset=cohort.dataset, stage_name=self.name)
        return prefix / f'{cohort.id}_qc_report.csv'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Queue the QC report generation job for the cohort.
        """
        outputs: Path = self.expected_outputs(cohort=cohort)

        # Get the plink2_qc_prefix from the Plink2Qc stage's 'smiss' output
        plink_qc_smiss_path: Path = inputs.as_path(target=cohort, stage=Plink2Qc, key='smiss')
        plink_qc_prefix = str(plink_qc_smiss_path).removesuffix('.smiss')

        # Autosomal KING --ibdseg pairwise summary (REL_ID:KINSHIP:INFTYPE)
        king_seg_path: Path = inputs.as_path(target=cohort, stage=KingIbdseg, key='seg')

        # Resolve the BafRegress output for every super-cohort SG from Metamist (full membership,
        # not just the new plates), since BafRegress does not run as a phase-2 stage.
        # The map is sg_id -> plate-level file (one file per plate cohort), so dedupe before
        # passing it on: repeating a plate file once per SG would multiply that plate's rows
        # in the report's IID merge. Sorted for a reproducible job command.
        bafregress_map: dict[str, str] = resolve_bafregress_map(sg_ids=cohort.get_sequencing_group_ids())
        bafregress_paths: list[str] = sorted(set(bafregress_map.values()))

        # Call the Hail Batch job function
        j: BashJob = run_qc_report(
            plink_qc_prefix=plink_qc_prefix,
            king_seg_path=str(king_seg_path),
            bafregress_paths=bafregress_paths,
            output_path=str(outputs),
            job_name=f'QcReport_{cohort.name}',
        )

        # Return the expected outputs of this stage
        return self.make_outputs(cohort, data=outputs, jobs=[j])
