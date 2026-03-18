"""
This file exists to define all the Stages for the workflow.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import TYPE_CHECKING

from cpg_flow.stage import CohortStage, MultiCohortStage, stage
from cpg_utils.config import config_retrieve

from popgen_genotyping.jobs.baf_regress_job import run_bafregress
from popgen_genotyping.jobs.cohort_bcf_to_plink_job import run_cohort_bcf_to_plink
from popgen_genotyping.jobs.export_cohort_datasets_job import run_export_cohort_datasets
from popgen_genotyping.jobs.gtc_to_bcfs_job import run_gtc_to_bcfs
from popgen_genotyping.jobs.merge_cohort_plink_job import run_merge_plink
from popgen_genotyping.jobs.plink2_qc_job import run_plink2_qc
from popgen_genotyping.jobs.plink2_to_plink1_job import run_plink2_to_plink1
from popgen_genotyping.jobs.qc_report_job import run_qc_report
from popgen_genotyping.metamist_utils import (
    query_reported_sex,
    resolve_cohort_gtc_mapping,
    resolve_rolling_aggregate,
)
from popgen_genotyping.utils import get_output_prefix

if TYPE_CHECKING:
    from cpg_flow.stage import StageInput, StageOutput
    from cpg_flow.targets import Cohort, MultiCohort
    from cpg_utils import Path
    from hailtop.batch.job import BashJob


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
            'light_bcf': prefix / f'{cohort.id}.light.bcf',
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


@stage(required_stages=[GtcToBcfs])
class BafRegress(CohortStage):
    """
    Estimate sample contamination for the entire cohort using BAFRegress.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Define the expected BAFRegress output text file for the cohort.
        """
        prefix: Path = get_output_prefix(dataset=cohort.dataset, stage_name=self.name)
        return {
            'bafregress_txt': prefix / f'{cohort.id}.BAFRegress.txt',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Queue the cohort-level BAFRegress estimation job.
        """
        outputs: dict[str, Path] = self.expected_outputs(cohort)

        # Pull the cohort heavy BCF output from the previous stage
        heavy_bcf_path: Path = inputs.as_path(cohort, GtcToBcfs, 'heavy_bcf')

        # Population AF reference optional for BAFRegress
        af_ref_path: str | None = config_retrieve(['popgen_genotyping', 'references', 'af_ref_path'], default=None)

        j: BashJob = run_bafregress(
            bcf_path=str(heavy_bcf_path),
            output_path=str(outputs['bafregress_txt']),
            af_ref_path=af_ref_path,
            job_name=f'BafRegress_{cohort.name}',
        )

        return self.make_outputs(cohort, data=outputs, jobs=[j])


@stage(required_stages=[GtcToBcfs])
class CohortBcfToPlink(CohortStage):
    """
    Convert the cohort-level light BCF to PLINK 1.9 format.
    Output is stored in tmp.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Define the expected PLINK 1.9 binary fileset in temporary storage.
        """
        # Store in tmp per requirement
        prefix: Path = get_output_prefix(dataset=cohort.dataset, stage_name=self.name, tmp=True)
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


@stage(required_stages=[CohortBcfToPlink])
class MergeCohortPlink(MultiCohortStage):
    """
    Merge all cohort PLINK 1.9 datasets into a single unified dataset, with rolling aggregate.
    Output is stored in tmp.
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        """
        Define the expected multi-cohort PLINK 1.9 fileset in temporary storage.
        """
        # Store in tmp per requirement
        prefix: Path = get_output_prefix(dataset=multicohort.analysis_dataset, stage_name=self.name, tmp=True)
        return {
            'bed': prefix / 'merged_cohorts.bed',
            'bim': prefix / 'merged_cohorts.bim',
            'fam': prefix / 'merged_cohorts.fam',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        Queue the multi-cohort PLINK 1.9 merge job, incorporating rolling aggregates if configured.
        """
        outputs: dict[str, Path] = self.expected_outputs(multicohort)

        # 1. Gather cohort outputs
        all_cohort_outputs: dict[str, dict[str, Path]] = inputs.as_dict_by_target(stage=CohortBcfToPlink)
        cohort_plink_paths: list[dict[str, str]] = []
        for _cohort_id, cohort_outs in all_cohort_outputs.items():
            cohort_plink_paths.append(
                {
                    'bed': str(cohort_outs['bed']),
                    'bim': str(cohort_outs['bim']),
                    'fam': str(cohort_outs['fam']),
                }
            )

        # 2. Check for rolling aggregate
        prev_analysis_id: str | None = config_retrieve(
            ['popgen_genotyping', 'merge_cohort_plink', 'previous_analysis_id'], default=None
        )

        previous_aggregate_plink1_paths: dict[str, str] | None = None
        samples_to_remove: list[str] | None = None
        merge_job_dependencies: list[BashJob] = []

        if prev_analysis_id:
            # A previous aggregate exists in PLINK2 format, so we need to convert it to PLINK1.9

            previous_aggregate_plink2_paths, samples_to_remove = resolve_rolling_aggregate(
                prev_analysis_id=prev_analysis_id
            )

            # Define an output prefix for the converted PLINK1.9 files in tmp storage
            plink1_prefix = get_output_prefix(
                dataset=multicohort.analysis_dataset, stage_name='Plink2ToPlink1', tmp=True
            )

            conversion_job, converted_plink1_resource = run_plink2_to_plink1(
                pfile_prefix=previous_aggregate_plink2_paths,
                output_prefix=str(plink1_prefix),
                job_name='Plink2ToPlink1',
            )
            merge_job_dependencies.append(conversion_job)

            # The converted PLINK1.9 files become the previous aggregate for the merge job
            previous_aggregate_plink1_paths = {
                'bed': str(converted_plink1_resource.bed),
                'bim': str(converted_plink1_resource.bim),
                'fam': str(converted_plink1_resource.fam),
            }

        # 3. Call merge job
        j: BashJob = run_merge_plink(
            cohort_plink_paths=cohort_plink_paths,
            output_prefix=str(outputs['bed']).replace('.bed', ''),
            previous_aggregate_paths=previous_aggregate_plink1_paths,
            samples_to_remove=samples_to_remove,
            job_name='MergeCohortPlink',
        )

        if merge_job_dependencies:
            j.depends_on(*merge_job_dependencies)

        return self.make_outputs(multicohort, data=outputs, jobs=[j])


@stage(required_stages=[MergeCohortPlink])
class ExportCohortDatasets(MultiCohortStage):
    """
    Export the merged cohort to PLINK2 format for long-term storage.
    BCF output goes to tmp for analysis, as it is too large for long term storage.
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        """
        Define the expected PLINK2 outputs to long-term storage.
        BCF output goes to tmp for analysis, as it is too large for long term storage.
        """
        prefix: Path = get_output_prefix(dataset=multicohort.analysis_dataset, stage_name=self.name)
        tmp_bcf_prefix: Path = get_output_prefix(dataset=multicohort.analysis_dataset, stage_name=self.name, tmp=True)
        datestamp: str = datetime.now(tz=timezone.utc).strftime('%Y%m%d')
        return {
            'pgen': prefix / f'{datestamp}_cohort.pgen',
            'pvar': prefix / f'{datestamp}_cohort.pvar',
            'psam': prefix / f'{datestamp}_cohort.psam',
            'bcf': tmp_bcf_prefix / f'{datestamp}_cohort.bcf',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        Queue the dataset export job using PLINK2.
        """
        outputs: dict[str, Path] = self.expected_outputs(multicohort=multicohort)

        # 1. Pull input from MergeCohortPlink
        input_plink: dict[str, Path] = inputs.as_dict(target=multicohort, stage=MergeCohortPlink)

        # 2. Call export job
        j: BashJob = run_export_cohort_datasets(
            input_plink_prefix={
                'bed': str(input_plink['bed']),
                'bim': str(input_plink['bim']),
                'fam': str(input_plink['fam']),
            },
            output_prefix=str(outputs['pgen']).replace('.pgen', ''),
            job_name='ExportCohortDatasets',
        )

        return self.make_outputs(multicohort, data=outputs, jobs=[j])


@stage(required_stages=[ExportCohortDatasets])
class Plink2Qc(MultiCohortStage):
    """
    Run PLINK2 QC on a cohort object's pgen/pvar/psam files.
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        """
        Define the expected PLINK2 QC output files for the multi-cohort.
        """
        prefix: Path = get_output_prefix(dataset=multicohort.analysis_dataset, stage_name=self.name)
        output_base_name = 'plink2qc'
        return {
            'smiss': prefix / f'{output_base_name}.smiss',
            'vmiss': prefix / f'{output_base_name}.vmiss',
            'afreq': prefix / f'{output_base_name}.afreq',
            'hwe': prefix / f'{output_base_name}.hwe',
            'het': prefix / f'{output_base_name}.het',
            'sexcheck': prefix / f'{output_base_name}.sexcheck',
            'kin': prefix / f'{output_base_name}.kin',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        Queue the PLINK2 QC job for the multi-cohort.
        """
        outputs: dict[str, Path] = self.expected_outputs(multicohort=multicohort)

        # Get the input PGEN file path from the ExportCohortDatasets stage
        input_plink_pgen: Path = inputs.as_path(target=multicohort, stage=ExportCohortDatasets, key='pgen')

        # The outputs_path for the run_plink2_qc job is the base prefix for all QC files.
        output_plink2_prefix = str(outputs['smiss']).removesuffix('.smiss')

        # Call the Hail Batch job function
        j: BashJob = run_plink2_qc(
            pgen_path=str(input_plink_pgen),
            outputs_path=output_plink2_prefix,
            job_name=f'Plink2Qc_{multicohort.name}',
        )

        # Return the expected outputs of this stage, referencing the outputs generated by the job
        return self.make_outputs(multicohort, data=outputs, jobs=[j])


@stage(required_stages=[Plink2Qc, BafRegress])
class QcReport(MultiCohortStage):
    """
    Create the QC report for an input object.
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        """
        Define the expected QC report output file for the multi-cohort.
        """
        prefix: Path = get_output_prefix(dataset=multicohort.analysis_dataset, stage_name=self.name)
        datestamp: str = datetime.now(tz=timezone.utc).strftime('%Y%m%d')
        return {
            'qc_report_csv': prefix / f'{datestamp}_qc_report.csv',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        Queue the QC report generation job for the multi-cohort.
        """
        outputs: dict[str, Path] = self.expected_outputs(multicohort=multicohort)

        # Get the plink2_qc_prefix from the Plink2Qc stage's 'smiss' output
        plink_qc_smiss_path: Path = inputs.as_path(target=multicohort, stage=Plink2Qc, key='smiss')
        plink_qc_prefix = str(plink_qc_smiss_path).removesuffix('.smiss')

        # Get all bafregress output paths from all cohorts
        bafregress_outputs: dict[str, dict[str, Path]] = inputs.as_dict_by_target(stage=BafRegress)
        bafregress_paths: list[str] = [str(baf_out['bafregress_txt']) for baf_out in bafregress_outputs.values()]

        # Call the Hail Batch job function
        j: BashJob = run_qc_report(
            plink_qc_prefix=plink_qc_prefix,
            bafregress_paths=bafregress_paths,
            output_path=str(outputs['qc_report_csv']),
            job_name=f'QcReport_{multicohort.name}',
        )

        # Return the expected outputs of this stage
        return self.make_outputs(multicohort, data=outputs, jobs=[j])
