"""
This file exists to define all the Stages for the workflow.
"""

from typing import TYPE_CHECKING

from cpg_flow.stage import CohortStage, MultiCohortStage, stage
from cpg_utils.config import config_retrieve

from popgen_genotyping.jobs.bafregress_job import run_bafregress
from popgen_genotyping.jobs.cohort_bcf_to_plink_job import run_cohort_bcf_to_plink
from popgen_genotyping.jobs.export_cohort_datasets_job import run_export_cohort_datasets
from popgen_genotyping.jobs.gtc_to_bcfs_job import run_gtc_to_bcfs
from popgen_genotyping.jobs.merge_plink_job import run_merge_plink
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


@stage
class GtcToBcfs(CohortStage):
    """
    Convert all cohort GTC files to Heavy and Light multi-sample BCFs.
    """

    def expected_outputs(self, cohort: 'Cohort') -> dict[str, 'Path']:
        """
        Define the expected cohort-level multi-sample BCF and metadata outputs.
        """
        prefix = get_output_prefix(cohort.analysis_dataset, self.name)
        return {
            'heavy_bcf': prefix / 'cohort.heavy.bcf',
            'light_bcf': prefix / 'cohort.light.bcf',
            'metadata_tsv': prefix / 'cohort_gtc_metadata.tsv',
        }

    def queue_jobs(self, cohort: 'Cohort', _inputs: 'StageInput') -> 'StageOutput':
        """
        Queue the cohort-level multi-sample GTC to BCF conversion job.
        """
        outputs = self.expected_outputs(cohort)

        # Retrieve reference paths from config
        fasta_ref = config_retrieve(['popgen_genotyping', 'references', 'fasta_ref_path'])
        bpm_manifest = config_retrieve(['popgen_genotyping', 'references', 'bpm_manifest_path'])
        egt_cluster = config_retrieve(['popgen_genotyping', 'references', 'egt_cluster_path'])

        # Resolve GTC mapping (SG_ID -> {gtc, old_name})
        mapping_data = resolve_cohort_gtc_mapping(cohort)
        
        gtc_paths = [d['gtc'] for d in mapping_data.values()]
        # sample_mapping: barcode_pos -> SG_ID
        sample_mapping = {d['old_name']: sg_id for sg_id, d in mapping_data.items()}

        j = run_gtc_to_bcfs(
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

    def expected_outputs(self, cohort: 'Cohort') -> dict[str, 'Path']:
        """
        Define the expected BAFRegress output text file for the cohort.
        """
        prefix = get_output_prefix(cohort.analysis_dataset, self.name)
        return {
            'bafregress_txt': prefix / 'cohort.BAFRegress.txt',
        }

    def queue_jobs(self, cohort: 'Cohort', inputs: 'StageInput') -> 'StageOutput':
        """
        Queue the cohort-level BAFRegress estimation job.
        """
        outputs = self.expected_outputs(cohort)

        # Pull the cohort heavy BCF output from the previous stage
        heavy_bcf_path = inputs.as_path(cohort, GtcToBcfs, 'heavy_bcf')
        
        # Population AF reference mandatory for BAFRegress
        af_ref_path = config_retrieve(['popgen_genotyping', 'references', 'af_ref_path'])

        j = run_bafregress(
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

    def expected_outputs(self, cohort: 'Cohort') -> dict[str, 'Path']:
        """
        Define the expected PLINK 1.9 binary fileset in temporary storage.
        """
        # Store in tmp per requirement
        prefix = get_output_prefix(cohort.analysis_dataset, self.name, tmp=True)
        return {
            'bed': prefix / 'cohort.bed',
            'bim': prefix / 'cohort.bim',
            'fam': prefix / 'cohort.fam',
        }

    def queue_jobs(self, cohort: 'Cohort', inputs: 'StageInput') -> 'StageOutput':
        """
        Queue the cohort-level BCF to PLINK 1.9 conversion job.
        """
        outputs = self.expected_outputs(cohort)

        # Pull the cohort light BCF path from GtcToBcfs
        light_bcf_path = inputs.as_path(cohort, GtcToBcfs, 'light_bcf')

        # Fetch reported sex metadata for these SGs
        full_sex_mapping = query_reported_sex(cohort.analysis_dataset.name)
        cohort_sg_ids = set(cohort.get_sequencing_group_ids())
        sex_mapping = {sg_id: sex_code for sg_id, sex_code in full_sex_mapping.items() if sg_id in cohort_sg_ids}

        # Define the job via the job utility
        j = run_cohort_bcf_to_plink(
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

    def expected_outputs(self, multicohort: 'MultiCohort') -> dict[str, 'Path']:
        """
        Define the expected multi-cohort PLINK 1.9 fileset in temporary storage.
        """
        # Store in tmp per requirement
        prefix = get_output_prefix(multicohort.analysis_dataset, self.name, tmp=True)
        return {
            'bed': prefix / 'merged_cohorts.bed',
            'bim': prefix / 'merged_cohorts.bim',
            'fam': prefix / 'merged_cohorts.fam',
        }

    def queue_jobs(self, multicohort: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
        """
        Queue the multi-cohort PLINK 1.9 merge job, incorporating rolling aggregates if configured.
        """
        outputs = self.expected_outputs(multicohort)

        # 1. Gather cohort outputs
        all_cohort_outputs = inputs.as_dict_by_target(CohortBcfToPlink)
        cohort_plink_paths = []
        for _cohort_id, cohort_outs in all_cohort_outputs.items():
            cohort_plink_paths.append(
                {
                    'bed': str(cohort_outs['bed']),
                    'bim': str(cohort_outs['bim']),
                    'fam': str(cohort_outs['fam']),
                }
            )

        # 2. Check for rolling aggregate
        prev_analysis_id = config_retrieve(
            ['popgen_genotyping', 'rolling_aggregate', 'previous_analysis_id'], default=None
        )

        previous_aggregate_paths = None
        samples_to_remove = None

        if prev_analysis_id:
            previous_aggregate_paths, samples_to_remove = resolve_rolling_aggregate(prev_analysis_id)

        # 3. Call merge job
        j = run_merge_plink(
            cohort_plink_paths=cohort_plink_paths,
            output_prefix=str(outputs['bed']).replace('.bed', ''),
            previous_aggregate_paths=previous_aggregate_paths,
            samples_to_remove=samples_to_remove,
            job_name='MergeCohortPlink',
        )

        return self.make_outputs(multicohort, data=outputs, jobs=[j])


@stage(required_stages=[MergeCohortPlink])
class ExportCohortDatasets(MultiCohortStage):
    """
    Export the merged cohort to PLINK2 and BCF formats for long-term storage.
    """

    def expected_outputs(self, multicohort: 'MultiCohort') -> dict[str, 'Path']:
        """
        Define the expected PLINK2 and BCF outputs in long-term storage.
        """
        prefix = get_output_prefix(multicohort.analysis_dataset, self.name)
        return {
            'pgen': prefix / 'cohort.pgen',
            'pvar': prefix / 'cohort.pvar',
            'psam': prefix / 'cohort.psam',
            'bcf': prefix / 'cohort.bcf',
        }

    def queue_jobs(self, multicohort: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
        """
        Queue the dataset export job using PLINK2.
        """
        outputs = self.expected_outputs(multicohort)

        # 1. Pull input from MergeCohortPlink
        input_plink = inputs.as_dict(multicohort, MergeCohortPlink)

        # 2. Call export job
        j = run_export_cohort_datasets(
            input_plink_prefix={
                'bed': str(input_plink['bed']),
                'bim': str(input_plink['bim']),
                'fam': str(input_plink['fam']),
            },
            output_prefix=str(outputs['pgen']).replace('.pgen', ''),
            job_name='ExportCohortDatasets',
        )

        return self.make_outputs(multicohort, data=outputs, jobs=[j])
