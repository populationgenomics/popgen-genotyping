"""
This file exists to define all the Stages for the workflow.
"""

from typing import TYPE_CHECKING

from cpg_flow.stage import CohortStage, MultiCohortStage, SequencingGroupStage, stage
from cpg_utils.config import config_retrieve

from popgen_genotyping.jobs.bafregress_job import run_bafregress
from popgen_genotyping.jobs.cohort_bcf_to_plink_job import run_cohort_bcf_to_plink
from popgen_genotyping.jobs.export_cohort_datasets_job import run_export_cohort_datasets
from popgen_genotyping.jobs.gtc_to_bcfs_job import run_gtc_to_bcfs
from popgen_genotyping.jobs.merge_plink_job import run_merge_plink
from popgen_genotyping.metamist_utils import resolve_gtc_path, query_previous_aggregate, query_reported_sex
from popgen_genotyping.utils import get_output_prefix, parse_psam

if TYPE_CHECKING:
    from cpg_flow.stage import StageInput, StageOutput
    from cpg_flow.targets import Cohort, MultiCohort, SequencingGroup
    from cpg_utils import Path


@stage()
class GtcToBcfs(SequencingGroupStage):
    """
    Convert GTC to both Heavy and Light BCFs in a single job, and output metadata.
    """

    def expected_outputs(self, sequencing_group: 'SequencingGroup') -> dict[str, 'Path']:
        prefix_tmp = get_output_prefix(sequencing_group.dataset, self.name, tmp=True)
        prefix_main = get_output_prefix(sequencing_group.dataset, self.name)
        return {
            'heavy_bcf': prefix_tmp / f'{sequencing_group.id}.heavy.bcf',
            'light_bcf': prefix_tmp / f'{sequencing_group.id}.light.bcf',
            'metadata_tsv': prefix_main / f'{sequencing_group.id}_gtc_metadata.tsv',
        }

    def queue_jobs(self, sequencing_group: 'SequencingGroup', inputs: 'StageInput') -> 'StageOutput':  # noqa: ARG002
        outputs = self.expected_outputs(sequencing_group)
        fasta_ref_path = config_retrieve(['popgen_genotyping', 'references', 'fasta_ref_path'])
        bpm_manifest_path = config_retrieve(['popgen_genotyping', 'references', 'bpm_manifest_path'])
        egt_cluster_path = config_retrieve(['popgen_genotyping', 'references', 'egt_cluster_path'])

        # Resolve GTC path from Metamist manifest using the helper
        gtc_path = resolve_gtc_path(sequencing_group)

        j = run_gtc_to_bcfs(
            gtc_path=gtc_path,
            sample_id=sequencing_group.id,
            output_heavy_bcf_path=str(outputs['heavy_bcf']),
            output_light_bcf_path=str(outputs['light_bcf']),
            output_metadata_path=str(outputs['metadata_tsv']),
            bpm_manifest_path=bpm_manifest_path,
            egt_cluster_path=egt_cluster_path,
            fasta_ref_path=fasta_ref_path,
            job_name=f'GtcToBcfs_{sequencing_group.id}',
        )

        return self.make_outputs(sequencing_group, data=outputs, jobs=[j])


@stage(required_stages=[GtcToBcfs])
class BafRegress(SequencingGroupStage):
    """
    Run BAFRegress on a Heavy BCF file to estimate sample contamination.
    """

    def expected_outputs(self, sequencing_group: 'SequencingGroup') -> 'Path':
        return get_output_prefix(sequencing_group.dataset, self.name) / f'{sequencing_group.id}.BAFRegress.txt'

    def queue_jobs(self, sequencing_group: 'SequencingGroup', inputs: 'StageInput') -> 'StageOutput':
        outputs = self.expected_outputs(sequencing_group)
        # Pull heavy_bcf from GtcToBcfs
        bcf_path = inputs.as_dict(sequencing_group, GtcToBcfs)['heavy_bcf']

        j = run_bafregress(
            bcf_path=str(bcf_path),
            output_path=str(outputs),
            job_name=f'BafRegress_{sequencing_group.id}',
        )

        return self.make_outputs(sequencing_group, data=outputs, jobs=[j])


@stage(required_stages=[GtcToBcfs])
class CohortBcfToPlink(CohortStage):
    """
    Convert all light BCFs in a cohort to PLINK 1.9 format and merge them.
    Output is stored in tmp.
    """

    def expected_outputs(self, cohort: 'Cohort') -> dict[str, 'Path']:
        # Store in tmp per EDIT requirement
        prefix = get_output_prefix(cohort.analysis_dataset, self.name, tmp=True)
        return {
            'bed': prefix / 'cohort.bed',
            'bim': prefix / 'cohort.bim',
            'fam': prefix / 'cohort.fam',
        }

    def queue_jobs(self, cohort: 'Cohort', inputs: 'StageInput') -> 'StageOutput':
        outputs = self.expected_outputs(cohort)

        # 1. Pull light BCF paths from GtcToBcfs for all SGs in this cohort
        sg_outputs = inputs.as_dict_by_target(GtcToBcfs)
        cohort_sg_ids = set(cohort.get_sequencing_group_ids())

        bcf_paths = {
            sg_id: str(outs['light_bcf'])
            for sg_id, outs in sg_outputs.items()
            if sg_id in cohort_sg_ids
        }

        # 2. Fetch reported sex metadata for these SGs
        full_sex_mapping = query_reported_sex(cohort.analysis_dataset.name)
        # Filter to include only those in the current cohort
        sex_mapping = {
            sg_id: sex_code
            for sg_id, sex_code in full_sex_mapping.items()
            if sg_id in cohort_sg_ids
        }

        # Define the job via the job utility
        j = run_cohort_bcf_to_plink(
            bcf_paths=bcf_paths,
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
        # Store in tmp per EDIT requirement
        prefix = get_output_prefix(multicohort.analysis_dataset, self.name, tmp=True)
        return {
            'bed': prefix / 'merged_cohorts.bed',
            'bim': prefix / 'merged_cohorts.bim',
            'fam': prefix / 'merged_cohorts.fam',
        }

    def queue_jobs(self, multicohort: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
        outputs = self.expected_outputs(multicohort)

        # 1. Gather new cohort outputs
        all_cohort_outputs = inputs.as_dict_by_target(CohortBcfToPlink)
        cohort_plink_paths = []
        for _cohort_id, cohort_outs in all_cohort_outputs.items():
            cohort_plink_paths.append({
                'bed': str(cohort_outs['bed']),
                'bim': str(cohort_outs['bim']),
                'fam': str(cohort_outs['fam']),
            })

        # 2. Check for rolling aggregate
        prev_analysis_id = config_retrieve(
            ['popgen_genotyping', 'rolling_aggregate', 'previous_analysis_id'],
            default=None
        )

        previous_aggregate_paths = None
        samples_to_remove = None

        if prev_analysis_id:
            # Query Metamist for the previous analysis
            prev_outputs, active_sg_ids = query_previous_aggregate(int(prev_analysis_id))
            previous_aggregate_paths = {
                'bed': prev_outputs['bed'],
                'bim': prev_outputs['bim'],
                'fam': prev_outputs['fam'],
            }

            # Parse the previous PSAM/FAM to find all samples
            # Note: We might need a parse_fam utility if parse_psam is too specific to PLINK2
            prev_samples = parse_psam(previous_aggregate_paths['fam'])

            # Find samples that are no longer active
            samples_to_remove = list(set(prev_samples) - set(active_sg_ids))

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
        prefix = get_output_prefix(multicohort.analysis_dataset, self.name)
        return {
            'pgen': prefix / 'cohort.pgen',
            'pvar': prefix / 'cohort.pvar',
            'psam': prefix / 'cohort.psam',
            'bcf': prefix / 'cohort.bcf',
        }

    def queue_jobs(self, multicohort: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
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
