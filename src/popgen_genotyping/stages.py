"""
This file exists to define all the Stages for the workflow.
"""

from typing import TYPE_CHECKING

from cpg_flow.stage import CohortStage, MultiCohortStage, SequencingGroupStage, stage
from cpg_utils.config import config_retrieve

from popgen_genotyping.jobs.bafregress_job import run_bafregress
from popgen_genotyping.jobs.cohort_bcf_to_plink_job import run_cohort_bcf_to_plink
from popgen_genotyping.jobs.gtc_to_bcfs_job import run_gtc_to_bcfs
from popgen_genotyping.jobs.merge_plink_job import run_merge_plink
from popgen_genotyping.metamist_utils import resolve_gtc_path
from popgen_genotyping.utils import get_output_prefix

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
    Convert all light BCFs in a cohort to PLINK2 format and merge them.
    """

    def expected_outputs(self, cohort: 'Cohort') -> dict[str, 'Path']:
        prefix = get_output_prefix(cohort.analysis_dataset, self.name)
        return {
            'pgen': prefix / 'cohort.pgen',
            'pvar': prefix / 'cohort.pvar',
            'psam': prefix / 'cohort.psam',
        }

    def queue_jobs(self, cohort: 'Cohort', inputs: 'StageInput') -> 'StageOutput':
        outputs = self.expected_outputs(cohort)

        # Pull light BCF paths from GtcToBcfs for all SGs in this cohort
        sg_outputs = inputs.as_dict_by_target(GtcToBcfs)
        cohort_sg_ids = set(cohort.get_sequencing_group_ids())

        bcf_paths = {
            sg_id: str(outs['light_bcf'])
            for sg_id, outs in sg_outputs.items()
            if sg_id in cohort_sg_ids
        }

        # Define the job via the job utility
        j = run_cohort_bcf_to_plink(
            bcf_paths=bcf_paths,
            output_prefix=str(outputs['pgen']).replace('.pgen', ''),
            job_name=f'CohortBcfToPlink_{cohort.name}',
        )

        return self.make_outputs(cohort, data=outputs, jobs=[j])


@stage(required_stages=[CohortBcfToPlink])
class MergeCohortPlink(MultiCohortStage):
    """
    Merge all cohort PLINK2 datasets into a single unified dataset.
    """

    def expected_outputs(self, multicohort: 'MultiCohort') -> dict[str, 'Path']:
        prefix = get_output_prefix(multicohort.analysis_dataset, self.name)
        return {
            'pgen': prefix / 'merged_cohorts.pgen',
            'pvar': prefix / 'merged_cohorts.pvar',
            'psam': prefix / 'merged_cohorts.psam',
        }

    def queue_jobs(self, multicohort: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
        outputs = self.expected_outputs(multicohort)

        # Gather all PLINK2 output paths from CohortBcfToPlink stage
        # as_dict_by_target for a MultiCohort stage with CohortStage dependencies
        # returns dict[target_id, outputs]
        all_cohort_outputs = inputs.as_dict_by_target(CohortBcfToPlink)

        cohort_plink_paths = []
        for _cohort_id, cohort_outs in all_cohort_outputs.items():
            # cohort_outs is a dict with 'pgen', 'pvar', 'psam'
            cohort_plink_paths.append({
                'pgen': str(cohort_outs['pgen']),
                'pvar': str(cohort_outs['pvar']),
                'psam': str(cohort_outs['psam']),
            })

        j = run_merge_plink(
            cohort_plink_paths=cohort_plink_paths,
            output_prefix=str(outputs['pgen']).replace('.pgen', ''),
            job_name='MergeCohortPlink',
        )

        return self.make_outputs(multicohort, data=outputs, jobs=[j])
