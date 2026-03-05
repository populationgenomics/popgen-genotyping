"""
This file exists to define all the Stages for the workflow.
"""

import json
from typing import TYPE_CHECKING

from cpg_flow.stage import CohortStage, SequencingGroupStage, stage
from cpg_utils import to_path
from cpg_utils.config import config_retrieve

from popgen_genotyping.jobs.bafregress_job import run_bafregress
from popgen_genotyping.jobs.cohort_bcf_to_plink_job import run_cohort_bcf_to_plink
from popgen_genotyping.jobs.gtc_to_heavy_vcf_job import run_gtc_to_heavy_vcf
from popgen_genotyping.jobs.heavy_to_light_vcf_job import run_heavy_to_light_vcf
from popgen_genotyping.utils import get_output_prefix

if TYPE_CHECKING:
    from cpg_flow.stage import StageInput, StageOutput
    from cpg_flow.targets import Cohort, SequencingGroup
    from cpg_utils import Path


@stage()
class GtcToHeavyBcf(SequencingGroupStage):
    """
    Convert GTC to BCF (Heavy).
    """

    def expected_outputs(self, sequencing_group: 'SequencingGroup') -> 'Path':
        return get_output_prefix(sequencing_group.dataset, self.name, tmp=True) / f'{sequencing_group.id}.bcf'

    def queue_jobs(self, sequencing_group: 'SequencingGroup', inputs: 'StageInput') -> 'StageOutput':  # noqa: ARG002
        outputs = self.expected_outputs(sequencing_group)
        fasta_ref_path = config_retrieve(['popgen_genotyping', 'references', 'fasta_ref_path'])
        bpm_manifest_path = config_retrieve(['popgen_genotyping', 'references', 'bpm_manifest_path'])
        egt_cluster_path = config_retrieve(['popgen_genotyping', 'references', 'egt_cluster_path'])

        # TODO: Read in the manifest file for the cohort to extract the correct path for the GTC file instead of assuming it's in a fixed location with the sample ID as the filename
        gtc_path = f'gs://cpg-popgen-main/gtc/{sequencing_group.id}.gtc'

        j = run_gtc_to_heavy_vcf(
            gtc_path=gtc_path,
            sample_id=sequencing_group.id,
            output_bcf_path=str(outputs),
            bpm_manifest_path=bpm_manifest_path,
            egt_cluster_path=egt_cluster_path,
            fasta_ref_path=fasta_ref_path,
            job_name=f'GtcToHeavyBcf_{sequencing_group.id}',
        )

        return self.make_outputs(sequencing_group, data=outputs, jobs=[j])


@stage(required_stages=[GtcToHeavyBcf])
class BafRegress(SequencingGroupStage):
    """
    Run BAFRegress on a BCF file to estimate sample contamination.
    """

    def expected_outputs(self, sequencing_group: 'SequencingGroup') -> 'Path':
        return get_output_prefix(sequencing_group.dataset, self.name) / f'{sequencing_group.id}.BAFRegress.txt'

    def queue_jobs(self, sequencing_group: 'SequencingGroup', inputs: 'StageInput') -> 'StageOutput':
        outputs = self.expected_outputs(sequencing_group)
        bcf_path = inputs.as_path(sequencing_group, GtcToHeavyBcf)

        j = run_bafregress(
            bcf_path=str(bcf_path),
            output_path=str(outputs),
            job_name=f'BafRegress_{sequencing_group.id}',
        )

        return self.make_outputs(sequencing_group, data=outputs, jobs=[j])


@stage(required_stages=[GtcToHeavyBcf])
class HeavyToLightBcf(SequencingGroupStage):
    """
    Strip intensities from heavy BCF to produce a light BCF.
    """

    def expected_outputs(self, sequencing_group: 'SequencingGroup') -> 'Path':
        return get_output_prefix(sequencing_group.dataset, self.name, tmp=True) / f'{sequencing_group.id}.bcf'

    def queue_jobs(self, sequencing_group: 'SequencingGroup', inputs: 'StageInput') -> 'StageOutput':
        outputs = self.expected_outputs(sequencing_group)
        fasta_ref_path = config_retrieve(['popgen_genotyping', 'references', 'fasta_ref_path'])
        heavy_bcf_path = inputs.as_path(sequencing_group, GtcToHeavyBcf)

        j = run_heavy_to_light_vcf(
            input_bcf_path=str(heavy_bcf_path),
            output_bcf_path=str(outputs),
            fasta_ref_path=fasta_ref_path,
            job_name=f'HeavyToLightBcf_{sequencing_group.id}',
        )

        return self.make_outputs(sequencing_group, data=outputs, jobs=[j])


@stage(required_stages=[HeavyToLightBcf])
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

        # 1. Create manifest of {sg_id: bcf_path}
        bcf_paths = inputs.as_path_by_sg(HeavyToLightBcf)
        manifest_data = {'manifest': {sg.id: str(path) for sg, path in bcf_paths.items()}}

        # 2. Write manifest to a tmp location on cloud storage
        manifest_path = cohort.analysis_dataset.prefix(category='tmp') / self.name / 'manifest.json'
        with to_path(manifest_path).open('w') as f:
            json.dump(manifest_data, f)

        # 3. Define the job via the job utility
        j = run_cohort_bcf_to_plink(
            manifest_path=str(manifest_path),
            output_prefix=str(outputs['pgen']).replace('.pgen', ''),
            job_name=f'CohortBcfToPlink_{cohort.name}',
        )

        return self.make_outputs(cohort, data=outputs, jobs=[j])
