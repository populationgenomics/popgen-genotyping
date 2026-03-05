"""
This file exists to define all the Stages for the workflow.
"""

from typing import TYPE_CHECKING

from popgen_genotyping.jobs.gtc_to_heavy_vcf_job import run_gtc_to_heavy_vcf
from popgen_genotyping.jobs.bafregress_job import run_bafregress
from popgen_genotyping.jobs.heavy_to_light_vcf_job import run_heavy_to_light_vcf
from popgen_genotyping.utils import get_output_prefix

from cpg_utils.config import config_retrieve
from cpg_flow.workflow import get_workflow
from cpg_flow.stage import SequencingGroupStage, stage

if TYPE_CHECKING:
    from cpg_utils import Path
    from cpg_flow.targets import SequencingGroup, StageInput, StageOutput


@stage()
class GtcToHeavyVcf(SequencingGroupStage):
    """
    Convert GTC to VCF.
    """

    def expected_outputs(self, sequencing_group: 'SequencingGroup') -> 'Path':
        """
        Expected outputs for this stage.
        """
        return get_output_prefix(sequencing_group.dataset, self.name, tmp=True) / f'{sequencing_group.id}.bcf'

    def queue_jobs(self, sequencing_group: 'SequencingGroup', inputs: 'StageInput') -> 'StageOutput':  # noqa: ARG002
        """
        Queue jobs for this stage.
        """
        outputs = self.expected_outputs(sequencing_group)

        # Retrieve reference paths from config
        fasta_ref_path = config_retrieve(['popgen_genotyping', 'references', 'fasta_ref_path'])
        bpm_manifest_path = config_retrieve(['popgen_genotyping', 'references', 'bpm_manifest_path'])
        egt_cluster_path = config_retrieve(['popgen_genotyping', 'references', 'egt_cluster_path'])

        # TODO: Retrieve actual GTC and ID mapping paths from Metamist or metadata
        # Mocking for now
        gtc_path = f'gs://cpg-popgen-main/gtc/{sequencing_group.id}.gtc'
        id_mappings_path = f'gs://cpg-popgen-main/metadata/{sequencing_group.id}_id_map.txt'

        j = run_gtc_to_heavy_vcf(
            gtc_path=gtc_path,
            id_mappings_path=id_mappings_path,
            output_bcf_path=str(outputs),
            bpm_manifest_path=bpm_manifest_path,
            egt_cluster_path=egt_cluster_path,
            fasta_ref_path=fasta_ref_path,
            job_name=f'GtcToHeavyVcf_{sequencing_group.id}',
        )

        return self.make_outputs(sequencing_group, data=outputs, jobs=[j])


@stage(required_stages=[GtcToHeavyVcf])
class BafRegress(SequencingGroupStage):
    """
    Run BAFRegress on a BCF file to estimate sample contamination.
    """

    def expected_outputs(self, sequencing_group: 'SequencingGroup') -> 'Path':
        """
        Expected outputs for this stage.
        """
        return get_output_prefix(sequencing_group.dataset, self.name) / f'{sequencing_group.id}.BAFRegress.txt'

    def queue_jobs(self, sequencing_group: 'SequencingGroup', inputs: 'StageInput') -> 'StageOutput':
        """
        Queue jobs for this stage.
        """
        outputs = self.expected_outputs(sequencing_group)

        # Get the previous stage's output
        bcf_path = inputs.as_path(sequencing_group, GtcToHeavyVcf)

        j = run_bafregress(
            bcf_path=str(bcf_path),
            output_path=str(outputs),
            job_name=f'BafRegress_{sequencing_group.id}',
        )

        return self.make_outputs(sequencing_group, data=outputs, jobs=[j])


@stage(required_stages=[GtcToHeavyVcf])
class HeavyToLightVcf(SequencingGroupStage):
    """
    Strip intensities and FORMAT fields from a heavy BCF to produce a light BCF.
    """

    def expected_outputs(self, sequencing_group: 'SequencingGroup') -> 'Path':
        """
        Expected outputs for this stage.
        """
        return get_output_prefix(sequencing_group.dataset, self.name, tmp=True) / f'{sequencing_group.id}.bcf'

    def queue_jobs(self, sequencing_group: 'SequencingGroup', inputs: 'StageInput') -> 'StageOutput':
        """
        Queue jobs for this stage.
        """
        outputs = self.expected_outputs(sequencing_group)

        # Retrieve reference paths from config
        fasta_ref_path = config_retrieve(['popgen_genotyping', 'references', 'fasta_ref_path'])

        # Get the previous stage's output
        heavy_bcf_path = inputs.as_path(sequencing_group, GtcToHeavyVcf)

        # TODO: Retrieve actual ID mapping path from Metamist or metadata
        # Mocking for now
        id_mappings_path = f'gs://cpg-popgen-main/metadata/{sequencing_group.id}_id_map.txt'

        j = run_heavy_to_light_vcf(
            input_bcf_path=str(heavy_bcf_path),
            id_mappings_path=id_mappings_path,
            output_bcf_path=str(outputs),
            fasta_ref_path=fasta_ref_path,
            job_name=f'HeavyToLightVcf_{sequencing_group.id}',
        )

        return self.make_outputs(sequencing_group, data=outputs, jobs=[j])
