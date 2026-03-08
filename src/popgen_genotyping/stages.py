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
from popgen_genotyping.jobs.gtc_to_bcfs_job import run_gtc_to_bcfs
from popgen_genotyping.metamist_utils import query_genotyping_manifests, parse_genotyping_manifest
from popgen_genotyping.utils import get_output_prefix, get_sequencing_group_cohort

if TYPE_CHECKING:
    from cpg_flow.stage import StageInput, StageOutput
    from cpg_flow.targets import Cohort, SequencingGroup
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

        # Resolve GTC path from Metamist manifest
        # 1. Query all possible genotyping manifests
        all_manifests = query_genotyping_manifests(sequencing_group.dataset.name)

        # 2. Resolve the cohort for this sequencing group
        cohort = get_sequencing_group_cohort(sequencing_group)
        cohort_id = cohort.id

        # Find manifest with the cohort ID in its basename
        matching_manifests = [
            m for m in all_manifests if cohort_id in m.get('basename', '')
        ]

        if not matching_manifests:
            raise ValueError(
                f'No manifest found for cohort {cohort_id} in project {sequencing_group.dataset.name}'
            )

        # Sort by ID to get the latest if multiple exist for the same cohort
        matching_manifests.sort(key=lambda x: x.get('id', 0), reverse=True)
        latest_manifest = matching_manifests[0]
        manifest_path = latest_manifest.get('path')

        if not manifest_path:
            raise ValueError(f'Manifest for cohort {cohort_id} has no path output')

        # 3. Parse the chosen manifest
        mapping = parse_genotyping_manifest(manifest_path)

        if sequencing_group.id not in mapping:
            raise ValueError(
                f'Sequencing group {sequencing_group.id} not found in genotyping manifest '
                f'{manifest_path}'
            )
        gtc_path = mapping[sequencing_group.id]

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

        # Pull light BCFs from GtcToBcfs for all SGs
        # as_dict_by_sg returns dict[SequencingGroup, dict[str, Path]]
        sg_to_outputs = inputs.as_dict_by_sg(GtcToBcfs)
        manifest_data = {
            'manifest': {sg.id: str(outputs['light_bcf']) for sg, outputs in sg_to_outputs.items()}
        }

        # Write manifest to a tmp location on cloud storage
        manifest_path = cohort.analysis_dataset.prefix(category='tmp') / self.name / 'manifest.json'
        with to_path(manifest_path).open('w') as f:
            json.dump(manifest_data, f)

        # Define the job via the job utility
        j = run_cohort_bcf_to_plink(
            manifest_path=str(manifest_path),
            output_prefix=str(outputs['pgen']).replace('.pgen', ''),
            job_name=f'CohortBcfToPlink_{cohort.name}',
        )

        return self.make_outputs(cohort, data=outputs, jobs=[j])
