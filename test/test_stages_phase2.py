"""
Stage-level tests for the phase-2 stages: MergeCohortPlink (Metamist-resolved merge inputs,
mandatory --keep trim), RelatednessExclusions (BafRegress resolved for the full super-cohort
membership, required thresholds), QcReport (the three relatedness estimators and the exclusion
list), SnpQcReport (required thresholds), and the stage graph the relatedness stages form.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest
from cpg_utils.config import ConfigError
from cpg_utils.config import config_retrieve as real_config_retrieve

from popgen_genotyping.stages import (
    IbisRelatedness,
    KingIbdseg,
    KingRobust,
    MergeCohortPlink,
    Plink2Qc,
    QcReport,
    RelatednessExclusions,
    SnpQcReport,
)

# -- Helpers ------------------------------------------------------------------

SUPER_SG_IDS = ['CPG1', 'CPG2', 'CPG3']

MERGE_OUTPUTS = {
    'bed': Path('/merge/COH999_merged.bed'),
    'bim': Path('/merge/COH999_merged.bim'),
    'fam': Path('/merge/COH999_merged.fam'),
}


def _mock_merge_cohort() -> MagicMock:
    cohort = MagicMock()
    cohort.id = 'COH999'
    cohort.get_sequencing_group_ids.return_value = SUPER_SG_IDS
    return cohort


def _resolved(previous_aggregate_paths: str | None) -> dict:
    # Shape returned by resolve_merge_inputs; plate entries carry more keys than
    # the bed/bim/fam triple the merge job consumes.
    return {
        'plate_merge_list': [
            {
                'cohort_id': 'COHP1',
                'bed': 'gs://p1.bed',
                'bim': 'gs://p1.bim',
                'fam': 'gs://p1.fam',
                'new_count': len(SUPER_SG_IDS),
            },
        ],
        'previous_aggregate_paths': previous_aggregate_paths,
        'samples_to_remove': ['CPGX'],
        'super_cohort_size': len(SUPER_SG_IDS),
    }


# -- Tests: MergeCohortPlink.queue_jobs ----------------------------------------


class TestMergeCohortPlinkQueueJobs:
    """Wire-up between the Metamist-resolved plan and the merge job."""

    def test_rolling_aggregate_run(self) -> None:
        """With a previous aggregate: convert it to PLINK1, merge with it, trim to super."""
        mock_cohort = _mock_merge_cohort()
        mock_self = MagicMock()
        mock_self.expected_outputs.return_value = MERGE_OUTPUTS

        conversion_job = MagicMock(name='conversion_job')
        converted_resource = MagicMock(name='converted_resource')
        plink1_prefix = Path('/merge/plink2_to_plink1')

        with (
            patch('popgen_genotyping.stages.config_retrieve', return_value='COH123') as mock_config,
            patch(
                'popgen_genotyping.stages.resolve_merge_inputs',
                return_value=_resolved(previous_aggregate_paths='gs://prev/agg'),
            ) as mock_resolve,
            patch('popgen_genotyping.stages.get_output_prefix', return_value=plink1_prefix),
            patch(
                'popgen_genotyping.stages.run_plink2_to_plink1',
                return_value=(conversion_job, converted_resource),
            ) as mock_convert,
            patch('popgen_genotyping.stages.run_merge_plink') as mock_merge,
        ):
            MergeCohortPlink.queue_jobs(mock_self, mock_cohort, MagicMock())

        # The previous aggregate comes from config; the plan is resolved from the super
        # cohort's membership (the cohort is the source of truth, plates are derived).
        mock_config.assert_called_once_with(
            ['popgen_genotyping', 'merge_cohort_plink', 'previous_aggregate_cohort_id'], default=None
        )
        mock_resolve.assert_called_once_with(
            super_cohort_sg_ids=SUPER_SG_IDS,
            previous_aggregate_cohort_id='COH123',
        )

        # The PLINK2 aggregate is converted back to PLINK 1.9 at a cohort-keyed prefix.
        mock_convert.assert_called_once_with(
            pfile_prefix='gs://prev/agg',
            output_prefix=str(plink1_prefix / 'COH999_plink1'),
            job_name='Plink2ToPlink1',
        )

        # The merge receives only bed/bim/fam per plate, the converted previous aggregate,
        # and keep_samples = the full super-cohort membership (the mandatory trim).
        mock_merge.assert_called_once_with(
            cohort_plink_paths=[{'bed': 'gs://p1.bed', 'bim': 'gs://p1.bim', 'fam': 'gs://p1.fam'}],
            output_prefix='/merge/COH999_merged',
            keep_samples=SUPER_SG_IDS,
            previous_aggregate_resource=converted_resource,
            samples_to_remove=['CPGX'],
            job_name='MergeCohortPlink',
        )

        # The merge waits for the conversion, and the merge job is what the stage returns.
        mock_merge.return_value.depends_on.assert_called_once_with(conversion_job)
        mock_self.make_outputs.assert_called_once_with(mock_cohort, data=MERGE_OUTPUTS, jobs=[mock_merge.return_value])

    def test_bootstrap_run_has_no_previous_aggregate(self) -> None:
        """Without a previous aggregate: no PLINK2 conversion, no merge dependency."""
        mock_cohort = _mock_merge_cohort()
        mock_self = MagicMock()
        mock_self.expected_outputs.return_value = MERGE_OUTPUTS

        with (
            patch('popgen_genotyping.stages.config_retrieve', return_value=None),
            patch(
                'popgen_genotyping.stages.resolve_merge_inputs',
                return_value=_resolved(previous_aggregate_paths=None),
            ),
            patch('popgen_genotyping.stages.run_plink2_to_plink1') as mock_convert,
            patch('popgen_genotyping.stages.run_merge_plink') as mock_merge,
        ):
            MergeCohortPlink.queue_jobs(mock_self, mock_cohort, MagicMock())

        mock_convert.assert_not_called()
        assert mock_merge.call_args.kwargs['previous_aggregate_resource'] is None
        mock_merge.return_value.depends_on.assert_not_called()

    def test_merge_plan_is_logged_for_operator_confirmation(self) -> None:
        """The resolved plan reaches the driver log — the operator's pre-flight check."""
        mock_cohort = _mock_merge_cohort()
        mock_self = MagicMock()
        mock_self.expected_outputs.return_value = MERGE_OUTPUTS

        with (
            patch('popgen_genotyping.stages.config_retrieve', return_value=None),
            patch(
                'popgen_genotyping.stages.resolve_merge_inputs',
                return_value=_resolved(previous_aggregate_paths=None),
            ),
            patch('popgen_genotyping.stages.format_merge_plan', return_value='THE MERGE PLAN') as mock_plan,
            patch('popgen_genotyping.stages.logger') as mock_logger,
            patch('popgen_genotyping.stages.run_merge_plink'),
        ):
            MergeCohortPlink.queue_jobs(mock_self, mock_cohort, MagicMock())

        mock_plan.assert_called_once()
        mock_logger.info.assert_called_once_with('THE MERGE PLAN')


# -- Tests: QcReport.queue_jobs ------------------------------------------------


class TestQcReportQueueJobs:
    """Wire-up between upstream QC outputs, the BafRegress map, and the report job."""

    def test_passes_inputs_and_full_membership_bafregress(self) -> None:
        """BafRegress is resolved for every super-cohort SG, not just the new plates.

        The resolved map is sg_id -> plate-level file, so SGs sharing a plate map to the
        same file (CPG1/CPG2 below); the job must receive each plate file once — a repeat
        would multiply that plate's rows in the report's IID merge.
        """
        mock_cohort = MagicMock()
        mock_cohort.name = 'super_cohort'
        mock_cohort.get_sequencing_group_ids.return_value = SUPER_SG_IDS

        report_path = Path('/out/COH999_20260115_qc_report.csv')
        mock_self = MagicMock()
        mock_self.expected_outputs.return_value = report_path

        def as_path(target: object, stage: Any, key: str) -> Path:
            del target
            return {
                (Plink2Qc, 'smiss'): Path('/qc/COH999_20260115_qc.smiss'),
                (KingIbdseg, 'seg'): Path('/king/COH999_20260115_king.seg'),
                (KingRobust, 'kin0'): Path('/king_robust/COH999_20260115_king_robust.kin0.gz'),
                (IbisRelatedness, 'pairs'): Path('/ibis/COH999_20260115_ibis_relatedness_pairs.tsv.gz'),
                (RelatednessExclusions, 'excluded'): Path('/exclusions/COH999_20260115_relatedness.excluded.tsv'),
            }[(stage, key)]

        mock_inputs = MagicMock()
        mock_inputs.as_path.side_effect = as_path

        bafregress_map = {
            'CPG1': 'gs://baf/COHP1.BAFRegress.txt',
            'CPG2': 'gs://baf/COHP1.BAFRegress.txt',
            'CPG3': 'gs://baf/COHP2.BAFRegress.txt',
        }

        with (
            patch('popgen_genotyping.stages.resolve_bafregress_map', return_value=bafregress_map) as mock_baf,
            patch('popgen_genotyping.stages.run_qc_report') as mock_run,
        ):
            QcReport.queue_jobs(mock_self, mock_cohort, mock_inputs)

        mock_baf.assert_called_once_with(sg_ids=SUPER_SG_IDS)
        mock_run.assert_called_once_with(
            plink_qc_prefix='/qc/COH999_20260115_qc',
            king_seg_path='/king/COH999_20260115_king.seg',
            king_robust_kin0_path='/king_robust/COH999_20260115_king_robust.kin0.gz',
            ibis_pairs_path='/ibis/COH999_20260115_ibis_relatedness_pairs.tsv.gz',
            relatedness_excluded_path='/exclusions/COH999_20260115_relatedness.excluded.tsv',
            bafregress_paths=['gs://baf/COHP1.BAFRegress.txt', 'gs://baf/COHP2.BAFRegress.txt'],
            output_path=str(report_path),
            job_name='QcReport_super_cohort',
        )
        mock_self.make_outputs.assert_called_once_with(mock_cohort, data=report_path, jobs=[mock_run.return_value])


# -- Tests: RelatednessExclusions ----------------------------------------------


class TestRelatednessExclusionsExpectedOutputs:
    """Path layout that Metamist `analysis_keys=['excluded']` and the three estimators depend on."""

    def test_path_layout_and_naming(self) -> None:
        """The audit TSV and the plink2 remove list, cohort-keyed, under the standard prefix."""
        prefix = Path('/cohort/relatedness_exclusions')
        mock_cohort = MagicMock()
        mock_cohort.id = 'COH123'
        mock_self = MagicMock()
        mock_self.name = 'RelatednessExclusions'

        with patch('popgen_genotyping.stages.get_output_prefix', return_value=prefix) as mock_prefix:
            result = RelatednessExclusions.expected_outputs(mock_self, mock_cohort)

        assert result == {
            'excluded': prefix / 'COH123_relatedness.excluded_samples.tsv',
            'remove_list': prefix / 'COH123_relatedness.remove_list.txt',
        }
        mock_prefix.assert_called_once_with(dataset=mock_cohort.dataset, stage_name='RelatednessExclusions')


class TestRelatednessExclusionsQueueJobs:
    """Wire-up between Plink2Qc, the BafRegress map, the thresholds and the job factory."""

    def test_passes_smiss_full_membership_bafregress_and_thresholds(self) -> None:
        """.smiss maps through; BafRegress is resolved for all SGs and deduped; thresholds come from config."""
        mock_cohort = MagicMock()
        mock_cohort.name = 'my_cohort'
        mock_cohort.get_sequencing_group_ids.return_value = SUPER_SG_IDS

        def as_path(target: object, stage: object, key: str) -> Path:
            del target
            assert stage is Plink2Qc
            assert key == 'smiss'
            return Path('/qc/cohort.smiss')

        mock_inputs = MagicMock()
        mock_inputs.as_path.side_effect = as_path

        expected_outputs = {
            'excluded': Path('/out/COH1_relatedness.excluded_samples.tsv'),
            'remove_list': Path('/out/COH1_relatedness.remove_list.txt'),
        }
        mock_self = MagicMock()
        mock_self.expected_outputs.return_value = expected_outputs

        bafregress_map = {
            'CPG1': 'gs://baf/COHP1.BAFRegress.txt',
            'CPG2': 'gs://baf/COHP1.BAFRegress.txt',
            'CPG3': 'gs://baf/COHP2.BAFRegress.txt',
        }
        thresholds: dict[str, float] = {'contamination_max': 0.02, 'fmiss_max': 0.03}

        def config_retrieve(path: list[str]) -> float:
            # No default is passed: a missing threshold must raise.
            assert path[:2] == ['popgen_genotyping', 'relatedness_exclusions']
            return thresholds[path[2]]

        with (
            patch('popgen_genotyping.stages.resolve_bafregress_map', return_value=bafregress_map) as mock_baf,
            patch('popgen_genotyping.stages.config_retrieve', side_effect=config_retrieve),
            patch('popgen_genotyping.stages.run_relatedness_exclusions') as mock_run,
        ):
            RelatednessExclusions.queue_jobs(mock_self, mock_cohort, mock_inputs)

        # BafRegress is resolved for the full cohort membership and deduped before passing on.
        mock_baf.assert_called_once_with(sg_ids=SUPER_SG_IDS)

        mock_run.assert_called_once_with(
            bafregress_paths=['gs://baf/COHP1.BAFRegress.txt', 'gs://baf/COHP2.BAFRegress.txt'],
            smiss_path='/qc/cohort.smiss',
            contamination_max=0.02,
            fmiss_max=0.03,
            output_excluded_path='/out/COH1_relatedness.excluded_samples.tsv',
            output_remove_list_path='/out/COH1_relatedness.remove_list.txt',
            job_name='RelatednessExclusions_my_cohort',
        )
        mock_self.make_outputs.assert_called_once_with(mock_cohort, data=expected_outputs, jobs=[mock_run.return_value])

    def test_missing_threshold_raises(self) -> None:
        """A missing threshold is a ConfigError, not a default."""
        mock_cohort = MagicMock()
        mock_cohort.get_sequencing_group_ids.return_value = SUPER_SG_IDS
        mock_inputs = MagicMock()
        mock_inputs.as_path.return_value = Path('/qc/cohort.smiss')

        config = {'popgen_genotyping': {'relatedness_exclusions': {'contamination_max': 0.03}}}

        with (
            patch('popgen_genotyping.stages.resolve_bafregress_map', return_value={}),
            patch(
                'popgen_genotyping.stages.config_retrieve',
                side_effect=lambda key, **_kw: real_config_retrieve(key, config=config),
            ),
            patch('popgen_genotyping.stages.run_relatedness_exclusions') as mock_run,
            pytest.raises(ConfigError, match='fmiss_max'),
        ):
            RelatednessExclusions.queue_jobs(MagicMock(), mock_cohort, mock_inputs)
        mock_run.assert_not_called()


# -- Tests: SnpQcReport ---------------------------------------------------------


SNP_QC_REPORT_CONFIG = {
    'popgen_genotyping': {
        'snp_qc_report': {
            'egt_info_bcf_reference_key': 'illumina_microarray/egt_bcf',
            'egt_info_bcf_index_reference_key': 'illumina_microarray/egt_bcf_index',
            'thresholds': {
                'gentrain_min': 0.7,
                'cluster_sep_min': 0.4,
                'fmiss_max': 0.02,
                'hwe_p': 1e-5,
                'hwe_k': 0.001,
                'hwe_midp': True,
                'hwe_keep_fewhet': True,
            },
        },
    },
}


class TestSnpQcReportQueueJobs:
    """Every threshold and reference key is required: none has a code default."""

    def _queue(self, config: dict[str, Any]) -> MagicMock:
        mock_cohort = MagicMock()
        mock_cohort.name = 'my_cohort'
        mock_inputs = MagicMock()
        mock_inputs.as_path.return_value = Path('/agg/COH1.pgen')
        mock_self = MagicMock()
        mock_self.expected_outputs.return_value = {
            'audit_tsv': Path('/out/audit'),
            'inclusion_list': Path('/out/include'),
            'summary_tsv': Path('/out/summary'),
        }

        with (
            patch(
                'popgen_genotyping.stages.config_retrieve',
                side_effect=lambda key, **_kw: real_config_retrieve(key, config=config),
            ),
            patch('popgen_genotyping.stages.reference_path', side_effect=lambda key: f'gs://ref/{key}'),
            patch('popgen_genotyping.stages.run_snp_qc_report') as mock_run,
        ):
            SnpQcReport.queue_jobs(mock_self, mock_cohort, mock_inputs)
        return mock_run

    def test_configured_values_reach_the_job(self) -> None:
        """Thresholds come from config and the EGT files from the configured reference keys."""
        mock_run = self._queue(SNP_QC_REPORT_CONFIG)
        kwargs = mock_run.call_args.kwargs
        assert kwargs['gentrain_min'] == 0.7
        assert kwargs['fmiss_max'] == 0.02
        assert kwargs['hwe_midp'] is True
        assert kwargs['egt_info_bcf_path'] == 'gs://ref/illumina_microarray/egt_bcf'
        assert kwargs['egt_info_bcf_index_path'] == 'gs://ref/illumina_microarray/egt_bcf_index'

    @pytest.mark.parametrize(
        'missing',
        [
            ('thresholds', 'gentrain_min'),
            ('thresholds', 'cluster_sep_min'),
            ('thresholds', 'fmiss_max'),
            ('thresholds', 'hwe_p'),
            ('thresholds', 'hwe_k'),
            ('thresholds', 'hwe_midp'),
            ('thresholds', 'hwe_keep_fewhet'),
            ('egt_info_bcf_reference_key',),
            ('egt_info_bcf_index_reference_key',),
        ],
    )
    def test_missing_key_raises(self, missing: tuple[str, ...]) -> None:
        """Dropping any one key raises ConfigError instead of falling back to a default."""
        config = copy.deepcopy(SNP_QC_REPORT_CONFIG)
        section: Any = config['popgen_genotyping']['snp_qc_report']
        for part in missing[:-1]:
            section = section[part]
        del section[missing[-1]]
        with pytest.raises(ConfigError, match=missing[-1]):
            self._queue(config)


# -- Tests: the relatedness stage graph -------------------------------------------


def _instantiate(stage_decorator: Any) -> Any:
    """Build a stage object the way cpg-flow does, without a workflow or config behind it."""
    with (
        patch('cpg_flow.stage.get_workflow'),
        patch('cpg_flow.stage.get_config', return_value={'workflow': {}}),
    ):
        return stage_decorator()


# stage -> (required stages, Metamist analysis type, analysis keys)
STAGE_WIRING: dict[Any, tuple[list[str], str, list[str] | None]] = {
    RelatednessExclusions: (['Plink2Qc'], 'array_relatedness_exclusions', ['excluded']),
    KingIbdseg: (['MergeCohortPlink', 'RelatednessExclusions'], 'array_relatedness_ibdseg', ['seg', 'seg_x']),
    KingRobust: (
        ['ExportCohortDatasets', 'RelatednessExclusions'],
        'array_relatedness_king_robust',
        ['kin0'],
    ),
    IbisRelatedness: (
        ['ExportCohortDatasets', 'RelatednessExclusions', 'KingIbdseg'],
        'array_relatedness_ibis',
        ['pairs'],
    ),
    QcReport: (
        ['Plink2Qc', 'RelatednessExclusions', 'KingIbdseg', 'KingRobust', 'IbisRelatedness'],
        'array_qc_report',
        None,
    ),
}


class TestRelatednessStageGraph:
    """The dependencies and Metamist registration of the relatedness stages."""

    @pytest.mark.parametrize('stage_decorator', list(STAGE_WIRING), ids=lambda s: s.__name__)
    def test_required_stages_and_analysis(self, stage_decorator: Any) -> None:
        """Each stage requires exactly the stages whose outputs it consumes, and registers its analysis."""
        required, analysis_type, analysis_keys = STAGE_WIRING[stage_decorator]
        stage = _instantiate(stage_decorator)
        assert [dependency.__name__ for dependency in stage.required_stages_classes] == required
        assert stage.analysis_type == analysis_type
        assert stage.analysis_keys == analysis_keys

    def test_analysis_keys_are_expected_outputs(self) -> None:
        """Every registered analysis key is one of the stage's expected outputs, or cpg-flow fails the stage."""
        cohort = MagicMock()
        cohort.id = 'COH1'
        for stage_decorator, (_required, _analysis_type, analysis_keys) in STAGE_WIRING.items():
            if analysis_keys is None:
                continue
            stage = _instantiate(stage_decorator)
            with patch('popgen_genotyping.stages.get_output_prefix', return_value=Path('/p')):
                outputs = stage.expected_outputs(cohort)
            assert set(analysis_keys) <= set(outputs), stage_decorator.__name__

    def test_kingibdseg_no_longer_exposes_an_excluded_output(self) -> None:
        """The excluded-samples TSV moved to RelatednessExclusions."""
        stage = _instantiate(KingIbdseg)
        cohort = MagicMock()
        cohort.id = 'COH1'
        with patch('popgen_genotyping.stages.get_output_prefix', return_value=Path('/p')):
            assert 'excluded' not in stage.expected_outputs(cohort)

    def test_expected_outputs_do_not_collide_across_stages(self) -> None:
        """No two relatedness stages write the same path, since they share a cohort ID and version."""
        cohort = MagicMock()
        cohort.id = 'COH1'
        paths: list[Path] = []
        for stage_decorator in (RelatednessExclusions, KingIbdseg, KingRobust, IbisRelatedness):
            stage = _instantiate(stage_decorator)
            with patch(
                'popgen_genotyping.stages.get_output_prefix',
                side_effect=lambda stage_name, **_kw: Path(f'/prefix/{stage_name}'),
            ):
                paths.extend(stage.expected_outputs(cohort).values())
        assert len(paths) == len(set(paths))
