"""
Stage-level tests for the two phase-2 stages rewritten for the two-phase run:
MergeCohortPlink (Metamist-resolved merge inputs, mandatory --keep trim) and
QcReport (BafRegress resolved for the full super-cohort membership).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

from popgen_genotyping.stages import KingIbdseg, MergeCohortPlink, Plink2Qc, QcReport

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
        """BafRegress is resolved for every super-cohort SG, not just the new plates."""
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
            bafregress_paths=list(bafregress_map.values()),
            output_path=str(report_path),
            job_name='QcReport_super_cohort',
        )
        mock_self.make_outputs.assert_called_once_with(mock_cohort, data=report_path, jobs=[mock_run.return_value])
