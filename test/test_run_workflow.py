"""
Tests for the submission-time phase check in run_workflow.py.
"""

from unittest.mock import patch

import pytest

from popgen_genotyping import stages as stages_module
from popgen_genotyping.run_workflow import (
    PHASE_1_STAGES,
    PHASE_2_STAGES,
    validate_phase_stage_selection,
)
from popgen_genotyping.stages import CohortStage


def _validate_with_only_stages(only_stages: list[str]) -> None:
    with patch('popgen_genotyping.run_workflow.config_retrieve', return_value=only_stages):
        validate_phase_stage_selection()


class TestValidatePhaseStageSelection:
    """A submission must select stages from exactly one phase."""

    def test_phase_1_selection_passes(self) -> None:
        _validate_with_only_stages(['GtcToBcfs', 'BafRegress', 'CohortBcfToPlink'])

    def test_phase_2_selection_passes(self) -> None:
        _validate_with_only_stages(
            ['MergeCohortPlink', 'ExportCohortDatasets', 'Plink2Qc', 'KingIbdseg', 'SnpQcReport', 'QcReport']
        )

    def test_single_stage_rerun_passes(self) -> None:
        """Re-running one stage (e.g. just the QC report) is a valid phase-2 subset."""
        _validate_with_only_stages(['QcReport'])

    def test_missing_only_stages_raises(self) -> None:
        with pytest.raises(ValueError, match='only_stages is not set'):
            _validate_with_only_stages([])

    def test_mixed_phases_raises(self) -> None:
        with pytest.raises(ValueError, match='mixes phase-1'):
            _validate_with_only_stages(['CohortBcfToPlink', 'MergeCohortPlink'])

    def test_unknown_stage_raises(self) -> None:
        with pytest.raises(ValueError, match='unknown stages'):
            _validate_with_only_stages(['NotAStage'])

    def test_wrong_case_raises(self) -> None:
        """cpg-flow matches only_stages case-sensitively when skipping, so a
        lower-cased name would silently skip every stage — reject it here."""
        with pytest.raises(ValueError, match='unknown stages'):
            _validate_with_only_stages(['gtctobcfs'])


def test_phase_lists_cover_all_stages() -> None:
    """The phase split must cover every CohortStage the module defines, so a new
    stage cannot be added without assigning it to a phase."""
    # @stage wraps each class in a decorator function; the class survives as __wrapped__.
    defined = {
        name
        for name, obj in vars(stages_module).items()
        if isinstance(getattr(obj, '__wrapped__', None), type) and issubclass(obj.__wrapped__, CohortStage)
    }
    split = {cls.__name__ for cls in PHASE_1_STAGES + PHASE_2_STAGES}
    assert split == defined
