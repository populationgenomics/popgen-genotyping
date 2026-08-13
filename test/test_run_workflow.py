"""
Tests for the submission-time phase check in run_workflow.py.
"""

import pytest

from popgen_genotyping import stages as stages_module
from popgen_genotyping.run_workflow import (
    PHASE_1_STAGES,
    PHASE_2_STAGES,
    validate_phase_stage_selection,
)
from popgen_genotyping.stages import CohortStage

PHASE_1_NAMES = ['GtcToBcfs', 'BafRegress', 'CohortBcfToPlink']
PHASE_2_NAMES = ['MergeCohortPlink', 'ExportCohortDatasets', 'Plink2Qc', 'KingIbdseg', 'SnpQcReport', 'QcReport']


class TestValidatePhaseStageSelection:
    """A submission must select stages from exactly one phase."""

    def test_phase_1_selection_passes(self) -> None:
        validate_phase_stage_selection(only_stages=PHASE_1_NAMES, input_cohorts=['COH101', 'COH102'])

    def test_phase_2_selection_passes(self) -> None:
        validate_phase_stage_selection(only_stages=PHASE_2_NAMES, input_cohorts=['COH200'])

    def test_single_stage_rerun_passes(self) -> None:
        """Re-running one stage (e.g. just the QC report) is a valid phase-2 subset."""
        validate_phase_stage_selection(only_stages=['QcReport'], input_cohorts=['COH200'])

    def test_missing_only_stages_raises(self) -> None:
        with pytest.raises(ValueError, match='only_stages is not set'):
            validate_phase_stage_selection(only_stages=[], input_cohorts=['COH200'])

    def test_mixed_phases_raises(self) -> None:
        with pytest.raises(ValueError, match='mixes phase-1'):
            validate_phase_stage_selection(
                only_stages=['CohortBcfToPlink', 'MergeCohortPlink'], input_cohorts=['COH200']
            )

    def test_unknown_stage_raises(self) -> None:
        with pytest.raises(ValueError, match='unknown stages'):
            validate_phase_stage_selection(only_stages=['NotAStage'], input_cohorts=['COH200'])

    def test_wrong_case_raises(self) -> None:
        """cpg-flow accepts wrong-cased only_stages names but skips stages by exact
        match, so in a mixed-case list the wrong-cased stage is silently skipped while
        the rest run (an all-lowercase list at least dies with 'No stages to run').
        Rejecting on exact case closes the silent-partial-skip case."""
        with pytest.raises(ValueError, match='unknown stages'):
            validate_phase_stage_selection(only_stages=['gtctobcfs'], input_cohorts=['COH101'])


class TestPhase2CohortCardinality:
    """Phase 2 runs against exactly one cohort: the super cohort."""

    def test_multiple_phase_2_cohorts_raises(self) -> None:
        """Two cohorts would each be treated as 'the super cohort', registering two
        aggregates that both claim the same previous-aggregate lineage."""
        with pytest.raises(ValueError, match='exactly one cohort'):
            validate_phase_stage_selection(only_stages=PHASE_2_NAMES, input_cohorts=['COH200', 'COH201'])

    def test_no_phase_2_cohorts_raises(self) -> None:
        with pytest.raises(ValueError, match='exactly one cohort'):
            validate_phase_stage_selection(only_stages=PHASE_2_NAMES, input_cohorts=[])

    def test_multiple_phase_1_cohorts_pass(self) -> None:
        """Phase 1 legitimately processes many plate cohorts in one submission."""
        validate_phase_stage_selection(only_stages=PHASE_1_NAMES, input_cohorts=['COH101', 'COH102', 'COH103'])


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
