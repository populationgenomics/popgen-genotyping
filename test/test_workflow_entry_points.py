"""
Tests for the submission-time checks in the phase entry points
(first_workflow / second_workflow).
"""

import pytest

from popgen_genotyping import stages as stages_module
from popgen_genotyping.first_workflow import PHASE_1_STAGES
from popgen_genotyping.second_workflow import PHASE_2_STAGES, validate_phase2_cohorts
from popgen_genotyping.stages import CohortStage, MultiCohortStage
from popgen_genotyping.utils import validate_only_stages

PHASE_2_NAMES = ['MergeCohortPlink', 'ExportCohortDatasets', 'Plink2Qc', 'KingIbdseg', 'SnpQcReport', 'QcReport']


class TestValidateOnlyStages:
    """only_stages may only name stages of the submitting entry point's phase."""

    def test_full_phase_selection_passes(self) -> None:
        validate_only_stages(only_stages=PHASE_2_NAMES, phase_stages=PHASE_2_STAGES, entry_point='second_workflow')

    def test_empty_selection_passes(self) -> None:
        """No only_stages runs the whole phase — the entry point pins the stage list."""
        validate_only_stages(only_stages=[], phase_stages=PHASE_1_STAGES, entry_point='first_workflow')

    def test_single_stage_rerun_passes(self) -> None:
        """Re-running one stage (e.g. just the QC report) is a valid within-phase subset."""
        validate_only_stages(only_stages=['QcReport'], phase_stages=PHASE_2_STAGES, entry_point='second_workflow')

    def test_other_phase_stage_raises(self) -> None:
        """A phase-1 stage in a phase-2 submission would be silently skipped by cpg-flow."""
        with pytest.raises(ValueError, match='does not run'):
            validate_only_stages(
                only_stages=['CohortBcfToPlink', 'MergeCohortPlink'],
                phase_stages=PHASE_2_STAGES,
                entry_point='second_workflow',
            )

    def test_unknown_stage_raises(self) -> None:
        with pytest.raises(ValueError, match='does not run'):
            validate_only_stages(only_stages=['NotAStage'], phase_stages=PHASE_1_STAGES, entry_point='first_workflow')

    def test_wrong_case_raises(self) -> None:
        """cpg-flow accepts wrong-cased only_stages names but skips stages by exact
        match, so in a mixed-case list the wrong-cased stage is silently skipped while
        the rest run (an all-lowercase list at least dies with 'No stages to run').
        Rejecting on exact case closes the silent-partial-skip case."""
        with pytest.raises(ValueError, match='does not run'):
            validate_only_stages(only_stages=['gtctobcfs'], phase_stages=PHASE_1_STAGES, entry_point='first_workflow')


class TestPhase2CohortCardinality:
    """Phase 2 runs against exactly one cohort: the super cohort."""

    def test_single_cohort_passes(self) -> None:
        validate_phase2_cohorts(input_cohorts=['COH200'])

    def test_multiple_cohorts_raises(self) -> None:
        """Two cohorts would each be treated as 'the super cohort', registering two
        aggregates that both claim the same previous-aggregate lineage."""
        with pytest.raises(ValueError, match='exactly one cohort'):
            validate_phase2_cohorts(input_cohorts=['COH200', 'COH201'])

    def test_no_cohorts_raises(self) -> None:
        with pytest.raises(ValueError, match='exactly one cohort'):
            validate_phase2_cohorts(input_cohorts=[])


def test_phase_lists_cover_all_stages() -> None:
    """The two entry points must cover every stage the module defines, so a new
    stage cannot be added without assigning it to a phase."""
    # @stage wraps each class in a decorator function; the class survives as __wrapped__.
    defined = {
        name
        for name, obj in vars(stages_module).items()
        if isinstance(getattr(obj, '__wrapped__', None), type)
        and issubclass(obj.__wrapped__, (CohortStage, MultiCohortStage))
    }
    split = {cls.__name__ for cls in PHASE_1_STAGES + PHASE_2_STAGES}
    assert split == defined
