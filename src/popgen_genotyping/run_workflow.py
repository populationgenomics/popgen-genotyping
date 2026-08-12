#!/usr/bin/env python3

"""
This is the main entry point for the workflow.
"""

from argparse import ArgumentParser

from cpg_flow.workflow import run_workflow
from cpg_utils.config import config_retrieve

from popgen_genotyping.stages import (
    BafRegress,
    CohortBcfToPlink,
    ExportCohortDatasets,
    GtcToBcfs,
    KingIbdseg,
    MergeCohortPlink,
    Plink2Qc,
    QcReport,
    SnpQcReport,
)

# The pipeline runs in two phases that must not share a submission (see the README):
# phase 1 processes the new plate cohorts; phase 2 aggregates against the manually
# created super cohort. Every stage is a CohortStage, so nothing in the stage graph
# separates them — the split below is what validate_phase_stage_selection enforces.
PHASE_1_STAGES: list = [GtcToBcfs, BafRegress, CohortBcfToPlink]
PHASE_2_STAGES: list = [MergeCohortPlink, ExportCohortDatasets, Plink2Qc, KingIbdseg, SnpQcReport, QcReport]


def validate_phase_stage_selection() -> None:
    """
    Require ``workflow.only_stages`` to select stages from exactly one phase.

    Without this, a phase-2 submission against the super cohort would also run the
    per-plate stages on it (registering plate-level outputs for a super cohort), and a
    phase-1 submission would run the aggregate stages once per plate. Failing here means
    a mismatched config dies at submission, before any job is queued.

    Raises:
        ValueError: If ``workflow.only_stages`` is missing, names an unknown stage, or
            mixes phase-1 and phase-2 stages.
    """
    phase_1_names = {cls.__name__ for cls in PHASE_1_STAGES}
    phase_2_names = {cls.__name__ for cls in PHASE_2_STAGES}

    only_stages: list[str] = config_retrieve(['workflow', 'only_stages'], default=[])
    if not only_stages:
        raise ValueError(
            'workflow.only_stages is not set. This pipeline runs in two phases that must not '
            f'share a submission: phase 1 ({sorted(phase_1_names)}) runs per-plate against the '
            f'new plate cohorts, phase 2 ({sorted(phase_2_names)}) aggregates against the super '
            'cohort. Start from config_phase1.toml or config_phase2.toml.'
        )

    selected = set(only_stages)
    unknown = selected - phase_1_names - phase_2_names
    if unknown:
        raise ValueError(
            f'workflow.only_stages names unknown stages {sorted(unknown)}; '
            f'known stages are {sorted(phase_1_names | phase_2_names)} (exact case).'
        )
    if selected & phase_1_names and selected & phase_2_names:
        raise ValueError(
            f'workflow.only_stages mixes phase-1 stages {sorted(selected & phase_1_names)} with '
            f'phase-2 stages {sorted(selected & phase_2_names)}. The phases run against different '
            'cohorts (new plates vs the super cohort) and must be submitted separately — start '
            'from config_phase1.toml or config_phase2.toml.'
        )


def cli_main() -> None:
    """
    Command line entry point for the genotyping pipeline.
    """
    parser = ArgumentParser(description='Genotyping microarray pipeline')
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    validate_phase_stage_selection()

    # The workflow name is derived from the package name
    workflow_name: str = __package__ or 'popgen_genotyping'
    run_workflow(name=workflow_name, stages=PHASE_1_STAGES + PHASE_2_STAGES, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
