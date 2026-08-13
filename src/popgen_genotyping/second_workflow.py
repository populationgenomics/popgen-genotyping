#!/usr/bin/env python3

"""
Phase-2 entry point: rolling merge, export and QC against the super cohort.

Run against the super cohort (``input_cohorts = [super]``). Normally submitted
automatically by phase 1's ``SubmitPhase2`` stage; can also be launched manually against a
hand-made super cohort.
"""

from argparse import ArgumentParser

from cpg_flow.workflow import run_workflow
from cpg_utils.config import config_retrieve

from popgen_genotyping.stages import (
    ExportCohortDatasets,
    KingIbdseg,
    MergeCohortPlink,
    Plink2Qc,
    QcReport,
    SnpQcReport,
)
from popgen_genotyping.utils import validate_only_stages

# The aggregate stages this entry point submits. The per-plate stages live in
# first_workflow.py: the two phases run against different cohorts (new plates vs the
# super cohort) and must not share a submission — see the README.
PHASE_2_STAGES: list = [MergeCohortPlink, ExportCohortDatasets, Plink2Qc, KingIbdseg, SnpQcReport, QcReport]


def validate_phase2_cohorts(input_cohorts: list[str]) -> None:
    """
    Require exactly one input cohort: the super cohort.

    Every phase-2 stage is a CohortStage, so each listed cohort would otherwise get its
    own run, each treating its cohort as "the super cohort" and rolling forward the same
    ``previous_aggregate_cohort_id`` — registering multiple aggregates that claim the
    same lineage.

    Args:
        input_cohorts (list[str]): The ``workflow.input_cohorts`` config value.

    Raises:
        ValueError: If anything other than exactly one cohort is listed.
    """
    if len(input_cohorts) != 1:
        raise ValueError(
            f'Phase 2 runs against exactly one cohort (the super cohort), but '
            f'workflow.input_cohorts has {len(input_cohorts)}: {input_cohorts}. '
            'Create the super cohort first and list only its ID.'
        )


def cli_main() -> None:
    """
    Command line entry point for phase 2 of the genotyping pipeline.
    """
    parser = ArgumentParser(description='Genotyping microarray pipeline: phase 2 (aggregate)')
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    validate_only_stages(
        only_stages=config_retrieve(['workflow', 'only_stages'], default=[]),
        phase_stages=PHASE_2_STAGES,
        entry_point='second_workflow (phase 2)',
    )
    validate_phase2_cohorts(input_cohorts=config_retrieve(['workflow', 'input_cohorts'], default=[]))

    # The workflow name is derived from the package name
    workflow_name: str = __package__ or 'popgen_genotyping'
    run_workflow(name=workflow_name, stages=PHASE_2_STAGES, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
