#!/usr/bin/env python3

"""
Phase-1 entry point: per-plate processing, then super-cohort creation and phase-2 submission.

Run against the new plate cohorts (``input_cohorts = [new plates]``). Registers the durable
per-plate outputs, then the final stage creates the super cohort and submits
``second_workflow`` against it — no manual Swagger step or second launch needed. To stop
before the automatic hand-off (e.g. accumulating plates across several phase-1 runs), set
``workflow.last_stages = ['BafRegress', 'CohortBcfToPlink']``.
"""

from argparse import ArgumentParser

from cpg_flow.workflow import run_workflow
from popgen_genotyping.stages import (
    BafRegress,
    CohortBcfToPlink,
    GtcToBcfs,
    SubmitPhase2,
)


def cli_main() -> None:
    """
    Command line entry point for phase 1 of the genotyping pipeline.
    """
    parser = ArgumentParser(description='Genotyping microarray pipeline: phase 1 (per-plate)')
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    # The workflow name is derived from the package name
    workflow_name: str = __package__ or 'popgen_genotyping'
    stages: list = [
        GtcToBcfs,
        BafRegress,
        CohortBcfToPlink,
        SubmitPhase2,
    ]

    run_workflow(name=workflow_name, stages=stages, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
