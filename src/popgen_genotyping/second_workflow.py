#!/usr/bin/env python3

"""
Phase-2 entry point: rolling merge, export and QC against the super cohort.

Run against the super cohort (``input_cohorts = [super]``). Normally submitted
automatically by phase 1's ``SubmitPhase2`` stage; can also be launched manually against a
hand-made super cohort.
"""

from argparse import ArgumentParser

from cpg_flow.workflow import run_workflow
from popgen_genotyping.stages import (
    ExportCohortDatasets,
    KingIbdseg,
    MergeCohortPlink,
    Plink2Qc,
    QcReport,
    SnpQcReport,
)


def cli_main() -> None:
    """
    Command line entry point for phase 2 of the genotyping pipeline.
    """
    parser = ArgumentParser(description='Genotyping microarray pipeline: phase 2 (aggregate)')
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    # The workflow name is derived from the package name
    workflow_name: str = __package__ or 'popgen_genotyping'
    stages: list = [
        MergeCohortPlink,
        ExportCohortDatasets,
        Plink2Qc,
        KingIbdseg,
        SnpQcReport,
        QcReport,
    ]

    run_workflow(name=workflow_name, stages=stages, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
