#!/usr/bin/env python3

"""
This is the main entry point for the workflow.
"""

from argparse import ArgumentParser

from cpg_flow.workflow import run_workflow
from popgen_genotyping.stages import (
    BafRegress,
    CohortBcfToPlink,
    ExportCohortDatasets,
    GtcToBcfs,
    MergeCohortPlink,
    Plink2Qc,
)


def cli_main() -> None:
    """
    Command line entry point for the genotyping pipeline.
    """
    parser = ArgumentParser(description='Genotyping microarray pipeline')
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args: ArgumentParser = parser.parse_args()

    # The workflow name is derived from the package name
    workflow_name: str = __package__ or 'popgen_genotyping'
    stages: list = [GtcToBcfs, BafRegress, CohortBcfToPlink, MergeCohortPlink, ExportCohortDatasets, Plink2Qc]

    run_workflow(name=workflow_name, stages=stages, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
