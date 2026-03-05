#!/usr/bin/env python3

"""
This is the main entry point for the workflow.
"""

from argparse import ArgumentParser

from popgen_genotyping.stages import GtcToHeavyVcf, BafRegress, HeavyToLightVcf

from cpg_flow.workflow import run_workflow


def cli_main():
    """
    CLI entrypoint - starts up the workflow
    """
    parser = ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    # The workflow name is derived from the package name
    workflow_name = __package__ or 'popgen_genotyping'
    stages = [GtcToHeavyVcf, BafRegress, HeavyToLightVcf]

    run_workflow(name=workflow_name, stages=stages, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
