#!/usr/bin/env python3

"""Ad hoc script to list registered array aggregate analyses in a Metamist project.

Run locally; it only needs Metamist auth (e.g. ``gcloud auth application-default login``)
and does NOT read any cpg-flow pipeline config or run via analysis-runner:

    python scripts/list_aggregates.py
    python scripts/list_aggregates.py --analysis-type array_aggregate_pgen

Use it to pick the ``previous_aggregate_analysis_id`` for a rolling-aggregate (phase-2) run.

Aggregates register against a cohort (with empty sequencing_group_ids), so the sample
count is read from the cohort, not the analysis. An aggregate normally maps to a single
cohort; legacy aggregates registered against multiple input plate cohorts map to several
(shown with a ``[legacy]`` marker).
"""

from __future__ import annotations

import argparse
from typing import Any

from metamist.graphql import gql, query

# The genotyping pipeline always runs against the ourdna dataset.
PROJECT = 'ourdna'

QUERY_AGGREGATE_ANALYSES = gql(
    """
    query AggregateAnalysesQuery($project: String!) {
      project(name: $project) {
        cohorts {
          id
          name
          sequencingGroups {
            id
          }
          analyses {
            id
            type
            timestampCompleted
            outputs
          }
        }
      }
    }
    """
)


def list_previous_aggregates(
    project: str = PROJECT,
    analysis_type: str = 'array_aggregate_pgen',
) -> list[dict[str, Any]]:
    """
    List registered aggregate analyses in a Metamist project, grouped by analysis ID.

    Args:
        project (str): Metamist project name. Defaults to 'ourdna'.
        analysis_type (str): Analysis type to list. Defaults to 'array_aggregate_pgen'.

    Returns:
        list[dict[str, Any]]: One entry per aggregate analysis, newest first, each with
            keys ``analysis_id``, ``created`` (timestampCompleted), ``output`` (pgen path),
            and ``cohorts`` (list of ``{'cohort_id', 'cohort_name', 'num_sgs'}``).
    """
    result: dict[str, Any] = query(QUERY_AGGREGATE_ANALYSES, {'project': project})
    cohorts: list[dict[str, Any]] = (result.get('project') or {}).get('cohorts', [])

    by_analysis: dict[int, dict[str, Any]] = {}
    for cohort in cohorts:
        num_sgs: int = len(cohort.get('sequencingGroups') or [])
        for analysis in cohort.get('analyses', []):
            if analysis.get('type') != analysis_type:
                continue
            analysis_id: int = analysis.get('id')
            outputs = analysis.get('outputs')
            output_path = outputs.get('path') if isinstance(outputs, dict) else outputs
            entry = by_analysis.setdefault(
                analysis_id,
                {
                    'analysis_id': analysis_id,
                    'created': analysis.get('timestampCompleted'),
                    'output': output_path,
                    'cohorts': [],
                },
            )
            entry['cohorts'].append(
                {'cohort_id': cohort.get('id'), 'cohort_name': cohort.get('name'), 'num_sgs': num_sgs}
            )

    aggregates: list[dict[str, Any]] = list(by_analysis.values())
    # Newest first; fall back to analysis ID when timestamps are equal/missing.
    aggregates.sort(key=lambda a: (a['created'] or '', a['analysis_id'] or 0), reverse=True)
    return aggregates


def _format_table(aggregates: list[dict[str, Any]]) -> str:
    """
    Render the aggregate listing as an aligned plain-text table.

    Args:
        aggregates (list[dict[str, Any]]): Output of ``list_previous_aggregates``.

    Returns:
        str: The formatted table (header + one row per aggregate).
    """
    headers = ('ANALYSIS_ID', 'CREATED', 'COHORT', 'N_SGS', 'OUTPUT')
    rows: list[tuple[str, str, str, str, str]] = []
    for agg in aggregates:
        cohorts = agg['cohorts']
        if len(cohorts) == 1:
            cohort_str = cohorts[0]['cohort_name'] or cohorts[0]['cohort_id'] or '-'
            sgs_str = str(cohorts[0]['num_sgs'])
        else:
            cohort_str = ', '.join((c['cohort_name'] or c['cohort_id'] or '?') for c in cohorts) + ' [legacy]'
            sgs_str = ', '.join(str(c['num_sgs']) for c in cohorts)
        rows.append(
            (
                str(agg['analysis_id']),
                str(agg['created'] or '-'),
                cohort_str,
                sgs_str,
                str(agg['output'] or '-'),
            )
        )

    widths = [len(h) for h in headers]
    for row in rows:
        widths = [max(w, len(cell)) for w, cell in zip(widths, row, strict=True)]

    def _line(cells: tuple[str, ...]) -> str:
        return '  '.join(cell.ljust(w) for cell, w in zip(cells, widths, strict=True))

    return '\n'.join([_line(headers), *[_line(row) for row in rows]])


def main() -> None:
    """Parse arguments, query Metamist, and print the aggregate listing."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--analysis-type',
        default='array_aggregate_pgen',
        help='Analysis type to list (default: array_aggregate_pgen).',
    )
    args = parser.parse_args()

    aggregates = list_previous_aggregates(analysis_type=args.analysis_type)
    if not aggregates:
        print(f'No {args.analysis_type} analyses found in project {PROJECT}.')
        return

    print(_format_table(aggregates))


if __name__ == '__main__':
    main()
