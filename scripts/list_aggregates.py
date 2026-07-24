#!/usr/bin/env python3

"""Ad hoc script to list registered array aggregate cohorts in a Metamist project.

Run locally; it only needs Metamist auth (e.g. ``gcloud auth application-default login``)
and does NOT read any cpg-flow pipeline config or run via analysis-runner:

    python scripts/list_aggregates.py
    python scripts/list_aggregates.py --analysis-type array_aggregate_pgen

Use it to pick the ``previous_aggregate_cohort_id`` for a rolling-aggregate (phase-2) run.

Aggregates register against a cohort (with empty sequencing_group_ids), so discovery is
cohort-first: the sample count and the rolling-forward handle both come from the cohort,
not the analysis. The analysis ID / output path are shown for reference only. An aggregate
normally maps to a single cohort; a legacy aggregate registered against multiple input
plate cohorts appears once per cohort, flagged ``[legacy]``.
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
    List cohorts carrying a registered aggregate analysis in a Metamist project.

    Args:
        project (str): Metamist project name. Defaults to 'ourdna'.
        analysis_type (str): Analysis type to list. Defaults to 'array_aggregate_pgen'.

    Returns:
        list[dict[str, Any]]: One entry per cohort that has an ``analysis_type`` analysis,
            newest first, each with keys ``cohort_id``, ``cohort_name``, ``num_sgs``,
            ``created`` (latest analysis timestamp), ``analysis_id`` and ``output`` (pgen
            path, both for reference), and ``legacy`` (True if that analysis is shared
            across multiple cohorts — an old-style multi-cohort registration).
    """
    result: dict[str, Any] = query(QUERY_AGGREGATE_ANALYSES, {'project': project})
    cohorts: list[dict[str, Any]] = (result.get('project') or {}).get('cohorts', [])

    # First pass: how many cohorts each aggregate analysis spans (for legacy detection).
    analysis_cohort_counts: dict[int, int] = {}
    for cohort in cohorts:
        for analysis in cohort.get('analyses', []):
            if analysis.get('type') == analysis_type:
                aid = analysis.get('id')
                analysis_cohort_counts[aid] = analysis_cohort_counts.get(aid, 0) + 1

    # Second pass: one row per cohort, using its latest aggregate analysis.
    entries: list[dict[str, Any]] = []
    for cohort in cohorts:
        aggregates = [a for a in cohort.get('analyses', []) if a.get('type') == analysis_type]
        if not aggregates:
            continue

        latest = max(aggregates, key=lambda a: (a.get('timestampCompleted') or '', a.get('id') or 0))
        outputs = latest.get('outputs')
        output_path = outputs.get('path') if isinstance(outputs, dict) else outputs

        entries.append(
            {
                'cohort_id': cohort.get('id'),
                'cohort_name': cohort.get('name'),
                'num_sgs': len(cohort.get('sequencingGroups') or []),
                'created': latest.get('timestampCompleted'),
                'analysis_id': latest.get('id'),
                'output': output_path,
                'legacy': analysis_cohort_counts.get(latest.get('id'), 0) > 1,
            }
        )

    # Newest first; fall back to cohort ID when timestamps are equal/missing.
    entries.sort(key=lambda e: (e['created'] or '', e['cohort_id'] or ''), reverse=True)
    return entries


def _format_table(aggregates: list[dict[str, Any]]) -> str:
    """
    Render the aggregate listing as an aligned plain-text table.

    Args:
        aggregates (list[dict[str, Any]]): Output of ``list_previous_aggregates``.

    Returns:
        str: The formatted table (header + one row per cohort).
    """
    headers = ('COHORT_ID', 'COHORT_NAME', 'N_SGS', 'CREATED', 'ANALYSIS_ID', 'OUTPUT')
    rows: list[tuple[str, str, str, str, str, str]] = []
    for agg in aggregates:
        cohort_name = agg['cohort_name'] or agg['cohort_id'] or '-'
        if agg['legacy']:
            cohort_name += ' [legacy]'
        rows.append(
            (
                str(agg['cohort_id'] or '-'),
                cohort_name,
                str(agg['num_sgs']),
                str(agg['created'] or '-'),
                str(agg['analysis_id'] or '-'),
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
