#!/usr/bin/env python3

"""Generate and register a genotyping-array manifest for a Metamist cohort.

Builds the manifest CSV that ``resolve_cohort_gtc_mapping`` expects (one row per
sequencing group, mapping SG ID to its GTC path and sentrix barcode/position),
uploads it next to the GTC data, and registers it as a ``manifest`` analysis so
the pipeline's manifest discovery (cohort ID in the basename) finds it.

In production the transfer pipeline produces this manifest; this script mimics it
for cohorts assembled by hand (e.g. test plate cohorts in ourdna-test), deriving
every field from the sequencing groups' Metamist meta
(``sentrix_barcode_a`` / ``sentrix_position_a`` / ``sample_plate``).

Run locally; it needs Metamist auth plus read access to the GTC bucket:

    python scripts/generate_genotyping_manifest.py --cohort COH123 --project ourdna-test
    python scripts/generate_genotyping_manifest.py --cohort COH123 --project ourdna-test --dry-run

Every GTC path is checked for existence before anything is written; a missing
GTC, missing SG meta, or an SG that is not a genotyping array fails the whole
run. An existing manifest for the cohort is never overwritten without --force.
"""

from __future__ import annotations

import argparse
import csv
import io
from typing import Any

from cpg_utils import to_path
from metamist.apis import AnalysisApi
from metamist.graphql import gql, query
from metamist.model.analysis import Analysis
from metamist.model.analysis_status import AnalysisStatus

# Column layout of the transfer pipeline's manifest. parse_genotyping_manifest_for_reheader
# only reads cpg_sequencing_group_id / cpg_gcp_filepath / sentrix_barcode_a /
# sentrix_position_a, but we keep the full header so generated manifests are
# drop-in identical in shape to production ones.
MANIFEST_COLUMNS = (
    'file_name',
    'sample_sheet_id',
    'sample_id',
    'md5sum',
    'bbv_barcode',
    'sentrix_barcode_a',
    'sentrix_position_a',
    'sample_plate',
    'sample_well',
    'sample_plate_position',
    'cpg_sample_id_internal',
    'cpg_gcp_filepath',
    'cpg_sequencing_group_id',
    'cpg_cohort_id',
)

# SG meta keys a row cannot be built without.
REQUIRED_META_KEYS = ('sentrix_barcode_a', 'sentrix_position_a', 'sample_plate')

QUERY_COHORT_SGS = gql(
    """
    query CohortSequencingGroups($project: String!, $cohortId: String!) {
      project(name: $project) {
        cohorts(id: {eq: $cohortId}) {
          id
          name
          sequencingGroups {
            id
            type
            meta
            sample {
              id
              externalId
            }
          }
        }
      }
    }
    """
)


def fetch_cohort_sequencing_groups(project: str, cohort_id: str) -> list[dict[str, Any]]:
    """
    Fetch the cohort's sequencing groups (with meta) from Metamist.

    Args:
        project (str): Metamist project name.
        cohort_id (str): Cohort ID (COH...).

    Returns:
        list[dict[str, Any]]: The cohort's sequencing group records.

    Raises:
        ValueError: If the cohort does not exist in the project or has no
            sequencing groups.
    """
    result: dict[str, Any] = query(QUERY_COHORT_SGS, {'project': project, 'cohortId': cohort_id})
    cohorts: list[dict[str, Any]] = (result.get('project') or {}).get('cohorts', [])
    if not cohorts:
        raise ValueError(f'Cohort {cohort_id} not found in project {project}')

    sequencing_groups: list[dict[str, Any]] = cohorts[0].get('sequencingGroups') or []
    if not sequencing_groups:
        raise ValueError(f'Cohort {cohort_id} in project {project} has no sequencing groups')
    return sequencing_groups


def build_manifest_rows(
    sequencing_groups: list[dict[str, Any]],
    cohort_id: str,
    gtc_prefix: str,
) -> list[dict[str, str]]:
    """
    Build one manifest row per sequencing group from its Metamist meta.

    Args:
        sequencing_groups (list[dict[str, Any]]): SG records from
            ``fetch_cohort_sequencing_groups``.
        cohort_id (str): Cohort ID written into every row.
        gtc_prefix (str): GCS prefix holding the per-plate GTC directories
            (e.g. ``gs://cpg-ourdna-test/gtc_genotyping_array``).

    Returns:
        list[dict[str, str]]: Manifest rows keyed by ``MANIFEST_COLUMNS``,
            sorted by (plate, barcode, position).

    Raises:
        ValueError: If any SG is not a genotyping array, is missing required
            meta, or shares a sentrix well with another SG in the cohort.
    """
    wrong_type: list[str] = [sg['id'] for sg in sequencing_groups if sg.get('type') != 'genotypingarray']
    if wrong_type:
        raise ValueError(f'Cohort {cohort_id} contains non-genotypingarray sequencing groups: {sorted(wrong_type)}')

    missing_meta: dict[str, list[str]] = {}
    for sg in sequencing_groups:
        meta = sg.get('meta') or {}
        absent = [key for key in REQUIRED_META_KEYS if not meta.get(key)]
        if absent:
            missing_meta[sg['id']] = absent
    if missing_meta:
        raise ValueError(f'Sequencing groups missing required meta {REQUIRED_META_KEYS}: {missing_meta}')

    rows: list[dict[str, str]] = []
    well_to_sg: dict[str, str] = {}
    for sg in sequencing_groups:
        meta = sg['meta']
        barcode: str = meta['sentrix_barcode_a']
        position: str = meta['sentrix_position_a']
        well = f'{barcode}_{position}'

        if well in well_to_sg:
            raise ValueError(f'Sentrix well {well} is shared by sequencing groups {well_to_sg[well]} and {sg["id"]}')
        well_to_sg[well] = sg['id']

        rows.append(
            {
                'file_name': f'{well}.gtc',
                'sample_sheet_id': '',
                'sample_id': '',
                'md5sum': '',
                'bbv_barcode': meta.get('bbv_barcode', ''),
                'sentrix_barcode_a': barcode,
                'sentrix_position_a': position,
                'sample_plate': meta['sample_plate'],
                'sample_well': meta.get('sample_well', well),
                'sample_plate_position': meta.get('sample_plate_position', ''),
                'cpg_sample_id_internal': (sg.get('sample') or {}).get('id', ''),
                'cpg_gcp_filepath': f'{gtc_prefix}/{meta["sample_plate"]}/{well}.gtc',
                'cpg_sequencing_group_id': sg['id'],
                'cpg_cohort_id': cohort_id,
            }
        )

    rows.sort(key=lambda r: (r['sample_plate'], r['sentrix_barcode_a'], r['sentrix_position_a']))
    return rows


def check_gtcs_exist(rows: list[dict[str, str]]) -> None:
    """
    Verify every row's GTC file exists in GCS and is non-empty.

    Args:
        rows (list[dict[str, str]]): Manifest rows from ``build_manifest_rows``.

    Raises:
        ValueError: If any GTC file is missing or zero bytes, listing every bad path.
    """
    missing: list[str] = []
    empty: list[str] = []
    for row in rows:
        path = to_path(row['cpg_gcp_filepath'])
        if not path.exists():
            missing.append(row['cpg_gcp_filepath'])
        elif path.stat().st_size == 0:
            empty.append(row['cpg_gcp_filepath'])

    problems: list[str] = []
    if missing:
        problems.append(f'{len(missing)} GTC file(s) not found:\n' + '\n'.join(missing))
    if empty:
        problems.append(f'{len(empty)} GTC file(s) are zero bytes:\n' + '\n'.join(empty))
    if problems:
        raise ValueError('\n'.join(problems))


def render_manifest_csv(rows: list[dict[str, str]]) -> str:
    """
    Render manifest rows as CSV text.

    Args:
        rows (list[dict[str, str]]): Manifest rows from ``build_manifest_rows``.

    Returns:
        str: CSV content including the header line.
    """
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=MANIFEST_COLUMNS, lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def register_manifest_analysis(project: str, cohort_id: str, manifest_path: str) -> int:
    """
    Register the uploaded manifest as a completed ``manifest`` analysis.

    Args:
        project (str): Metamist project name.
        cohort_id (str): Cohort the manifest belongs to.
        manifest_path (str): GCS path of the uploaded manifest CSV.

    Returns:
        int: The created analysis ID.
    """
    analysis = Analysis(
        type='manifest',
        status=AnalysisStatus('completed'),
        output=manifest_path,
        sequencing_group_ids=[],
        cohort_ids=[cohort_id],
        meta={},
    )
    return AnalysisApi().create_analysis(project, analysis)


def main() -> None:
    """Parse arguments, build the manifest, and (unless --dry-run) upload and register it."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--cohort', required=True, help='Cohort ID (COH...) to build the manifest for.')
    parser.add_argument('--project', required=True, help='Metamist project name (e.g. ourdna-test).')
    parser.add_argument(
        '--gtc-prefix',
        default=None,
        help='GCS prefix of the per-plate GTC directories (default: gs://cpg-<project>/gtc_genotyping_array).',
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print the manifest CSV instead of uploading and registering it.',
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite an existing manifest CSV for this cohort.',
    )
    args = parser.parse_args()

    gtc_prefix: str = (args.gtc_prefix or f'gs://cpg-{args.project}/gtc_genotyping_array').rstrip('/')
    manifest_path: str = f'{gtc_prefix}/manifests/genotyping_array_manifest_cohort_{args.cohort}.csv'

    sequencing_groups = fetch_cohort_sequencing_groups(args.project, args.cohort)
    rows = build_manifest_rows(sequencing_groups, args.cohort, gtc_prefix)
    check_gtcs_exist(rows)
    content = render_manifest_csv(rows)

    plates = sorted({row['sample_plate'] for row in rows})
    print(f'Cohort {args.cohort}: {len(rows)} sequencing groups across plate(s) {", ".join(plates)}')

    if args.dry_run:
        print(f'[dry-run] would upload to {manifest_path} and register a manifest analysis\n')
        print(content)
        return

    target = to_path(manifest_path)
    if target.exists() and not args.force:
        raise ValueError(f'{manifest_path} already exists; re-run with --force to overwrite')

    target.write_text(content)
    print(f'Uploaded {manifest_path}')

    analysis_id = register_manifest_analysis(args.project, args.cohort, manifest_path)
    print(f'Registered manifest analysis {analysis_id} against cohort {args.cohort} in {args.project}')


if __name__ == '__main__':
    main()
