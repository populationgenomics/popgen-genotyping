"""Merge PLINK2 QC data files and bafregress results into a single CSV summary."""

import argparse
import csv
import sys

# KING kinship coefficient thresholds for relatedness degree classification.
# Ref: Manichaikul et al. (2010), doi:10.1093/bioinformatics/btq559
KINSHIP_MZ: float = 0.354
KINSHIP_1ST: float = 0.177
KINSHIP_2ND: float = 0.0884
KINSHIP_3RD: float = 0.0442


def read_whitespace_file(filepath: str) -> list[dict[str, str]]:
    """Read a whitespace-delimited file into a list of row dicts.

    Args:
        filepath: Path to the whitespace-delimited file.

    Returns:
        List of dictionaries, one per row, keyed by column name.
    """
    rows: list[dict[str, str]] = []
    with open(filepath) as f:
        header: list[str] | None = None
        for line in f:
            fields = line.strip().split()
            if not fields:
                continue
            if header is None:
                header = [col.lstrip('#') for col in fields]
            else:
                rows.append(dict(zip(header, fields, strict=False)))
    return rows


def index_by(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    """Build a lookup dict keyed by a single column value.

    Args:
        rows: List of row dicts.
        key: Column name to index by.

    Returns:
        Dict mapping key values to row dicts.
    """
    return {row[key]: row for row in rows}


def left_merge(
    left: list[dict[str, str]],
    right: list[dict[str, str]],
    key: str,
    suffix: str = '',
) -> list[dict[str, str]]:
    """Left-join two lists of row dicts on a single key column.

    Args:
        left: Left-side rows.
        right: Right-side rows.
        key: Column name to join on.
        suffix: Suffix to append to duplicate column names from right.

    Returns:
        Merged list of row dicts.
    """
    right_idx = index_by(right, key)
    left_cols: set[str] = {col for row in left for col in row}
    merged: list[dict[str, str]] = []
    for row in left:
        new_row = dict(row)
        match = right_idx.get(row[key], {})
        for col, val in match.items():
            if col == key:
                continue
            out_col = f'{col}{suffix}' if col in left_cols else col
            new_row[out_col] = val
        merged.append(new_row)
    return merged


def get_degree(kinship: float) -> str | None:
    """Categorise a kinship coefficient into a relatedness degree.

    Args:
        kinship: KING kinship coefficient.

    Returns:
        Degree string or None if below 3rd degree threshold.
    """
    if kinship >= KINSHIP_MZ:
        return 'RELATED_MZ'
    if kinship >= KINSHIP_1ST:
        return 'RELATED_1ST'
    if kinship >= KINSHIP_2ND:
        return 'RELATED_2ND'
    if kinship >= KINSHIP_3RD:
        return 'RELATED_3RD'
    return None


DEGREE_COLS: list[str] = ['RELATED_MZ', 'RELATED_1ST', 'RELATED_2ND', 'RELATED_3RD']


def process_kinship(kin0_path: str) -> dict[str, dict[str, str]]:
    """Process a PLINK2 .kin0 file into per-sample relatedness columns.

    Each degree column contains semicolon-separated ``REL_ID:KINSHIP:IBS0`` strings.

    Args:
        kin0_path: Path to the .kin0 file.

    Returns:
        Dict keyed by IID, with values being dicts of degree columns.
    """
    rows = read_whitespace_file(kin0_path)
    if not rows:
        return {}

    required = {'IID1', 'IID2', 'KINSHIP', 'IBS0'}
    if not required.issubset(rows[0].keys()):
        print(
            'Warning: Expected IID1, IID2, KINSHIP, and IBS0 columns in .kin0 file but column(s) not found.',
        )
        return {}

    # Build two-way pairs filtered to >= 3rd degree
    pairs: list[tuple[str, str, float, float]] = []
    for row in rows:
        kinship = float(row['KINSHIP'])
        if kinship < KINSHIP_3RD:
            continue
        ibs0 = float(row['IBS0'])
        pairs.append((row['IID1'], row['IID2'], kinship, ibs0))
        pairs.append((row['IID2'], row['IID1'], kinship, ibs0))

    if not pairs:
        return {}

    # Group by (IID, degree)
    grouped: dict[str, dict[str, list[str]]] = {}
    for iid, rel_id, kinship, ibs0 in pairs:
        degree = get_degree(kinship)
        if degree is None:
            continue
        grouped.setdefault(iid, {}).setdefault(degree, []).append(
            f'{rel_id}:{kinship:.4f}:{ibs0:.4f}',
        )

    # Flatten to per-sample dicts
    result: dict[str, dict[str, str]] = {}
    for iid, degrees in grouped.items():
        result[iid] = {deg: ';'.join(degrees.get(deg, [])) for deg in DEGREE_COLS}
    return result


def process_bafregress(paths: list[str]) -> list[dict[str, str]]:
    """Read and concatenate bafregress output files.

    Args:
        paths: List of paths to bafregress output files.

    Returns:
        Combined list of row dicts with sample_id renamed to IID.
    """
    all_rows: list[dict[str, str]] = []
    for path in paths:
        try:
            rows = read_whitespace_file(path)
            for row in rows:
                if 'sample_id' in row:
                    row['IID'] = row.pop('sample_id')
                all_rows.append(row)
        except (OSError, ValueError) as e:
            print(f'Warning: Could not read bafregress file {path}. Error: {e}')
    return all_rows


def write_csv(rows: list[dict[str, str]], output_path: str) -> None:
    """Write a list of row dicts to a CSV file.

    Args:
        rows: List of row dicts to write.
        output_path: Path to the output CSV file.
    """
    if not rows:
        print('Warning: No data to write.')
        return

    # Preserve column order from first row, then append any extras
    fieldnames: list[str] = list(rows[0].keys())
    seen: set[str] = set(fieldnames)
    for row in rows[1:]:
        for key in row:
            if key not in seen:
                fieldnames.append(key)
                seen.add(key)

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    """Parse arguments and run the QC merge pipeline."""
    parser = argparse.ArgumentParser(
        description='Merge PLINK2 QC files and bafregress results into a CSV.',
    )
    parser.add_argument('--sexcheck', required=True, help='Path to .sexcheck file')
    parser.add_argument('--het', required=True, help='Path to .het file')
    parser.add_argument('--smiss', required=True, help='Path to .smiss file')
    parser.add_argument('--kin0', required=True, help='Path to .kin0 file')
    parser.add_argument('--output', required=True, help='Output CSV path')
    parser.add_argument(
        '--bafregress',
        nargs='*',
        default=[],
        help='Paths to bafregress output files',
    )
    args = parser.parse_args()

    # Read PLINK2 QC files
    try:
        smiss = read_whitespace_file(args.smiss)
        het = read_whitespace_file(args.het)
        sexcheck = read_whitespace_file(args.sexcheck)
    except (OSError, ValueError) as e:
        print(f'Error reading QC files: {e}')
        sys.exit(1)

    # Merge smiss, het, sexcheck on IID
    merged = left_merge(smiss, het, key='IID', suffix='_het')
    merged = left_merge(merged, sexcheck, key='IID', suffix='_sex')

    # Process kinship data
    try:
        kin_data = process_kinship(args.kin0)
        if kin_data:
            for row in merged:
                kin_row = kin_data.get(row.get('IID', ''), {})
                for col in DEGREE_COLS:
                    row[col] = kin_row.get(col, '')
        else:
            for row in merged:
                for col in DEGREE_COLS:
                    row[col] = ''
    except (OSError, ValueError, KeyError) as e:
        print(f'Warning: Could not process kinship file: {e}')

    # Process bafregress files
    if args.bafregress:
        baf_rows = process_bafregress(args.bafregress)
        if baf_rows:
            merged = left_merge(merged, baf_rows, key='IID', suffix='_baf')

    write_csv(merged, args.output)


if __name__ == '__main__':
    main()
