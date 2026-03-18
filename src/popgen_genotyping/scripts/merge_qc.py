"""Merge PLINK2 QC data files and bafregress results into a single CSV summary."""

import argparse
import sys

import pandas as pd

# KING kinship coefficient thresholds for relatedness degree classification.
# Ref: Manichaikul et al. (2010), doi:10.1093/bioinformatics/btq559
KINSHIP_MZ: float = 0.354
KINSHIP_1ST: float = 0.177
KINSHIP_2ND: float = 0.0884
KINSHIP_3RD: float = 0.0442

DEGREE_COLS: list[str] = ['RELATED_MZ', 'RELATED_1ST', 'RELATED_2ND', 'RELATED_3RD']


def read_qc_file(filepath: str) -> pd.DataFrame:
    """Read a whitespace-delimited PLINK2 QC file.

    Strips leading ``#`` from column headers.

    Args:
        filepath: Path to the whitespace-delimited file.

    Returns:
        DataFrame with cleaned column names.
    """
    df = pd.read_csv(filepath, sep=r'\s+', engine='python')
    df.columns = df.columns.str.lstrip('#')
    return df


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


def process_kinship(kin0_path: str) -> pd.DataFrame:
    """Process a PLINK2 .kin0 file into per-sample relatedness columns.

    Each degree column contains semicolon-separated ``REL_ID:KINSHIP:IBS0`` strings.

    Args:
        kin0_path: Path to the .kin0 file.

    Returns:
        DataFrame with IID and RELATED_MZ/1ST/2ND/3RD columns, or an empty
        DataFrame with those columns if no relationships are found.
    """
    empty = pd.DataFrame(columns=['IID', *DEGREE_COLS])

    kin_df = read_qc_file(kin0_path)
    if kin_df.empty:
        return empty

    required = {'IID1', 'IID2', 'KINSHIP', 'IBS0'}
    if not required.issubset(kin_df.columns):
        print(
            'Warning: Expected IID1, IID2, KINSHIP, and IBS0 columns in .kin0 file but column(s) not found.',
        )
        return empty

    # Filter to >= 3rd degree
    kin_df = kin_df[kin_df['KINSHIP'] >= KINSHIP_3RD]
    if kin_df.empty:
        return empty

    # Create two-way mapping
    kin_1 = kin_df[['IID1', 'IID2', 'KINSHIP', 'IBS0']].copy()
    kin_1.columns = pd.Index(['IID', 'REL_ID', 'KINSHIP', 'IBS0'])

    kin_2 = kin_df[['IID2', 'IID1', 'KINSHIP', 'IBS0']].copy()
    kin_2.columns = pd.Index(['IID', 'REL_ID', 'KINSHIP', 'IBS0'])

    kin_stacked = pd.concat([kin_1, kin_2], ignore_index=True)

    # Format the string payload: REL_ID:KINSHIP:IBS0
    kin_stacked['REL_STR'] = (
        kin_stacked['REL_ID'].astype(str)
        + ':'
        + kin_stacked['KINSHIP'].round(4).astype(str)
        + ':'
        + kin_stacked['IBS0'].round(4).astype(str)
    )

    # Categorise degrees
    kin_stacked['DEGREE'] = kin_stacked['KINSHIP'].apply(get_degree)

    # Group and pivot
    kin_grouped = (
        kin_stacked.groupby(['IID', 'DEGREE'])['REL_STR']
        .apply(';'.join)
        .reset_index()
    )
    kin_pivoted = kin_grouped.pivot_table(
        index='IID', columns='DEGREE', values='REL_STR', aggfunc='first',
    ).reset_index()

    # Guarantee all degree columns exist
    for col in DEGREE_COLS:
        if col not in kin_pivoted.columns:
            kin_pivoted[col] = pd.NA

    return kin_pivoted[['IID', *DEGREE_COLS]]


def process_bafregress(paths: list[str]) -> pd.DataFrame:
    """Read and concatenate bafregress output files.

    Args:
        paths: List of paths to bafregress output files.

    Returns:
        Combined DataFrame with sample_id renamed to IID, or empty DataFrame.
    """
    dfs: list[pd.DataFrame] = []
    for path in paths:
        try:
            df = pd.read_csv(path, sep=r'\s+', engine='python')
            dfs.append(df)
        except (OSError, ValueError) as e:
            print(f'Warning: Could not read bafregress file {path}. Error: {e}')

    if not dfs:
        return pd.DataFrame()

    return pd.concat(dfs, ignore_index=True).rename(columns={'sample_id': 'IID'})


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
        smiss = read_qc_file(args.smiss)
        het = read_qc_file(args.het)
        sexcheck = read_qc_file(args.sexcheck)
    except (OSError, ValueError) as e:
        print(f'Error reading QC files: {e}')
        sys.exit(1)

    # Identify merge keys
    merge_cols = [col for col in ['FID', 'IID'] if col in smiss.columns]

    # Merge QC dataframes
    merged = smiss.merge(het, on=merge_cols, suffixes=('', '_het'))
    merged = merged.merge(sexcheck, on=merge_cols, suffixes=('', '_sex'))

    # Process kinship data
    try:
        kin_pivoted = process_kinship(args.kin0)
        merged = merged.merge(kin_pivoted, on='IID', how='left')
    except (OSError, ValueError, KeyError) as e:
        print(f'Warning: Could not process kinship file: {e}')

    # Ensure degree columns exist even if kinship processing was skipped
    for col in DEGREE_COLS:
        if col not in merged.columns:
            merged[col] = pd.NA

    # Process bafregress files
    if args.bafregress:
        baf_df = process_bafregress(args.bafregress)
        if not baf_df.empty:
            merged = merged.merge(baf_df, on='IID', how='left', suffixes=('', '_baf'))

    merged.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
