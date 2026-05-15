"""Merge PLINK2 QC data files, KING ``--ibdseg`` results and bafregress outputs into a single CSV summary."""

import argparse
import sys

import pandas as pd

# Map KING --ibdseg `InfType` labels to the report's relatedness-degree columns.
# Ref: Manichaikul et al. (2010) doi:10.1093/bioinformatics/btq559;
# Chen et al. (2024) for the IBD-segment-based inference used by `--ibdseg`.
INFTYPE_TO_DEGREE: dict[str, str] = {
    'Dup/MZ': 'RELATED_MZ',
    'PO': 'RELATED_1ST',
    'FS': 'RELATED_1ST',
    '2nd': 'RELATED_2ND',
    '3rd': 'RELATED_3RD',
}

DEGREE_COLS: list[str] = ['RELATED_MZ', 'RELATED_1ST', 'RELATED_2ND', 'RELATED_3RD']

IBDSEG_REQUIRED_COLS: set[str] = {'ID1', 'ID2', 'IBD1Seg', 'IBD2Seg', 'PropIBD', 'InfType'}


def read_qc_file(filepath: str) -> pd.DataFrame:
    """Read a whitespace-delimited PLINK2/KING QC file.

    Strips leading ``#`` from column headers.

    Args:
        filepath: Path to the whitespace-delimited file.

    Returns:
        DataFrame with cleaned column names.
    """
    df = pd.read_csv(filepath, sep=r'\s+', engine='python')
    df.columns = df.columns.str.lstrip('#')
    return df


def process_ibdseg(seg_path: str) -> pd.DataFrame:
    """Process a KING ``--ibdseg`` .seg file into per-sample relatedness columns.

    Each degree column contains semicolon-separated
    ``REL_ID:IBD0:IBD1:IBD2:PropIBD`` strings, where IBD0 is derived as
    ``max(0.0, 1.0 - IBD1Seg - IBD2Seg)``.

    Args:
        seg_path: Path to the KING ``--ibdseg`` .seg file.

    Returns:
        DataFrame with IID and RELATED_MZ/1ST/2ND/3RD columns, or an empty
        DataFrame with those columns if no related pairs are found.
    """
    empty = pd.DataFrame(columns=['IID', *DEGREE_COLS])

    seg_df = read_qc_file(seg_path)
    if seg_df.empty:
        return empty

    if not IBDSEG_REQUIRED_COLS.issubset(seg_df.columns):
        print(
            f'Warning: KING --ibdseg .seg file missing expected columns; '
            f'have {sorted(seg_df.columns)}, need {sorted(IBDSEG_REQUIRED_COLS)}.',
        )
        return empty

    # Restrict to relationships we map to a degree column (drops UN, 4th, etc.).
    seg_df = seg_df[seg_df['InfType'].isin(INFTYPE_TO_DEGREE)].copy()
    if seg_df.empty:
        return empty

    seg_df['DEGREE'] = seg_df['InfType'].map(INFTYPE_TO_DEGREE)
    seg_df['IBD0Seg'] = (1.0 - seg_df['IBD1Seg'] - seg_df['IBD2Seg']).clip(lower=0.0)

    # Two-way stack so each individual appears as the focal IID exactly once per pair.
    seg_a = seg_df[['ID1', 'ID2', 'IBD0Seg', 'IBD1Seg', 'IBD2Seg', 'PropIBD', 'DEGREE']].copy()
    seg_a.columns = pd.Index(['IID', 'REL_ID', 'IBD0', 'IBD1', 'IBD2', 'PropIBD', 'DEGREE'])

    seg_b = seg_df[['ID2', 'ID1', 'IBD0Seg', 'IBD1Seg', 'IBD2Seg', 'PropIBD', 'DEGREE']].copy()
    seg_b.columns = pd.Index(['IID', 'REL_ID', 'IBD0', 'IBD1', 'IBD2', 'PropIBD', 'DEGREE'])

    stacked = pd.concat([seg_a, seg_b], ignore_index=True)

    # Format the payload: REL_ID:IBD0:IBD1:IBD2:PropIBD
    stacked['REL_STR'] = (
        stacked['REL_ID'].astype(str)
        + ':'
        + stacked['IBD0'].round(4).astype(str)
        + ':'
        + stacked['IBD1'].round(4).astype(str)
        + ':'
        + stacked['IBD2'].round(4).astype(str)
        + ':'
        + stacked['PropIBD'].round(4).astype(str)
    )

    grouped = stacked.groupby(['IID', 'DEGREE'])['REL_STR'].apply(';'.join).reset_index()
    pivoted = grouped.pivot_table(
        index='IID',
        columns='DEGREE',
        values='REL_STR',
        aggfunc='first',
    ).reset_index()

    for col in DEGREE_COLS:
        if col not in pivoted.columns:
            pivoted[col] = pd.NA

    return pivoted[['IID', *DEGREE_COLS]]


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
        description='Merge PLINK2 QC files, KING --ibdseg and bafregress results into a CSV.',
    )
    parser.add_argument('--sexcheck', required=True, help='Path to .sexcheck file')
    parser.add_argument('--het', required=True, help='Path to .het file')
    parser.add_argument('--smiss', required=True, help='Path to .smiss file')
    parser.add_argument('--ibdseg', required=True, help='Path to KING --ibdseg .seg file')
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

    # Process KING --ibdseg relatedness data
    try:
        ibdseg_pivoted = process_ibdseg(args.ibdseg)
        merged = merged.merge(ibdseg_pivoted, on='IID', how='left')
    except (OSError, ValueError, KeyError) as e:
        print(f'Warning: Could not process KING --ibdseg file: {e}')

    # Ensure degree columns exist even if ibdseg processing was skipped
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
