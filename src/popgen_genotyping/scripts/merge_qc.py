"""Merge PLINK2 QC data files, KING --ibdseg relatedness, and bafregress results into a CSV summary.

Also folds in the set of samples excluded from KING --ibdseg for contamination or high
missingness (flagged RELATEDNESS_EXCLUDED with the exclusion reason), so they keep their
QC columns without appearing to have relatives.
"""

import argparse
import sys

import pandas as pd

DEGREE_COLS: list[str] = ['RELATED_MZ', 'RELATED_1ST', 'RELATED_2ND', 'RELATED_3RD']

# Map KING --ibdseg InfType values to RELATED_* degree columns.
# Ref: Manichaikul et al. (2010), doi:10.1093/bioinformatics/btq559;
# Chen et al. (2024) KING --ibdseg. KING emits InfType in {Dup/MZ, PO, FS, 2nd, 3rd}
# under --degree 3; PO and FS both fall in the 1st-degree bin.
INFTYPE_TO_DEGREE: dict[str, str] = {
    'Dup/MZ': 'RELATED_MZ',
    'PO': 'RELATED_1ST',
    'FS': 'RELATED_1ST',
    '2nd': 'RELATED_2ND',
    '3rd': 'RELATED_3RD',
}


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


def process_seg(seg_path: str) -> pd.DataFrame:
    """Process a KING ``--ibdseg`` ``.seg`` file into per-sample relatedness columns.

    Each degree column contains semicolon-separated ``REL_ID:KINSHIP:INFTYPE``
    strings. Kinship is the IBD-based KING kinship coefficient
    ``IBD1Seg/4 + IBD2Seg/2`` (equivalent to ``PropIBD/2``). InfType is KING's
    relationship inference (``PO``, ``FS``, ``2nd``, ``3rd``, ``Dup/MZ``); any
    InfType not in that set raises ``ValueError``.

    Args:
        seg_path: Path to the autosomal ``.seg`` file emitted by KING
            ``--ibdseg``.

    Returns:
        DataFrame with IID and RELATED_MZ/1ST/2ND/3RD columns, or an empty
        DataFrame with those columns if no relationships are found.

    Raises:
        ValueError: If the ``.seg`` file is missing any of the required
            columns or contains an unrecognised ``InfType``.
    """
    empty = pd.DataFrame(columns=['IID', *DEGREE_COLS])

    seg_df = read_qc_file(seg_path)
    if seg_df.empty:
        return empty

    required = {'ID1', 'ID2', 'IBD1Seg', 'IBD2Seg', 'InfType'}
    missing = required - set(seg_df.columns)
    if missing:
        raise ValueError(
            f'.seg file {seg_path} is missing required column(s): {sorted(missing)}',
        )

    unknown = set(seg_df['InfType']) - set(INFTYPE_TO_DEGREE)
    if unknown:
        raise ValueError(
            f'.seg file {seg_path} contains unrecognised InfType value(s): {sorted(unknown)}. '
            f'Expected one of {sorted(INFTYPE_TO_DEGREE)}.',
        )

    seg_df = seg_df.assign(
        KINSHIP=(seg_df['IBD1Seg'].astype(float) / 4.0 + seg_df['IBD2Seg'].astype(float) / 2.0).round(4),
        DEGREE=seg_df['InfType'].map(INFTYPE_TO_DEGREE),
    )
    if seg_df.empty:
        return empty

    rel_1 = seg_df[['ID1', 'ID2', 'KINSHIP', 'InfType', 'DEGREE']].copy()
    rel_1.columns = pd.Index(['IID', 'REL_ID', 'KINSHIP', 'INFTYPE', 'DEGREE'])

    rel_2 = seg_df[['ID2', 'ID1', 'KINSHIP', 'InfType', 'DEGREE']].copy()
    rel_2.columns = pd.Index(['IID', 'REL_ID', 'KINSHIP', 'INFTYPE', 'DEGREE'])

    rel_stacked = pd.concat([rel_1, rel_2], ignore_index=True)

    # Format the string payload: REL_ID:KINSHIP:INFTYPE
    rel_stacked['REL_STR'] = (
        rel_stacked['REL_ID'].astype(str)
        + ':'
        + rel_stacked['KINSHIP'].astype(str)
        + ':'
        + rel_stacked['INFTYPE'].astype(str)
    )

    rel_grouped = rel_stacked.groupby(['IID', 'DEGREE'])['REL_STR'].apply(';'.join).reset_index()
    rel_pivoted = rel_grouped.pivot_table(
        index='IID',
        columns='DEGREE',
        values='REL_STR',
        aggfunc='first',
    ).reset_index()

    for col in DEGREE_COLS:
        if col not in rel_pivoted.columns:
            rel_pivoted[col] = pd.NA

    return rel_pivoted[['IID', *DEGREE_COLS]]


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


def process_relatedness_excluded(path: str, report_iids: pd.Series) -> pd.DataFrame:
    """Flag every report IID excluded from KING ``--ibdseg`` with its exclusion reason.

    Args:
        path: Path to the excluded-samples TSV emitted by ``relatedness_exclusions.py``
            (columns ``IID``, ``REASON``, ``BAF_REGRESS``, ``F_MISS``; header-only
            when nothing was excluded).
        report_iids: The full set of IIDs in the QC report.

    Returns:
        DataFrame with ``IID`` and ``RELATEDNESS_EXCLUDED`` for every IID in
        ``report_iids``: the ``REASON`` string (``contamination``,
        ``missingness``, or ``contamination;missingness``) for excluded
        samples, ``NA`` otherwise.

    Raises:
        ValueError: If an excluded IID is not present in ``report_iids``.
    """
    excluded_df = read_qc_file(path)

    missing = sorted(set(excluded_df['IID']) - set(report_iids))
    if missing:
        raise ValueError(f'Relatedness-excluded sample(s) not found in QC report: {missing}')

    reason_by_iid = excluded_df.set_index('IID')['REASON']
    return pd.DataFrame(
        {
            'IID': report_iids,
            'RELATEDNESS_EXCLUDED': report_iids.map(reason_by_iid),
        },
    )


def main() -> None:
    """Parse arguments and run the QC merge pipeline."""
    parser = argparse.ArgumentParser(
        description=(
            'Merge PLINK2 QC files, KING --ibdseg relatedness (plus its contamination/'
            'missingness-based exclusion list), and bafregress into a CSV.'
        ),
    )
    parser.add_argument('--sexcheck', required=True, help='Path to .sexcheck file')
    parser.add_argument('--het', required=True, help='Path to .het file')
    parser.add_argument('--smiss', required=True, help='Path to .smiss file')
    parser.add_argument('--seg', required=True, help='Path to KING --ibdseg autosomal .seg file')
    parser.add_argument(
        '--relatedness-excluded',
        required=True,
        help='Path to the TSV of samples excluded from KING --ibdseg for contamination or missingness',
    )
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

    # Process KING IBD-segment relatedness. Errors from process_seg propagate
    # so the Batch job fails fast on a malformed .seg rather than emitting a
    # silently-incomplete QC report.
    seg_pivoted = process_seg(args.seg)
    merged = merged.merge(seg_pivoted, on='IID', how='left')

    # Flag samples excluded from KING --ibdseg for contamination or missingness. KING never
    # saw them, so they must not carry any relatedness — a non-null RELATED_* value here
    # would mean the exclusion list and the .seg file disagree about which samples went
    # into KING.
    relatedness_excluded = process_relatedness_excluded(args.relatedness_excluded, merged['IID'])
    merged = merged.merge(relatedness_excluded, on='IID', how='left')

    excluded_rows = merged[merged['RELATEDNESS_EXCLUDED'].notna()]
    has_relatives = excluded_rows[DEGREE_COLS].notna().any(axis=1)
    if has_relatives.any():
        bad_iids = sorted(excluded_rows.loc[has_relatives, 'IID'])
        raise ValueError(f'Relatedness-excluded sample(s) unexpectedly have KING relatives: {bad_iids}')

    # Place RELATEDNESS_EXCLUDED immediately after RELATED_3RD.
    cols = [c for c in merged.columns if c != 'RELATEDNESS_EXCLUDED']
    cols.insert(cols.index('RELATED_3RD') + 1, 'RELATEDNESS_EXCLUDED')
    merged = merged[cols]

    # Process bafregress files
    if args.bafregress:
        baf_df = process_bafregress(args.bafregress)
        if not baf_df.empty:
            merged = merged.merge(baf_df, on='IID', how='left', suffixes=('', '_baf'))

    merged.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
