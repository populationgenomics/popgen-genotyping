"""Select samples to exclude from KING --ibdseg for contamination or high missingness.

A contaminated sample's genotype noise, or a sample with unusually high per-sample
missingness, generates thousands of spurious pairwise relationships under KING
`--ibdseg`, which in turn contaminate every clean sample's relatedness list. This
script reads BAFRegress output and PLINK2 `.smiss` per-sample missingness, selects
the samples that exceed either configured threshold, and writes both an audit TSV
(IID, REASON, BAF_REGRESS, F_MISS) and a plink2 `--remove` sample-ID list.
"""

import argparse
import sys

import pandas as pd

CONTAMINATION_REASON: str = 'contamination'
MISSINGNESS_REASON: str = 'missingness'


def read_bafregress(paths: list[str]) -> pd.DataFrame:
    """Read and concatenate BAFRegress output files.

    Each file is one plate cohort's whitespace-delimited BAFRegress output,
    with a header row and columns ``sample_id``, ``baf_regress``, ``Nhom``.

    Args:
        paths: Paths to BAFRegress output files.

    Returns:
        The concatenated DataFrame across all files.

    Raises:
        ValueError: If any file lacks a ``sample_id`` or ``baf_regress`` column.
    """
    dfs: list[pd.DataFrame] = []
    for path in paths:
        df = pd.read_csv(path, sep=r'\s+')
        missing = {'sample_id', 'baf_regress'} - set(df.columns)
        if missing:
            raise ValueError(f'BAFRegress file {path} is missing required column(s): {sorted(missing)}')
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)


def read_smiss(path: str) -> pd.DataFrame:
    """Read a PLINK2 `.smiss` per-sample missingness file.

    Strips leading ``#`` from column headers, mirroring ``merge_qc.read_qc_file``.

    Args:
        path: Path to the ``.smiss`` file.

    Returns:
        DataFrame with cleaned column names.

    Raises:
        ValueError: If the file lacks an ``IID`` or ``F_MISS`` column.
    """
    df = pd.read_csv(path, sep=r'\s+', engine='python')
    df.columns = df.columns.str.lstrip('#')

    missing = {'IID', 'F_MISS'} - set(df.columns)
    if missing:
        raise ValueError(f'.smiss file {path} is missing required column(s): {sorted(missing)}')

    return df


def select_excluded_samples(
    bafregress: pd.DataFrame,
    smiss: pd.DataFrame,
    contamination_max: float,
    fmiss_max: float,
) -> pd.DataFrame:
    """Select samples to exclude from KING `--ibdseg` for contamination or missingness.

    The ``.smiss`` sample set defines the required universe: every merged-cohort
    sample must have a BAFRegress estimate, since ``.smiss`` is computed on the
    merged pgen and BAFRegress is resolved for the same membership.

    Args:
        bafregress: DataFrame with ``sample_id`` and ``baf_regress`` columns, as
            returned by :func:`read_bafregress`.
        smiss: DataFrame with ``IID`` and ``F_MISS`` columns, as returned by
            :func:`read_smiss`.
        contamination_max: Samples with ``baf_regress`` strictly greater than
            this value are excluded for contamination.
        fmiss_max: Samples with ``F_MISS`` strictly greater than this value are
            excluded for missingness.

    Returns:
        DataFrame with columns ``IID``, ``REASON``, ``BAF_REGRESS``, ``F_MISS``,
        sorted by ``IID``. ``REASON`` is ``contamination``, ``missingness``, or
        ``contamination;missingness``. Empty (but correctly columned) when no
        sample exceeds either threshold.

    Raises:
        ValueError: If any ``.smiss`` IID has no BAFRegress estimate.
    """
    missing = sorted(set(smiss['IID']) - set(bafregress['sample_id']))
    if missing:
        raise ValueError(f'No BAFRegress estimate found for sequencing group(s): {missing}')

    merged = smiss[['IID', 'F_MISS']].merge(
        bafregress[['sample_id', 'baf_regress']].rename(columns={'sample_id': 'IID', 'baf_regress': 'BAF_REGRESS'}),
        on='IID',
        how='left',
    )

    contaminated = merged['BAF_REGRESS'] > contamination_max
    high_missingness = merged['F_MISS'] > fmiss_max
    excluded_mask = contaminated | high_missingness

    def _reason(is_contaminated: bool, is_high_missingness: bool) -> str:
        reasons = [
            reason
            for reason, flag in ((CONTAMINATION_REASON, is_contaminated), (MISSINGNESS_REASON, is_high_missingness))
            if flag
        ]
        return ';'.join(reasons)

    result = merged.loc[excluded_mask, ['IID', 'BAF_REGRESS', 'F_MISS']].copy()
    result['REASON'] = [
        _reason(c, m) for c, m in zip(contaminated[excluded_mask], high_missingness[excluded_mask], strict=True)
    ]
    result = result[['IID', 'REASON', 'BAF_REGRESS', 'F_MISS']]
    return result.sort_values('IID').reset_index(drop=True)


def main() -> None:
    """Parse arguments, select excluded samples, and write the TSV + remove-list."""
    parser = argparse.ArgumentParser(
        description=(
            'Select samples to exclude from KING --ibdseg for BAFRegress contamination or '
            'high per-sample missingness, and write an audit TSV plus a plink2 --remove list.'
        ),
    )
    parser.add_argument('--bafregress', nargs='+', required=True, help='Paths to BAFRegress output files')
    parser.add_argument('--smiss', required=True, help='Path to the merged-cohort .smiss file')
    parser.add_argument('--contamination-max', type=float, required=True, help='BAFRegress exclusion threshold')
    parser.add_argument('--fmiss-max', type=float, required=True, help='Per-sample F_MISS exclusion threshold')
    parser.add_argument('--output-tsv', required=True, help='Output path for the excluded-samples audit TSV')
    parser.add_argument('--output-remove-list', required=True, help='Output path for the plink2 --remove list')
    args = parser.parse_args()

    try:
        bafregress = read_bafregress(args.bafregress)
        smiss = read_smiss(args.smiss)
    except (OSError, ValueError) as e:
        print(f'Error reading input files: {e}')
        sys.exit(1)

    excluded = select_excluded_samples(bafregress, smiss, args.contamination_max, args.fmiss_max)

    excluded.to_csv(args.output_tsv, sep='\t', index=False)

    with open(args.output_remove_list, 'w') as f:
        f.write('#IID\n')
        for iid in excluded['IID']:
            f.write(f'{iid}\n')

    # Avoid pandas' `.str` accessor: an empty REASON column (nothing excluded) has no
    # inferrable string dtype and raises AttributeError.
    contamination_count = sum(CONTAMINATION_REASON in reason for reason in excluded['REASON'])
    missingness_count = sum(MISSINGNESS_REASON in reason for reason in excluded['REASON'])
    print(
        f'Excluded {len(excluded)} sample(s) from KING --ibdseg: '
        f'{contamination_count} for contamination, {missingness_count} for missingness.',
    )


if __name__ == '__main__':
    main()
