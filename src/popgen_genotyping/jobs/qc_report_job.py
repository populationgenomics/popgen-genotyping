"""
Job to merge the QC data files into a csv.
"""

from hailtop.batch.job import BashJob
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from popgen_genotyping.utils import register_job


def run_qc_report(
    plink_qc_prefix: str,
    bafregress_paths: list[str],
    output_path: str,
    job_name: str = 'qc_report',
) -> 'BashJob':
    """
    Merge the QC files into a single csv QC summary file.

    Args:
        cohort_qc_paths (str): Cloud paths to Plink QC files prefix.
        bafregress_paths (list[str]): Cloud paths to bafregress files for merging.
        output_path (str): Cloud path to output QC summary csv.
        job_name (str): Name for the Hail Batch Job.

    Returns:
        Job: A Hail Batch job object.
    """
    b = get_batch()
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'merge_qc'],
        image=config_retrieve(['workflow', 'plink_image']),  # TODO maybe swap image to something generic
        default_cpu=2,
        default_storage='20G',
    )

    # Read in PLINK2 QC files
    sexcheck_file = b.read_input(f'{plink_qc_prefix}.sexcheck')
    het_file = b.read_input(f'{plink_qc_prefix}.het')
    smiss_file = b.read_input(f'{plink_qc_prefix}.smiss')
    kin0_file = b.read_input(f'{plink_qc_prefix}.kin0')

    # Read in all bafregress files, then convert to
    # space-separated string of files as input to script.
    bafregress_files = [b.read_input(path) for path in bafregress_paths]
    bafregress_bash_args = ' '.join([str(f) for f in bafregress_files])

    # Define output file
    j.declare_resource_group(
        output_csv={'qc_report': '{root}.csv'}  # TODO double check syntax here
    )

    j.command(
        f"""
        set -ex
        cat <<EOF > merge_qc.py
import pandas as pd
import sys

try:
    smiss = pd.read_csv("{smiss_file}", sep=r'\s+')
    het = pd.read_csv("{het_file}", sep=r'\s+')
    sexcheck = pd.read_csv("{sexcheck_file}", sep=r'\s+')
except Exception as e:
    print(f"Error reading QC files: {{e}}")
    sys.exit(1)

# Strip the '#' from PLINK's #FID or #IID headers
for df in [smiss, het, sexcheck]:
    df.columns = df.columns.str.lstrip('#')

# Identify merge keys safely
merge_cols = [col for col in ['FID', 'IID'] if col in smiss.columns]

# Merge the dataframes sequentially
merged_df = smiss.merge(het, on=merge_cols, suffixes=('', '_het'))
merged_df = merged_df.merge(sexcheck, on=merge_cols, suffixes=('', '_sex'))

# Process PLINK2 Kinship data
expected_kin_cols = ['RELATED_MZ', 'RELATED_1ST', 'RELATED_2ND', 'RELATED_3RD']

try:
    kin_df = pd.read_csv("{kin0_file}", sep=r'\s+')
    kin_df.columns = kin_df.columns.str.lstrip('#')

    if all(c in kin_df.columns for c in ['IID1', 'IID2', 'KINSHIP', 'IBS0']):
        # Filter out anything less than 3rd degree (0.0442) immediately
        kin_df = kin_df[kin_df['KINSHIP'] >= 0.0442]

        if not kin_df.empty:
            # Create the two-way mapping for grouping, including IBS0
            kin_1 = kin_df[['IID1', 'IID2', 'KINSHIP', 'IBS0']].copy()
            kin_1.columns = ['IID', 'REL_ID', 'KINSHIP', 'IBS0']

            kin_2 = kin_df[['IID2', 'IID1', 'KINSHIP', 'IBS0']].copy()
            kin_2.columns = ['IID', 'REL_ID', 'KINSHIP', 'IBS0']

            kin_stacked = pd.concat([kin_1, kin_2], ignore_index=True)

            # Format the string payload: ID:KINSHIP:IBS0
            kin_stacked['REL_STR'] = (
                kin_stacked['REL_ID'].astype(str) + ':' +
                kin_stacked['KINSHIP'].round(4).astype(str) + ':' +
                kin_stacked['IBS0'].round(4).astype(str)
            )

            # Categorize degrees with mutually exclusive >= boundaries
            def get_degree(k):
                if k >= 0.354: return 'RELATED_MZ'
                elif k >= 0.177: return 'RELATED_1ST'
                elif k >= 0.0884: return 'RELATED_2ND'
                elif k >= 0.0442: return 'RELATED_3RD'
                return None

            kin_stacked['DEGREE'] = kin_stacked['KINSHIP'].apply(get_degree)

            # Group by sample and degree category, join strings with semicolons
            kin_grouped = kin_stacked.groupby(['IID', 'DEGREE'])['REL_STR'].apply(lambda x: ';'.join(x)).reset_index()

            # Pivot the dataframe so degrees become distinct columns
            kin_pivoted = kin_grouped.pivot(index='IID', columns='DEGREE', values='REL_STR').reset_index()

            # Guarantee all degree columns exist
            for col in expected_kin_cols:
                if col not in kin_pivoted.columns:
                    kin_pivoted[col] = pd.NA

            # Merge into the main dataframe
            merged_cols = ['IID'] + expected_kin_cols
            merged_df = merged_df.merge(kin_pivoted[merged_cols], on='IID', how='left')
        else:
            # No relationships found in the whole cohort
            for col in expected_kin_cols:
                merged_df[col] = pd.NA
    else:
        print("Warning: Expected IID1, IID2, KINSHIP, and IBS0 columns in .kin0 file but didn't find them.")
except Exception as e:
    print(f"Warning: Could not process kinship file: {{e}}")# Concatenate and merge the bafregress files

# Process bafregress
bafregress_paths = sys.argv[1:]
if bafregress_paths:
    bafregress_dfs = []
    for path in bafregress_paths:
        try:
            df = pd.read_csv(path, sep=r'\s+')
            bafregress_dfs.append(df)
        except Exception as e:
            print(f'Warning: Could not read bafregress file {{path}}. Error: {{e}}')

    if bafregress_dfs:
        bafregress_combined = pd.concat(bafregress_dfs, ignore_index=True)
        bafregress_combined.rename(columns={{'sample_id': 'IID'}}, inplace=True)

        # Merge with the plink2 QC data
        merged_df = merged_df.merge(bafregress_combined, on='IID', how='left', suffixes=('', '_baf'))

# Write directly to the Hail Batch output resource path
merged_df.to_csv("{j.output_csv.qc_report}", index=False)
EOF

        python3 merge_qc.py {bafregress_bash_args}
        """
    )

    b.write_output(j.output_csv.qc_report, output_path)

    return j
