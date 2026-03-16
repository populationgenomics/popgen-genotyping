# Research: Update Rolling Aggregate to use PLINK2 format

## Objective
Update the `MergeCohortPlink` stage to handle previous rolling aggregates stored in PLINK2 format (`.pgen`, `.pvar`, `.psam`) instead of PLINK 1.9 format.

## Current State
- `resolve_rolling_aggregate` in `metamist_utils.py` expects and returns PLINK 1.9 paths (`bed`, `bim`, `fam`).
- `MergeCohortPlink` stage directly passes these paths to the `run_merge_plink` job, which expects PLINK 1.9 format.

## Proposed Changes
1.  **Update `metamist_utils.py`**:
    - Modify `resolve_rolling_aggregate` to query for `pgen`, `pvar`, and `psam` files from the previous analysis `outputs`.
    - The function will now return these PLINK2 paths.

2.  **Create a New Conversion Job**:
    - Create `src/popgen_genotyping/jobs/plink2_to_plink1_job.py`.
    - This will contain a function `run_plink2_to_plink1` that takes a PLINK2 fileset prefix and converts it to a PLINK 1.9 fileset using `plink2 --make-bed`.

3.  **Update `MergeCohortPlink` Stage**:
    - In `stages.py`, if a `previous_analysis_id` is provided, the `MergeCohortPlink` stage will first queue the new `run_plink2_to_plink1` conversion job.
    - The existing `run_merge_plink` job will be made dependent on this new job.
    - The `previous_aggregate_paths` passed to `run_merge_plink` will be the *output* of the new conversion job, ensuring it receives the PLINK 1.9 files it expects.
    - If no `previous_analysis_id` is provided, this conversion step will be skipped, and the merge job will run as it currently does.

This approach correctly handles the format conversion within the workflow, isolating the logic and maintaining the integrity of the existing merge job, while also adhering to `cpg-flow`'s stage dependency requirements.
