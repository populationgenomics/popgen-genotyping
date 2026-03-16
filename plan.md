# Plan: Update Rolling Aggregate to PLINK2 format

## Objective
Implement the changes required for the `MergeCohortPlink` stage to process previous aggregates stored in PLINK2 format, with conditional logic to handle the conversion.

## Implementation Details
1.  **Modify `src/popgen_genotyping/metamist_utils.py`**:
    - Update the `resolve_rolling_aggregate` function to retrieve `pgen`, `pvar`, and `psam` paths. The sample parsing logic will also be updated to use the `.psam` file.
2.  **Create `src/popgen_genotyping/jobs/plink2_to_plink1_job.py`**:
    - Implement the `run_plink2_to_plink1` job function. This job will accept a PLINK2 prefix, run `plink2 --pfile <prefix> --make-bed --out <output>`, and return the resulting PLINK 1.9 `ResourceGroup`.
3.  **Modify `src/popgen_genotyping/stages.py`**:
    - In the `MergeCohortPlink` stage, add conditional logic:
        - If `previous_analysis_id` is present, queue the `run_plink2_to_plink1` job and make the `run_merge_plink` job dependent on it.
        - If not, proceed with the merge as normal.
    - The output of the conversion job will be passed to `run_merge_plink`.

## Testing Strategy
- Run existing unit tests to check for regressions.
- Add a new unit test for the `run_plink2_to_plink1` job.
- Manually review the conditional logic in the `MergeCohortPlink` stage.

## Definition of Done
- The `MergeCohortPlink` stage correctly handles both scenarios (with and without a previous PLINK2 aggregate).
- The PLINK2 to PLINK 1.9 conversion is successfully and conditionally integrated.
- All new and modified code passes linting checks.
- Existing tests pass, and new tests are added where appropriate.
