# Jira Ticket: Update Rolling Aggregate to PLINK2 format

## Ticket Details
- **Project:** POPGEN (Board 291)
- **Type:** Task
- **Title:** Update Rolling Aggregate to Handle PLINK2 Format
- **Description:** 
The `MergeCohortPlink` stage needs to be updated to handle previous rolling aggregates stored in PLINK2 format (`.pgen`, `.pvar`, `.psam`). This requires adding a new conversion step to the workflow.

Key changes:
1. Update `resolve_rolling_aggregate` to retrieve PLINK2 file paths from Metamist.
2. Create a new job (`plink2_to_plink1`) to convert a PLINK2 fileset to PLINK 1.9 format.
3. Integrate this new conversion job into the `MergeCohortPlink` stage, making it a dependency for the merge job.
- **Story Points:** 3
- **Epic:** POPGEN-942
- **Sprint:** 1453
- **Definition of Done:**
  - The pipeline successfully processes rolling aggregates stored in PLINK2 format.
  - The PLINK2 to PLINK 1.9 conversion is correctly integrated.
  - All code changes are tested and pass linting.
