# Project Context: popgen-genotyping

## Project Overview
This repository contains the genotyping microarray pipeline used by Population Genomics. It is built on `cpg-flow` and orchestrates the conversion of raw Illumina GTC files into standardized cohort-level datasets (PLINK2 and BCF).

## Core Pipeline Architecture
The pipeline follows a specific sequence to ensure stability and data integrity:
1.  **GtcToBcfs**: Converts raw GTC files to "Heavy" (intensities) and "Light" (GT/GQ only) BCFs using `bcftools +gtc2vcf`.
2.  **BafRegress**: Estimates sample contamination from the Light BCFs.
3.  **CohortBcfToPlink**: Converts Light BCFs to individual PLINK 1.9 binary files.
    *   *Note*: Uses PLINK2 internally for initial conversion to leverage `--max-alleles 2` and `--split-par hg38` for clean biallelic filtering.
4.  **MergeCohortPlink**: Merges individual samples into a single PLINK 1.9 cohort using `plink --merge-list`. This provides a stable merging path that avoids limitations in current PLINK2 alpha builds.
5.  **ExportCohortDatasets**: Converts the merged PLINK 1.9 cohort into modern formats for long-term storage:
    *   **PLINK2 (PGEN/PVAR/PSAM)**
    *   **BCF (with index)**

## Standardized Tooling & Docker
The pipeline uses official CPG common images. No local tool builds are required.

### Official Images
- **BCFtools**: `australia-southeast1-docker.pkg.dev/cpg-common/images/bcftools:1.23-1`
- **PLINK (Multi-version)**: `australia-southeast1-docker.pkg.dev/cpg-common/images/plink:1.9-20250819-PLINK-2.0-20260228-1`
    *   This image contains both `plink1.9` and `plink2` binaries, routed via a wrapper.

### Usage Instructions
1.  **Start Docker Desktop**: Mandatory before running any local tests.
    ```bash
    open -a Docker
    ```
2.  **Set Environment Variables**:
    ```bash
    export PATH="/Applications/Docker.app/Contents/Resources/bin:$PATH"
    ```

## Local Development & Testing
### Data Generation
You can generate realistic synthetic GTC files for testing using:
```bash
python3 test/scripts/generate_synthetic_gtc.py <output_path> --num 100000 --name GDA-8v1-0_D2.bpm
```

### Full Pipeline Validation
The entire workflow (GTC -> BCF -> PLINK 1.9 -> PLINK2/BCF) can be verified locally using the reproduction script:
```bash
python3 test/scripts/reproduce_full_pipeline.py
```

## Engineering Standards
- **Linting**: Zero `ruff` errors are mandatory. Always run `ruff check src/ test/` before committing.
- **Documentation**: All functions, classes, and stages must have PEP 257 standardized docstrings (Description, Args, Returns, Types).
- **Strict Typing**: All variables and function signatures must have explicit type annotations.
- **Job Factory**: Use `popgen_genotyping.utils.register_job` to initialize all Hail Batch jobs to ensure consistent configuration (CPU, Memory, Image) from the project config.

## Operational Protocols (Optional)
**Note**: The following protocols are managed via the agent's global memory. Check global memory first for specific configurations (e.g., Sprint IDs, Board IDs).

### 1. Jira Integration (Phase 0)
- All work must be linked to a Jira ticket.
- Every task begins with documenting details in `ticket.md` and obtaining approval before posting.
- Tickets must be associated with the current active sprint (ID found in global memory).

### 2. Source Control
- Explicit written permission is required from the user before performing any `git commit` or `git push`.
- Never stage or commit changes automatically.
