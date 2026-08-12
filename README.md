# Population Genomics Genotyping Pipeline

A `cpg-flow` pipeline to process genotyping microarray data, from raw Illumina GTC files to cohort-level PLINK and BCF datasets.

## Project Overview
This pipeline is designed to automate the conversion of dense, single-sample GTC files into analysis-ready, multi-sample datasets. It handles data conversion, quality control, and formatting, producing standard outputs suitable for downstream genetic analysis.

## Pipeline Architecture
The pipeline is composed of several sequential stages, orchestrated by `cpg-flow`. Each stage is responsible for a specific part of the data processing workflow.

![Pipeline DAG](pipeline_dag.png)

### Stages
- **GtcToBcfs**: Converts raw GTC files into two BCF formats: a "Heavy" BCF containing full intensity data and a "Light" BCF containing only genotype calls (GT) and quality scores (GQ).
- **BafRegress**: Estimates sample contamination by analyzing B-Allele Frequencies (BAF) against a population reference. If no reference is provided, it will estimate AF from the cohort. Output is written to durable, version-independent storage and registered as an `array_bafregress` Metamist analysis (one per plate cohort).
- **CohortBcfToPlink**: Converts the Light BCF into PLINK 1.9 binary format (`.bed`, `.bim`, `.fam`), preparing it for merging. Output is written to durable, version-independent storage and registered as an `array_cohort_bed` Metamist analysis (one per plate cohort). See [Per-plate outputs are immutable](#per-plate-outputs-are-immutable).
- **MergeCohortPlink**: Merges PLINK files from multiple cohorts into a single, unified dataset. This stage also supports a "rolling aggregate" workflow, where new samples are added to a previously generated aggregate. See [Rolling aggregate & two-phase run](#rolling-aggregate--two-phase-run).
- **ExportCohortDatasets**: Converts the merged PLINK 1.9 dataset into PLINK2 (`.pgen`) format for long-term storage and analysis, and `.bcf` format in temporary storage for ancestry analysis.
- **Plink2Qc**: Performs a standard suite of quality control checks on the final PLINK2 dataset, including sample/variant missingness, allele frequency, HWE, heterozygosity, and kinship.
- **KingIbdseg**: Runs KING `--ibdseg --degree 3` against the merged PLINK 1.9 dataset to call pairwise IBD segments. Emits autosomal `.seg` / `.segments.gz` and (when chrX SNPs are present) X-chr companions `X.seg` / `X.segments.gz`, plus the captured KING log. Outputs land in long-term storage and are registered as an `array_relatedness_ibdseg` Metamist analysis; folding the pairwise summary into the QC CSV is tracked as a follow-up.

### Per-plate outputs are immutable
`CohortBcfToPlink` and `BafRegress` write to durable, **version-independent** paths
(`<default_prefix>/<workflow>/<StageName>/<cohort_id>...`, with no `workflow.version`
segment). This decouples a plate's outputs from the pipeline version — a plate is processed
once and reused across runs, and phase 2 no longer has to run at the same `workflow.version`
before `tmp` is garbage-collected to find them.

Because there is no version segment, re-processing a plate **overwrites in place**. cpg-flow
skips a stage whose `expected_outputs` already exist, so the default is reuse (not silent
overwrite); an overwrite only happens if the outputs are deleted and regenerated — which
would change the bytes any prior aggregate was built from. **Treat per-plate outputs as
immutable once an aggregate references them.** To legitimately re-process a plate, prefer a
new cohort rather than overwriting an existing one.

### Rolling aggregate & two-phase run
Aggregate datasets are registered against a **super cohort** (previous aggregate SGs + new
plate SGs) so downstream consumers (e.g. the genomic atlas) can query array data by cohort.
cpg-flow cannot create or validate a cohort mid-run, so the pipeline runs in two phases:

1. **Phase 1** — run against the **new plate cohorts** (`input_cohorts=[new plates]`,
   `config_phase1.toml`). Produces and registers the per-plate `array_cohort_bed` and
   `array_bafregress` outputs.
2. **Create the super cohort** manually in Swagger (previous aggregate SGs ∪ new plate SGs).
3. **Phase 2** — run against the **super cohort** (`input_cohorts=[super]`,
   `config_phase2.toml`). Rolls the previous aggregate forward and merges only the new plates.

Every stage is a `CohortStage`, so nothing in the stage graph itself separates the phases: a
phase-2 submission would otherwise also run the per-plate stages on the super cohort, and a
phase-1 submission would run the aggregate stages once per plate. Each phase config therefore
pins `workflow.only_stages` to its phase's stages, and `run_workflow.py` refuses to submit when
`only_stages` is missing or mixes stages from both phases — a mismatched config fails at
submission, before any job is queued.

The previous aggregate is selected explicitly by **cohort ID** (`previous_aggregate_cohort_id`).
The new plates are **not listed** in phase-2 config — they are **derived**
(`NEW = super − previous aggregate`), and each new SG is resolved to its plate via the registered
`array_cohort_bed` analysis. The resolved plan (the contributing plate cohorts and their new-SG
counts) is printed to the driver log at submission: **confirm the plates match what you ran in
phase 1** before letting the run proceed. Use `scripts/list_aggregates.py` to pick the previous
aggregate cohort.

**Why derive the new plates instead of reusing the phase-1 plate list?** It might look simpler to
just carry the plate cohorts forward from phase 1, but deriving from the super cohort's membership
keeps that cohort the single source of truth and is robust to things a hand-carried list gets wrong:
- **Withdrawn / inactivated SGs** — an SG dropped since phase 1 is simply absent from the super
  cohort, so it is excluded automatically; a static plate list would silently re-include it.
- **Plates accumulated across several phase-1 runs** — phase 1 can run many times before a single
  phase-2 aggregation; you never have to collate every plate cohort ID from all those runs.
- **Custom / partial selection** — the super cohort can be any hand-picked SG set (a subset of a
  plate, or spanning plates); per-SG resolution handles this, whereas a plate list cannot.
- **No config drift** — the merged output cannot disagree with the cohort it is registered against;
  the post-`--keep` kept-sample-count assert plus the per-SG coverage check catch any plate
  never run through phase 1, failing loudly instead of producing a silently short aggregate.

The plan printout is what makes the derivation trustworthy: you get to eyeball the derived plate
set against your phase-1 runs rather than trusting it blind.

The final `plink --keep` trims the merged fileset to super-cohort membership (`merged ⊆ super`) and
asserts the kept-sample count equals the super cohort (`super ⊆ merged`) before the aggregate is
registered, so the released dataset cannot silently disagree with the cohort it registers against.

All phase-2 output filenames embed the super-cohort ID (e.g. `<cohort_id>_merged.bed`,
`<cohort_id>_<date>.pgen`). Every rolling aggregate gets a new super cohort, so successive
aggregates at the same `workflow.version` land on distinct paths — cpg-flow's skip-if-exists
can therefore never reuse a previous super cohort's merge for a new one.

## Prerequisites
Before running the pipeline, ensure you have the following tools installed and configured:
- **`cpg-flow`**: The core workflow management system.
- **`analysis-runner`**: The command-line tool used to launch `cpg-flow` pipelines in the cloud.
- **Docker**: Required for running the local reproduction scripts.

## Configuration
The pipeline is configured using a TOML file, one per phase: start from
`src/popgen_genotyping/config_phase1.toml` (per-plate processing) or
`src/popgen_genotyping/config_phase2.toml` (aggregation against the super cohort).

### Key Parameters
- `[workflow]`:
    - `dataset`: The analysis dataset for the output.
    - `input_cohorts`: A list of cohort IDs to include in the run — the new plate cohorts in
      phase 1, exactly the super cohort in phase 2.
    - `only_stages`: The stages belonging to the phase being run. Mandatory; a submission
      that omits it or mixes stages from both phases is rejected (see
      [Rolling aggregate & two-phase run](#rolling-aggregate--two-phase-run)).
    - `sequencing_type`: Must be set to `array`.
    - `driver_image`: The Docker image for the main `cpg-flow` driver.
    - `bcftools_image`, `plink_image`, `king_image`: Docker images for the respective tools.
- `[popgen_genotyping.references]`:
    - `fasta_ref_path`: Path to the human genome FASTA reference.
    - `bpm_manifest_path`: Path to the Illumina BPM manifest file.
    - `egt_cluster_path`: Path to the Illumina EGT cluster file.
    - `af_ref_path` (optional): Path to a VCF containing population allele frequencies for `BafRegress`.
- `[popgen_genotyping.merge_cohort_plink]`:
    - `previous_aggregate_cohort_id` (optional): The Metamist **cohort ID** of a previous aggregate to roll forward. Omit for a from-scratch (bootstrap) build. Use `scripts/list_aggregates.py` to list registered aggregate cohorts and pick one. See [Rolling aggregate & two-phase run](#rolling-aggregate--two-phase-run).

## Execution
To run the pipeline, use the `analysis-runner` command. You will need to specify the path to your configuration file, the output directory, and the script to execute.

```bash
analysis-runner
    --dataset <your-dataset>
    --output-dir <output-directory>
    --config config_phase1.toml  # or config_phase2.toml
    run_workflow.py
```

## Local Development & Testing
This repository includes scripts for local development and testing.

- **`test/scripts/reproduce_full_pipeline.py`**: A modular script that reproduces the entire pipeline workflow locally using Docker. This is useful for verifying changes and understanding the pipeline's behavior.
- **`test/scripts/reproduce_bafregress_production.py`**: A specialized test for `BafRegress` using a production-sized cohort (94 samples) and internal AF estimation.

To run the local tests, ensure Docker is running and execute the scripts directly. For example:
```bash
python3 test/scripts/reproduce_full_pipeline.py --samples 5 --snps 10000
```
