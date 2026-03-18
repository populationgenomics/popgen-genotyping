# Local Testing and Synthetic Data

This project uses a deterministic synthetic data system and a local Docker-based testing harness to validate the genotyping pipeline stages.

## Prerequisites
- **Docker Desktop**: Must be running.
- **BCFTools**: Installed locally for verification (optional).
- **Reference Files**: Local copies of BPM, EGT, and FNA must be in `test/local/reference/`.

## 1. Prepare Synthetic Data
Generate deterministic GTC files and a matching manifest CSV:
```bash
python3 scripts/prepare_test_data.py
```
This will create:
- 3 deterministic GTC files in `test/local/data/gtc/`.
- A matching manifest in `test/local/data/synthetic_manifest.csv`.

## 2. Run Local Stage Tests
Execute the local testing harness to convert GTCs to BCFs using the production Docker container:
```bash
python3 scripts/local_gtc_to_bcfs.py
```
This script mimics the `GtcToBcfs` pipeline stage, including:
- `gtc2vcf` conversion.
- `bcftools sort` and `bcftools norm`.
- Sample reheadering and indexing.
- Generation of both Heavy (full) and Light (stripped) BCFs.

## 3. Verify Results
Outputs are located in `test/local/output/`. You can verify them using local `bcftools`:
```bash
# Check sample IDs
bcftools query -l test/local/output/CPGSYN001.heavy.bcf

# Count variants
bcftools view -H test/local/output/CPGSYN001.heavy.bcf | wc -l
```

## Maintenance
- **Updating Loci Count**: If the BPM changes, update the `NUM_SNPS` constant in `scripts/prepare_test_data.py`.
- **Relocated Scripts**: Note that `generate_synthetic_manifest.py` has been moved to `test/scripts/` to separate test utilities from production code.
