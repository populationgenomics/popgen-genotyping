"""
Job logic for extracting observed per-SNP cluster statistics from a cohort Heavy BCF.

Streams FORMAT/{GT, THETA, R} through a stdlib aggregator that groups by genotype
cluster (AA / AB / BB) and sex subset, emitting a bgzip-compressed, tabix-indexed
TSV with per-cluster counts and (mean, var) of THETA and R per variant.
"""

from __future__ import annotations

from importlib.resources import files
from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


# Trailing extension produced by the aggregator + tabix; used to derive the
# Hail Batch resource-group root from a full output path.
OUTPUT_SUFFIX: str = '.cluster_stats.tsv.gz'


def run_cohort_cluster_stats(
    bcf_path: str,
    output_stats_path: str,
    sex_mapping: dict[str, str],
    job_name: str = 'cohort_cluster_stats',
) -> BashJob:
    """
    Build the Hail Batch job that produces per-cohort cluster statistics.

    Args:
        bcf_path: Cloud path to the cohort Heavy BCF (FORMAT/THETA and R present).
        output_stats_path: Cloud path for the bgzipped TSV. Must end in
            `.cluster_stats.tsv.gz`; the `.tbi` is written alongside it.
        sex_mapping: SG ID -> PLINK sex code ('1'=male, '2'=female). Unknown
            codes are accepted but those samples are excluded from the
            sex-stratified subsets on chrX/Y/MT.
        job_name: Name shown in Hail Batch UI / logs.

    Returns:
        BashJob whose declared resource group writes both the `.tsv.gz`
        and `.tsv.gz.tbi` outputs.
    """
    if not output_stats_path.endswith(OUTPUT_SUFFIX):
        raise ValueError(f'output_stats_path must end with {OUTPUT_SUFFIX!r}, got {output_stats_path!r}')

    b = get_batch()
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'cohort_cluster_stats'],
        image=config_retrieve(['workflow', 'bcftools_image']),
        default_cpu=2,
        default_storage='30G',
    )

    # Stage Heavy BCF + index.
    bcf_file = b.read_input_group(bcf=bcf_path, csi=f'{bcf_path}.csi')

    # Stage the aggregator script. Bundled as package data; resolved via
    # importlib.resources so the job sees a real filesystem path.
    aggregator_script_path = str(files('popgen_genotyping.scripts').joinpath('aggregate_cluster_stats.py'))
    aggregator_script = b.read_input(aggregator_script_path)

    # Output resource group: tsv.gz + tabix index, both rooted at the same prefix.
    j.declare_resource_group(
        output={
            'tsv_gz': '{root}.cluster_stats.tsv.gz',
            'tbi': '{root}.cluster_stats.tsv.gz.tbi',
        },
    )

    # Inline sex mapping as a TSV literal. Cohorts run a few hundred samples
    # at most, so heredoc is comfortable; for larger cohorts, switch to a
    # staged file via to_path/read_input as merge_cohort_plink_job does.
    sex_lines = [f'{sg_id}\t{sex_code}' for sg_id, sex_code in sorted(sex_mapping.items())]
    sex_tsv_content = '\\n'.join(sex_lines)

    j.command(
        f"""
        set -euxo pipefail

        # 1. Sample order from BCF header (drives the per-sample column mapping).
        bcftools query -l {bcf_file.bcf} > sample_order.tsv

        # 2. Reported-sex mapping.
        printf "{sex_tsv_content}\\n" > sex.tsv

        # 3. Stream FORMAT/{{GT, THETA, R}} through the aggregator (bgzips inline).
        bcftools query \\
            -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID[\\t%GT,%THETA,%R]\\n' \\
            {bcf_file.bcf} \\
          | python3 {aggregator_script} \\
              --sample-order sample_order.tsv \\
              --sex-tsv sex.tsv \\
              --output {j.output.tsv_gz}

        # 4. Tabix index. -S 1 skips the (non-#-prefixed) header line.
        tabix -s 1 -b 2 -e 2 -S 1 {j.output.tsv_gz}
        """,
    )

    # Write both group members back to cloud, rooted at the prefix without suffix.
    output_root = output_stats_path.removesuffix(OUTPUT_SUFFIX)
    b.write_output(j.output, output_root)

    return j
