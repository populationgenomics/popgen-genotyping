"""
Job logic for filtering a merged PLINK 1.9 fileset to exclude samples whose
BAFRegress contamination estimate fails the threshold, ahead of the KING
``--ibdseg`` relatedness calculation.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob
    from hailtop.batch.resource import ResourceGroup


# Column name carrying the BAFRegress contamination estimate in the
# ``bcftools +BAFregress`` plugin output (alongside ``sample_id`` and ``Nhom``).
BAFREGRESS_ESTIMATE_COL: str = 'baf_regress'

# Maximum BAFRegress contamination estimate accepted before a sample is excluded
# from the IBD calculation. The 3% cutoff follows the BAFRegress convention
# (Jun et al. 2012, doi:10.1016/j.ajhg.2012.09.004); samples whose estimate is
# above this, non-numeric (e.g. ``-nan`` from a failed fit), or absent from any
# BAFRegress output are dropped.
BAFREGRESS_THRESHOLD: float = 0.03


def run_plink_filter_for_king(
    bed_path: str,
    bim_path: str,
    fam_path: str,
    bafregress_paths: list[str],
    job_name: str = 'plink_filter_for_king',
) -> tuple[BashJob, ResourceGroup]:
    """
    Filter a merged PLINK 1.9 fileset to drop BAFRegress contamination failures.

    The job builds a remove-list of samples whose BAFRegress ``baf_regress``
    column is above ``BAFREGRESS_THRESHOLD``, non-numeric (e.g. ``-nan``), or
    absent from every per-cohort BAFRegress output. It then runs
    ``plink1.9 --remove ... --make-bed`` to materialise a contamination-filtered
    ``.bed``/``.bim``/``.fam`` fileset that is consumed by the downstream KING
    job through Hail Batch in-memory wiring (no GCS round-trip).

    The ``--make-bed`` re-pass has a useful side-effect on small synthetic
    cohorts: it normalises the major-allele assignment, which suppresses KING's
    "Too many first alleles as major" warning seen in the local Docker harness.

    Args:
        bed_path: Cloud path to the input PLINK 1.9 .bed file.
        bim_path: Cloud path to the input PLINK 1.9 .bim file.
        fam_path: Cloud path to the input PLINK 1.9 .fam file.
        bafregress_paths: Per-cohort BAFRegress output files. Samples whose
            ``baf_regress`` value exceeds ``BAFREGRESS_THRESHOLD``, is
            non-numeric, or is absent from every file are dropped.
        job_name: Name for the Hail Batch job.

    Returns:
        Tuple of the queued Hail Batch BashJob and the filtered PLINK 1.9
        ResourceGroup (with ``bed``/``bim``/``fam`` keys).
    """
    b = get_batch()
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'king_ibdseg'],
        image=config_retrieve(['workflow', 'plink_image']),
        default_cpu=2,
        default_storage='50G',
    )

    plink_input = b.read_input_group(bed=bed_path, bim=bim_path, fam=fam_path)
    bafregress_files = [b.read_input(p) for p in bafregress_paths]

    j.declare_resource_group(
        filtered_plink={
            'bed': '{root}.bed',
            'bim': '{root}.bim',
            'fam': '{root}.fam',
        },
    )

    # Build the remove list: (all .fam IIDs) \ (IIDs with a numeric BAFRegress
    # estimate <= threshold). The set-diff handles the absent-from-BAFRegress
    # case in the same operation.
    bafregress_inputs_block = ' '.join(bafregress_files) if bafregress_files else ''
    j.command(
        f"""
        set -euo pipefail

        # All FID/IID pairs in the merged PLINK fileset, sorted by IID.
        awk '{{print $1"\\t"$2}}' {plink_input.fam} | sort -k2,2 > all_samples.tsv

        # IIDs with a numeric BAFRegress estimate <= {BAFREGRESS_THRESHOLD}.
        : > good_iids.txt
        for baf in {bafregress_inputs_block}; do
            awk -v thresh={BAFREGRESS_THRESHOLD} -v est_name={BAFREGRESS_ESTIMATE_COL} '
                NR==1 {{
                    sid_col = 0; est_col = 0
                    for (i=1; i<=NF; i++) {{
                        if ($i == "sample_id") sid_col = i
                        else if ($i == est_name) est_col = i
                    }}
                    if (sid_col == 0 || est_col == 0) {{
                        print "ERROR: missing sample_id or " est_name " column in " FILENAME > "/dev/stderr"
                        exit 1
                    }}
                    next
                }}
                $est_col ~ /^-?[0-9]+\\.?[0-9]*([eE][-+]?[0-9]+)?$/ && ($est_col + 0) <= thresh {{
                    print $sid_col
                }}
            ' "$baf" >> good_iids.txt
        done
        sort -u good_iids.txt -o good_iids.txt

        # remove = samples in .fam whose IID is not in good_iids.
        awk 'NR==FNR{{good[$1]=1; next}} !($2 in good) {{print $1"\\t"$2}}' \\
            good_iids.txt all_samples.tsv > remove_samples.tsv

        n_total=$(wc -l < all_samples.tsv)
        n_remove=$(wc -l < remove_samples.tsv)
        echo "BAFRegress contamination filter: excluding $n_remove / $n_total samples" \\
             "(estimate > {BAFREGRESS_THRESHOLD}, non-numeric, or absent)"

        plink --bfile {plink_input} --allow-extra-chr \\
            --remove remove_samples.tsv \\
            --keep-allele-order \\
            --make-bed --out {j.filtered_plink}
        """,
    )

    return j, j.filtered_plink
