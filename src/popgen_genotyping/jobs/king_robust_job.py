"""
Job logic for pairwise KING-robust kinship with ``plink2 --make-king-table``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def run_king_robust(
    pgen_path: str,
    pvar_path: str,
    psam_path: str,
    remove_list_path: str,
    kinship_min: float,
    output_kin0_path: str,
    output_log_path: str,
    job_name: str = 'king_robust',
) -> BashJob:
    """
    Write the pairwise KING-robust kinship table of the merged cohort with plink2.

    ``plink2 --make-king-table`` reports the KING-robust between-family kinship estimator of
    Manichaikul et al. (2010), doi:10.1093/bioinformatics/btq559, for every pair at or above
    ``kinship_min``. The estimator only considers autosomes; ``--autosome`` is passed anyway so
    that is explicit rather than a plink2 behaviour to remember. It runs on all the merged
    cohort's variants, with no SNP-QC filtering, to match the input of KING ``--ibdseg``. The
    samples in the ``RelatednessExclusions`` list are removed first (a header-only list removes
    nothing, so ``--remove`` is always passed).

    The table columns are fixed with ``cols=`` to ``IID1``, ``IID2``, ``NSNP``, ``HETHET``,
    ``IBS0`` and ``KINSHIP``, so the layout does not depend on whether the ``.psam`` has a
    ``FID`` column. plink2 writes a header-only ``.kin0`` when no pair reaches ``kinship_min``,
    so a missing or empty ``.kin0`` or log means plink2 did not finish, and the job fails. plink2
    can only compress the table with zstd, which the driver image's Python cannot read without an
    extra package, so the job gzips it instead. Every output file is named in the command: Hail
    Batch refuses to write out a resource the command does not mention.

    Args:
        pgen_path: Cloud path to the merged cohort ``.pgen``.
        pvar_path: Cloud path to its ``.pvar``.
        psam_path: Cloud path to its ``.psam``.
        remove_list_path: Cloud path to the plink2 ``--remove`` list (``#IID`` header, one
            IID per line) written by ``RelatednessExclusions``.
        kinship_min: Minimum kinship of a pair to be written (``--king-table-filter``).
        output_kin0_path: Cloud path for the gzipped ``.kin0`` pairwise table.
        output_log_path: Cloud path for the plink2 log.
        job_name: Name for the Hail Batch job.

    Returns:
        The queued Hail Batch job.
    """
    b = get_batch()
    j = register_job(
        batch=b,
        job_name=job_name,
        config_path=['popgen_genotyping', 'king_robust'],
        image=config_retrieve(['workflow', 'plink_image']),
        default_cpu=8,
        default_storage='20G',
    )

    pgen_files = b.read_input_group(pgen=pgen_path, pvar=pvar_path, psam=psam_path)
    remove_list = b.read_input(remove_list_path)

    j.declare_resource_group(
        king_robust={
            'kin0': '{root}.kin0',
            'kin0_gz': '{root}.kin0.gz',
            'log': '{root}.log',
        },
    )

    j.command(
        f"""
        set -ex
        plink2 --pfile {pgen_files} \\
            --remove {remove_list} \\
            --autosome \\
            --make-king-table cols=id,nsnp,hethet,ibs0,kinship \\
            --king-table-filter {kinship_min} \\
            --out {j.king_robust}
        for output in {j.king_robust.kin0} {j.king_robust.log}; do
            if [ ! -s "$output" ]; then
                echo "plink2 wrote no $output" >&2
                exit 1
            fi
        done
        gzip -c {j.king_robust.kin0} > {j.king_robust.kin0_gz}
        """,
    )

    b.write_output(j.king_robust.kin0_gz, output_kin0_path)
    b.write_output(j.king_robust.log, output_log_path)

    return j
