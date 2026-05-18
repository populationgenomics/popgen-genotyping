"""
Stage: KingIbdseg reproduction.
"""

from pathlib import Path

from scripts.testing_utils import (
    DATA_DIR,
    KING_IMAGE,
    run_docker,
    to_container,
)

# Keep in sync with popgen_genotyping.jobs.king_ibdseg_job.BAFREGRESS_THRESHOLD.
BAFREGRESS_THRESHOLD: float = 0.03


def run_king_ibdseg(
    plink1_prefix_host: Path,
    bafregress_paths_host: list[Path],
    prefix: str = 'king',
) -> dict[str, Path]:
    """
    Run KING `--ibdseg` against the merged PLINK 1.9 dataset in Docker.

    Mirrors the production stage: builds a ``--remove`` list of samples whose
    BAFRegress ``estimate`` is non-numeric, missing, or above
    ``BAFREGRESS_THRESHOLD``, then runs ``king --ibdseg --degree 3``.

    KING 2.3.2 emits autosome outputs at ``{prefix}.seg`` / ``{prefix}.segments.gz``
    and, when the input contains chrX SNPs (as it will for any reproduction whose
    --snps spans BPM index 119,078+), additional X-chromosome outputs at
    ``{prefix}X.seg`` and ``{prefix}X.segments.gz`` (no dot before the ``X``).

    Args:
        plink1_prefix_host (Path): Host path prefix for the merged PLINK 1.9
            fileset (i.e. the prefix shared by .bed/.bim/.fam).
        bafregress_paths_host (list[Path]): Host paths of the BAFRegress
            ``.txt`` outputs, one per cohort.
        prefix (str): KING output prefix.

    Returns:
        dict[str, Path]: Host paths of the captured KING outputs.
    """
    print('\n>>> Stage: KingIbdseg (KING --ibdseg) <<<')
    king_dir: Path = DATA_DIR / 'king'
    king_dir.mkdir(parents=True, exist_ok=True)

    bed_int: str = to_container(plink1_prefix_host.with_suffix('.bed'))
    fam_int: str = to_container(plink1_prefix_host.with_suffix('.fam'))
    king_prefix_host: Path = king_dir / prefix
    king_prefix_int: str = to_container(king_prefix_host)
    log_int: str = to_container(king_dir / f'{prefix}.log')

    baf_inputs: str = ' '.join(to_container(p) for p in bafregress_paths_host)
    awk_program: str = (
        'NR==1 {'
        ' sid_col = 0; est_col = 0;'
        ' for (i=1; i<=NF; i++) {'
        '   if ($i == "sample_id") sid_col = i;'
        '   else if ($i == "estimate") est_col = i;'
        ' }'
        ' if (sid_col == 0 || est_col == 0) {'
        '   print "ERROR: missing sample_id or estimate column in " FILENAME > "/dev/stderr";'
        '   exit 1;'
        ' }'
        ' next;'
        '} '
        '$est_col ~ /^-?[0-9]+\\.?[0-9]*([eE][-+]?[0-9]+)?$/ && ($est_col + 0) <= thresh {'
        ' print $sid_col;'
        '}'
    )
    filter_block: str = (
        f'awk \'{{print $1"\\t"$2}}\' {fam_int} | sort -k2,2 > /tmp/all_samples.tsv && '
        f': > /tmp/good_iids.txt && '
        f'for baf in {baf_inputs}; do '
        f'  awk -v thresh={BAFREGRESS_THRESHOLD} \'{awk_program}\' "$baf" >> /tmp/good_iids.txt; '
        f'done && '
        f'sort -u /tmp/good_iids.txt -o /tmp/good_iids.txt && '
        f'awk \'NR==FNR{{good[$1]=1; next}} !($2 in good) {{print $1"\\t"$2}}\' '
        f'  /tmp/good_iids.txt /tmp/all_samples.tsv > /tmp/remove_samples.tsv && '
        f'n_total=$(wc -l < /tmp/all_samples.tsv) && '
        f'n_remove=$(wc -l < /tmp/remove_samples.tsv) && '
        f'echo "BAFRegress contamination filter: excluding $n_remove / $n_total samples '
        f'(estimate > {BAFREGRESS_THRESHOLD}, non-numeric, or absent)"'
    )
    king_block: str = (
        f'if [ "$n_remove" -gt 0 ]; then '
        f'  king -b {bed_int} --ibdseg --degree 3 --cpus $(nproc) '
        f'    --remove /tmp/remove_samples.tsv --prefix {king_prefix_int} 2>&1 | tee {log_int}; '
        f'else '
        f'  king -b {bed_int} --ibdseg --degree 3 --cpus $(nproc) '
        f'    --prefix {king_prefix_int} 2>&1 | tee {log_int}; '
        f'fi'
    )
    cmd: str = f"bash -c '{filter_block} && {king_block}'"
    run_docker(KING_IMAGE, cmd)

    outputs: dict[str, Path] = {
        'seg': king_prefix_host.with_suffix('.seg'),
        'segments': king_prefix_host.with_suffix('.segments.gz'),
        'seg_x': king_prefix_host.parent / f'{prefix}X.seg',
        'segments_x': king_prefix_host.parent / f'{prefix}X.segments.gz',
        'allsegs': king_prefix_host.parent / f'{prefix}allsegs.txt',
        'log': king_prefix_host.with_suffix('.log'),
    }
    return outputs
