"""
Stage: KingIbdseg reproduction.
"""

from pathlib import Path

from scripts.testing_utils import (
    DATA_DIR,
    KING_IMAGE,
    PLINK_IMAGE,
    run_docker,
    to_container,
)


def run_king_ibdseg(plink1_prefix_host: Path, prefix: str = 'king') -> dict[str, Path]:
    """
    Run KING `--ibdseg` against the merged PLINK 1.9 dataset in Docker.

    KING 2.3.2 emits autosome outputs at ``{prefix}.seg`` / ``{prefix}.segments.gz``
    and, when the input contains chrX SNPs (as it will for any reproduction whose
    --snps spans BPM index 119,078+), additional X-chromosome outputs at
    ``{prefix}X.seg`` and ``{prefix}X.segments.gz`` (no dot before the ``X``).

    Args:
        plink1_prefix_host (Path): Host path prefix for the merged PLINK 1.9
            fileset (i.e. the prefix shared by .bed/.bim/.fam).
        prefix (str): KING output prefix.

    Returns:
        dict[str, Path]: Host paths of the captured KING outputs.
    """
    print('\n>>> Stage: KingIbdseg (KING --ibdseg) <<<')
    king_dir: Path = DATA_DIR / 'king'
    king_dir.mkdir(parents=True, exist_ok=True)

    king_prefix_host: Path = king_dir / prefix
    king_prefix_int: str = to_container(king_prefix_host)
    log_int: str = to_container(king_dir / f'{prefix}.log')

    # KING only parses numeric chromosome codes (23=X, 24=Y, 26=MT); the merged
    # fileset is chr-prefixed, so recode to numeric codes with plink2 first.
    recode_prefix_host: Path = king_dir / f'{prefix}_recode'
    recode_cmd: str = (
        f"bash -c 'plink2 --bfile {to_container(plink1_prefix_host)} --allow-extra-chr "
        f"--output-chr 26 --make-bed --out {to_container(recode_prefix_host)}'"
    )
    run_docker(PLINK_IMAGE, recode_cmd)

    bed_int: str = to_container(recode_prefix_host.with_suffix('.bed'))
    cmd: str = (
        f"bash -c 'king -b {bed_int} --ibdseg --degree 3 --cpus $(nproc) "
        f"--prefix {king_prefix_int} 2>&1 | tee {log_int}'"
    )
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
