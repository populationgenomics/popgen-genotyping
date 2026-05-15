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

    bed_int: str = to_container(plink1_prefix_host.with_suffix('.bed'))
    king_prefix_host: Path = king_dir / prefix
    king_prefix_int: str = to_container(king_prefix_host)
    log_int: str = to_container(king_dir / f'{prefix}.log')

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
