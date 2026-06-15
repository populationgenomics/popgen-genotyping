#!/usr/bin/env python3

"""
Regenerate the BPM locus-index to chromosome cache from a real Illumina BPM.

The synthetic GTC generator needs to know which BPM indices fall on chrX (and
chrY/MT) so it can emit the correct ploidy per sample sex. The cache it reads
is `test/data/bpm_chrom.txt.gz`: one chromosome string per line, where line N
(1-indexed) corresponds to BPM locus index N.

This script regenerates the cache by invoking `bcftools +gtc2vcf -b <bpm>
--verbose` inside the production bcftools image, parsing the [Assay] CSV
section, and writing the gzipped one-column file. Run this only when the BPM
itself changes (i.e., a new array version is adopted).
"""

import argparse
import csv
import gzip
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

BCFTOOLS_IMAGE: str = 'australia-southeast1-docker.pkg.dev/cpg-common/images/bcftools:1.23-1'


def _dump_bpm(bpm_path: Path, work_dir: Path) -> Path:
    """
    Dump the BPM manifest to a text [Assay] table using bcftools +gtc2vcf.

    Args:
        bpm_path: Path to the Illumina BPM file on disk.
        work_dir: Working directory to mount into the bcftools container.

    Returns:
        Path to the verbose text dump (containing the [Assay] CSV table).
    """
    staged_bpm: Path = work_dir / bpm_path.name
    if not staged_bpm.exists():
        shutil.copy2(bpm_path, staged_bpm)

    dump_path: Path = work_dir / 'bpm_verbose.txt'
    cmd: list[str] = [
        'docker',
        'run',
        '--rm',
        '-v',
        f'{work_dir}:/data',
        '-w',
        '/data',
        BCFTOOLS_IMAGE,
        'bash',
        '-c',
        f'bcftools +gtc2vcf -b {staged_bpm.name} --verbose > {dump_path.name}',
    ]
    subprocess.run(cmd, check=True)  # noqa: S603
    return dump_path


def _parse_assay_table(dump_path: Path) -> list[str]:
    """
    Parse the [Assay] CSV table from the bcftools +gtc2vcf verbose dump.

    Args:
        dump_path: Path to the verbose text dump.

    Returns:
        Ordered list of chromosome strings, indexed by (BPM_index - 1).
    """
    chroms: list[str] = []
    header_seen: bool = False
    index_col: int = -1
    chr_col: int = -1
    expected_index: int = 1
    with open(dump_path, encoding='utf-8') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            if not header_seen:
                if row[0] == 'Index' and 'Chr' in row:
                    index_col = row.index('Index')
                    chr_col = row.index('Chr')
                    header_seen = True
                continue
            # After the [Assay] header, each row whose first field is a number
            # is a locus. Stop when we hit [Controls] or similar trailing sections.
            if not row[0].isdigit():
                break
            idx: int = int(row[index_col])
            if idx != expected_index:
                raise ValueError(
                    f'BPM Index column is not monotonic 1..N at row {expected_index}: '
                    f'expected {expected_index}, got {idx}',
                )
            chroms.append(row[chr_col])
            expected_index += 1

    if not chroms:
        raise ValueError(f'No assay rows parsed from {dump_path}; check bcftools output')
    return chroms


def _write_cache(chroms: list[str], cache_path: Path) -> None:
    """
    Write the chromosome list to a gzipped one-chromosome-per-line file.

    Args:
        chroms: Ordered list of chromosome strings (BPM_index 1..N).
        cache_path: Destination path (`.txt.gz`).
    """
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(cache_path, 'wt', encoding='utf-8') as f:
        for chrom in chroms:
            f.write(f'{chrom}\n')


def main() -> None:
    """
    CLI entry point: regenerate test/data/bpm_chrom.txt.gz from a BPM file.
    """
    repo_root: Path = Path(__file__).resolve().parents[2]
    default_cache: Path = repo_root / 'test' / 'data' / 'bpm_chrom.txt.gz'

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bpm', type=Path, help='Path to the Illumina BPM file')
    parser.add_argument(
        '--output',
        type=Path,
        default=default_cache,
        help=f'Output cache path (default: {default_cache})',
    )
    args: argparse.Namespace = parser.parse_args()

    if not args.bpm.exists():
        sys.exit(f'BPM not found: {args.bpm}')

    with tempfile.TemporaryDirectory() as tmpdir_str:
        tmpdir: Path = Path(tmpdir_str)
        dump_path: Path = _dump_bpm(args.bpm, tmpdir)
        chroms: list[str] = _parse_assay_table(dump_path)

    _write_cache(chroms, args.output)
    print(f'Wrote {len(chroms)} loci to {args.output}')


if __name__ == '__main__':
    main()
