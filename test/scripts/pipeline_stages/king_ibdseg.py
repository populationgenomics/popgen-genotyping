"""
Stage: KingIbdseg reproduction.

Mirrors the production two-job chain:
1. ``run_plink_filter_for_king`` — PLINK image; build a BAFRegress remove list,
   then ``plink1.9 --remove --make-bed`` to emit a contamination-filtered .bed/
   .bim/.fam fileset.
2. ``run_king_ibdseg`` — KING image; plain ``king --ibdseg --degree 3`` over the
   filtered fileset.
"""

from pathlib import Path

from scripts.testing_utils import (
    DATA_DIR,
    KING_IMAGE,
    PLINK_IMAGE,
    run_docker,
    to_container,
)

# Keep in sync with popgen_genotyping.jobs.plink_filter_for_king_job.
BAFREGRESS_THRESHOLD: float = 0.03
BAFREGRESS_ESTIMATE_COL: str = 'baf_regress'

# Header constants kept in sync with popgen_genotyping.jobs.king_ibdseg_job.
# KING omits the .seg / .segments.gz files entirely when fewer than 2 samples
# survive the filter (it falls back to --kinship), and omits the X-chr variants
# entirely when the input has no chrX SNPs. The test harness backfills the same
# header-only placeholders the production job does so the reproduction reports
# the full output set as present.
_AUTOSOME_SEG_HEADER: str = 'FID1\tID1\tFID2\tID2\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType'
_X_SEG_HEADER: str = 'FID1\tID1\tFID2\tID2\tSex1\tSex2\tMaxIBD1\tMaxIBD2\tIBD1Seg\tIBD2Seg\tPropIBD'


def _run_plink_filter_for_king(
    plink1_prefix_host: Path,
    bafregress_paths_host: list[Path],
) -> Path:
    """
    Build a BAFRegress remove list and run plink1.9 --remove --make-bed.

    The shell logic is materialised to a script file mounted into the container
    rather than passed via ``bash -c '...'``: the awk programs use single quotes
    and embedding them directly inside the outer single-quoted argument is
    brittle (collides with the surrounding shell quoting in run_docker).

    Args:
        plink1_prefix_host: Host path prefix shared by the input .bed/.bim/.fam.
        bafregress_paths_host: Per-cohort BAFRegress output files.

    Returns:
        Host path prefix for the filtered PLINK 1.9 fileset.
    """
    print('\n>>> Stage: KingIbdseg / PlinkFilterForKing (PLINK 1.9) <<<')
    filter_dir: Path = DATA_DIR / 'king_filter'
    filter_dir.mkdir(parents=True, exist_ok=True)

    input_prefix_int: str = to_container(plink1_prefix_host)
    fam_int: str = to_container(plink1_prefix_host.with_suffix('.fam'))
    filtered_prefix_host: Path = filter_dir / 'merged_filtered'
    filtered_prefix_int: str = to_container(filtered_prefix_host)
    baf_inputs: str = ' '.join(to_container(p) for p in bafregress_paths_host)

    script_host: Path = filter_dir / 'plink_filter_run.sh'
    script_int: str = to_container(script_host)
    script_host.write_text(
        f"""#!/bin/bash
set -euo pipefail

awk '{{print $1"\\t"$2}}' {fam_int} | sort -k2,2 > /tmp/all_samples.tsv

: > /tmp/good_iids.txt
for baf in {baf_inputs}; do
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
    ' "$baf" >> /tmp/good_iids.txt
done
sort -u /tmp/good_iids.txt -o /tmp/good_iids.txt

awk 'NR==FNR{{good[$1]=1; next}} !($2 in good) {{print $1"\\t"$2}}' \\
    /tmp/good_iids.txt /tmp/all_samples.tsv > /tmp/remove_samples.tsv

n_total=$(wc -l < /tmp/all_samples.tsv)
n_remove=$(wc -l < /tmp/remove_samples.tsv)
echo "BAFRegress contamination filter: excluding $n_remove / $n_total samples" \\
     "(estimate > {BAFREGRESS_THRESHOLD}, non-numeric, or absent)"

plink --bfile {input_prefix_int} --allow-extra-chr \\
    --remove /tmp/remove_samples.tsv \\
    --keep-allele-order \\
    --make-bed --out {filtered_prefix_int}
""",
    )
    run_docker(PLINK_IMAGE, f'bash {script_int}')

    return filtered_prefix_host


def run_king_ibdseg(
    plink1_prefix_host: Path,
    bafregress_paths_host: list[Path],
    prefix: str = 'king',
) -> dict[str, Path]:
    """
    Run the KingIbdseg two-job chain (PLINK filter → KING --ibdseg) in Docker.

    KING 2.3.2 emits autosome outputs at ``{prefix}.seg`` / ``{prefix}.segments.gz``
    and, when the input contains chrX SNPs (as it will for any reproduction whose
    --snps spans BPM index 119,078+), additional X-chromosome outputs at
    ``{prefix}X.seg`` and ``{prefix}X.segments.gz`` (no dot before the ``X``).

    Args:
        plink1_prefix_host: Host path prefix for the merged PLINK 1.9 fileset
            (i.e. the prefix shared by .bed/.bim/.fam).
        bafregress_paths_host: Host paths of the BAFRegress ``.txt`` outputs,
            one per cohort.
        prefix: KING output prefix.

    Returns:
        Host paths of the captured KING outputs.
    """
    filtered_prefix_host: Path = _run_plink_filter_for_king(plink1_prefix_host, bafregress_paths_host)

    print('\n>>> Stage: KingIbdseg / KING --ibdseg (KING image) <<<')
    king_dir: Path = DATA_DIR / 'king'
    king_dir.mkdir(parents=True, exist_ok=True)

    filtered_bed_int: str = to_container(filtered_prefix_host.with_suffix('.bed'))
    king_prefix_host: Path = king_dir / prefix
    king_prefix_int: str = to_container(king_prefix_host)
    log_host: Path = king_dir / f'{prefix}.log'
    log_int: str = to_container(log_host)

    seg_int: str = to_container(king_prefix_host.with_suffix('.seg'))
    segments_gz_int: str = to_container(king_prefix_host.with_suffix('.segments.gz'))
    seg_x_int: str = to_container(king_prefix_host.parent / f'{prefix}X.seg')
    segments_gz_x_int: str = to_container(king_prefix_host.parent / f'{prefix}X.segments.gz')

    script_host: Path = king_dir / 'king_run.sh'
    script_int: str = to_container(script_host)
    script_host.write_text(
        f"""#!/bin/bash
set -euo pipefail

king -b {filtered_bed_int} --ibdseg --degree 3 --cpus $(nproc) \\
    --prefix {king_prefix_int} 2>&1 | tee {log_int}

if [ ! -s {seg_int} ]; then
    printf '%s\\n' '{_AUTOSOME_SEG_HEADER}' > {seg_int}
fi
if [ ! -s {segments_gz_int} ]; then
    : | gzip -c > {segments_gz_int}
fi
if [ ! -s {seg_x_int} ]; then
    printf '%s\\n' '{_X_SEG_HEADER}' > {seg_x_int}
fi
if [ ! -s {segments_gz_x_int} ]; then
    : | gzip -c > {segments_gz_x_int}
fi
""",
    )
    run_docker(KING_IMAGE, f'bash {script_int}')

    outputs: dict[str, Path] = {
        'seg': king_prefix_host.with_suffix('.seg'),
        'segments': king_prefix_host.with_suffix('.segments.gz'),
        'seg_x': king_prefix_host.parent / f'{prefix}X.seg',
        'segments_x': king_prefix_host.parent / f'{prefix}X.segments.gz',
        'allsegs': king_prefix_host.parent / f'{prefix}allsegs.txt',
        'log': log_host,
    }
    return outputs
