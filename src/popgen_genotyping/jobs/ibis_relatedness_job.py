"""
Job logic for the per-pair relatedness table from IBIS cM segments.

KING --ibdseg can only filter segments by length in Mb. A population with background
relatedness (a bottleneck, say) gets "3rd-degree" KING pairs made of many short segments, which
look different from recent relatives, whose IBD is a few long segments. IBIS calls IBD segments
in genetic-map cM, and the table re-estimates each pair's kinship from IBD runs of at least
several cM lengths, so analysts can tell recent relatedness from background IBD.

Five jobs of four kinds run in the batch:

1. plink2: from the aggregate PLINK2 fileset, drop the samples ``RelatednessExclusions`` left out
   (contamination or missingness), then write two PLINK 1 filesets with deCODE 2019 centimorgan
   positions: the autosomes, and chromosome X (unknown-sex samples removed, because IBIS needs
   sex to read X). PAR1 and PAR2 variants are in neither. Chromosomes are written as 1-22 and X
   (``--output-chr MT``): IBIS recognises X only under the name "X" (it passes "X" to its reader),
   and a chr-prefixed or numeric (23) name would make it analyse X as an autosome.
2. prepare-maps (Python): plink2 extrapolates variants before a chromosome's first map row to a
   slightly negative cM, so shift each such chromosome up until its first marker is at 0 cM. The
   shift keeps every distance between markers. Clamping those cM to 0 instead breaks IBIS: with
   the first two markers of chromosome 1 both at 0 cM it takes the map to be missing and then
   exits with "some markers have non-zero genetic position". Then check cM is non-decreasing,
   every chromosome has a map, and measure the autosomal map length.
3. IBIS, autosomes and X as two separate jobs. The X run is kept apart so X never enters the
   autosomal kinship. IBIS reads X from the .fam sex, sets male heterozygous X calls to missing
   and otherwise treats X like an autosome; its X handling is undocumented, so treat the X
   columns as informative (for example to order a pedigree), not as validated kinship.
   IBIS runs with ``-min_l2 5`` rather than its 2 cM default: parent-offspring pairs, who cannot be
   IBD2, otherwise get a few false IBD2 calls of 2 to 3 cM each.
4. summarise (Python): the table (see ``ibis_relatedness_summary.build_relatedness_table``).

IBIS reports one continuous stretch of IBD as several abutting segments: IBD2 stretches are
separate segments, so every IBD1/IBD2 boundary starts a new one, and now and then two IBD1
segments abut with no gap. The summarise job joins the segments of a pair on a chromosome into
runs: a segment joins the previous one's run when the gap between the previous end and its start
is under ``join_gap_cm``, whatever their IBD types. A run is as long as its genetic span and its
kinship contribution is the sum over its segments of length x weight (IBD1 1/4, IBD2 1/2).
``min_seg_cm`` and ``cutoffs_cm`` apply to run lengths, and the ``IBIS_*`` and ``X_*`` columns
describe runs.

Outputs:

    pairs table (.tsv.gz)   one row per pair related by KING or by IBIS autosomal segments
    pairs sidecar (.json)   run parameters, input paths, the autosomal map length, counts
    autosomal IBIS run      .seg.gz, .coef.gz, .log
    X IBIS run              .seg.gz, .coef.gz, .log

IBIS .coef files list every pair (not only related ones), and their fourth and fifth columns are
swapped relative to the header: the data are kinship, segment count, IBD2 fraction, degree.

The genetic map is one deCODE 2019 file per chromosome (1-22 and X), resolved from the references
config. The plink2 image is the 2.0 beta 1 build (``workflow.plink2_cm_map_image``), because the
plink2 in the plink image has no ``--cm-map``.
"""

from __future__ import annotations

import json
from importlib.resources import files
from shlex import quote
from typing import TYPE_CHECKING, Any, Final

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, get_driver_image, reference_path
from cpg_utils.hail_batch import get_batch

from popgen_genotyping.utils import register_job

if TYPE_CHECKING:
    from hailtop.batch import Batch
    from hailtop.batch.job import BashJob

MAP_CHROMOSOMES: Final[tuple[str, ...]] = (*(str(chromosome) for chromosome in range(1, 23)), 'X')

# Name each genetic map is linked under in the plink2 job, so that plink2's '@' pattern (which
# stands for the chromosome code) finds one file per chromosome. The reference keys, which name the
# map in use, come from config.
GENETIC_MAP_LINK: Final = 'genetic_map_chr{chromosome}.txt.gz'

CONFIG_PATH: Final[list[str]] = ['popgen_genotyping', 'ibis_relatedness']
RESOURCE_KEYS: Final[tuple[str, ...]] = ('cpu', 'memory', 'storage')


def resolve_genetic_map_paths(reference_key_template: str) -> dict[str, str]:
    """
    The per-chromosome genetic map path for each of 1-22 and X, from the references config.

    Args:
        reference_key_template: References-config key of one chromosome's map, with a
            ``{chromosome}`` placeholder that is filled with 1-22 and X, for example
            ``genetic_maps_decode2019/sexavg_grch38_chr{chromosome}``.

    Returns:
        Chromosome name (``'1'``-``'22'``, ``'X'``) to the cloud path of its map.

    Raises:
        ValueError: If the template has no ``{chromosome}`` placeholder.
    """
    if '{chromosome}' not in reference_key_template:
        raise ValueError(f'genetic map reference key {reference_key_template!r} has no {{chromosome}} placeholder')
    return {
        chromosome: reference_path(reference_key_template.format(chromosome=chromosome))
        for chromosome in MAP_CHROMOSOMES
    }


def check_inputs_exist(paths: list[str]) -> None:
    """
    Raise if any input is missing, listing all of them, so a missing reference fails before any job is queued.
    """
    missing = [path for path in paths if not to_path(path).exists()]
    if missing:
        raise FileNotFoundError('missing input(s):\n' + '\n'.join(missing))


def check_parameters(min_seg_cm: float, join_gap_cm: float, cutoffs_cm: list[float]) -> None:
    """
    Raise if the segment length floor, join gap or cutoffs cannot be used.

    Args:
        min_seg_cm: IBIS ``-min_l``: shortest run kept (cM).
        join_gap_cm: Segments closer than this (cM) join the same run.
        cutoffs_cm: Run lengths (cM) the kinship is re-estimated at.

    Raises:
        ValueError: If there are no cutoffs, a cutoff is below ``min_seg_cm`` (IBIS never reports
            it) or repeated, or ``join_gap_cm`` is negative.
    """
    if not cutoffs_cm:
        raise ValueError('at least one cutoff is required')
    below_floor = [cutoff for cutoff in cutoffs_cm if cutoff < min_seg_cm]
    if below_floor:
        raise ValueError(f'cutoffs_cm {below_floor} are below min_seg_cm {min_seg_cm:g}, which IBIS never reports')
    if len(set(cutoffs_cm)) != len(cutoffs_cm):
        raise ValueError(f'cutoffs_cm {cutoffs_cm} repeats a value')
    if join_gap_cm < 0:
        raise ValueError(f'join_gap_cm must not be negative, got {join_gap_cm:g}')


def plink_command(pfile: Any, remove_list: Any, maps: dict[str, Any], autosomes: Any, x: Any) -> str:
    """
    Shell script for the plink2 job: drop the excluded samples and write the cM-mapped filesets.
    """
    link_maps = '\n'.join(
        f'ln -s {maps[chromosome]} maps/{GENETIC_MAP_LINK.format(chromosome=chromosome)}'
        for chromosome in MAP_CHROMOSOMES
    )
    # The '@' in --cm-map stands for the chromosome code. A genome-wide Eagle-style map would
    # be one file, but plink2 refuses it ("Chromosome 'X' is split") for a chromosome subset, and
    # the autosomes and X are written as separate filesets.
    return f"""
set -euo pipefail

# The remove list is a header, then the IID of every sample left out of relatedness estimation. It
# is header-only when nothing was excluded, which removes nothing. The sample counts below rely
# on the one header line.
if [ "$(head -n 1 {remove_list})" != "#IID" ]; then
    echo "unexpected remove-list header: $(head -n 1 {remove_list})" >&2
    exit 1
fi

mkdir maps
{link_maps}
cm_map='maps/{GENETIC_MAP_LINK.format(chromosome='@')}'

plink2 --pfile {pfile} --remove {remove_list} --chr 1-22 \\
    --cm-map "$cm_map" --output-chr MT --make-bed --out {autosomes}
plink2 --pfile {pfile} --remove {remove_list} --remove-nosex --chr X \\
    --cm-map "$cm_map" --output-chr MT --make-bed --out {x}

n_input=$(($(wc -l < {pfile.psam}) - 1))
n_excluded=$(($(wc -l < {remove_list}) - 1))
n_autosomal=$(wc -l < {autosomes.fam})
n_x=$(wc -l < {x.fam})
echo "samples: $n_input in the aggregate, $n_excluded excluded, $n_autosomal autosomal, $n_x on X"
if [ "$n_autosomal" -ne "$((n_input - n_excluded))" ]; then
    echo "expected $((n_input - n_excluded)) samples after --remove, got $n_autosomal" >&2
    exit 1
fi
if [ "$n_x" -eq 0 ] || [ "$n_x" -gt "$n_autosomal" ]; then
    echo "implausible X sample count $n_x (autosomal $n_autosomal)" >&2
    exit 1
fi
echo "variants: $(wc -l < {autosomes.bim}) autosomal, $(wc -l < {x.bim}) on X"
"""


def ibis_command(
    bed: Any,
    bim: Any,
    fam: Any,
    out: Any,
    min_seg_cm: float,
    min_ibd2_cm: float,
    threads: int,
    chromosome: str | None = None,
) -> str:
    """
    Shell script for an IBIS job; the -chr option restricts the run to one chromosome (X).
    """
    chromosome_option = f'-chr {chromosome}' if chromosome else ''
    # bed, bim and fam are IBIS's first three arguments so the corrected .bim can sit beside the
    # plink2 .bed and .fam under different names. -2 adds IBD2 segments, kept down to -min_l2;
    # -gzip compresses the .seg and .coef; the log, which holds "Total Genetic Length in use", is
    # kept for the summary step to cross-check against the .bim map length. pipefail: IBIS exits
    # non-zero on error, and tee would otherwise hide that. Each output file is named in the command
    # so Batch can write it out; an empty .gz still has a header, so a missing or zero-byte file
    # means IBIS did not finish.
    # -chr X is redundant on the X-only fileset, but IBIS then exits with "no markers to process"
    # if the .bim does not name the chromosome exactly X.
    return f"""
set -euo pipefail
ibis {bed} {bim} {fam} -min_l {min_seg_cm:g} -min_l2 {min_ibd2_cm:g} -2 -printCoef -gzip -t {threads} \\
    {chromosome_option} -f {out} 2>&1 \\
    | tee {out.log}
for output in {out.seg} {out.coef}; do
    if [ ! -s "$output" ]; then
        echo "IBIS wrote no $output" >&2
        exit 1
    fi
done
"""


def _register_job(batch: Batch, job_name: str, section: str, image: str) -> BashJob:
    """
    A job with the resources of ``[popgen_genotyping.ibis_relatedness.<section>]``, all required.

    ``register_job`` falls back to code defaults for a missing resource key; here a missing key
    raises, so every job's size is a deliberate config choice.

    Raises:
        ConfigError: If the section lacks ``cpu``, ``memory`` or ``storage``.
    """
    config_path = [*CONFIG_PATH, section]
    for key in RESOURCE_KEYS:
        config_retrieve([*config_path, key])
    return register_job(batch=batch, job_name=job_name, config_path=config_path, image=image)


def _register_ibis_job(batch: Batch, job_name: str, section: str, image: str) -> BashJob:
    """
    An IBIS job (see ``_register_job``) with its output group (seg, coef, log) declared; the caller adds the command.
    """
    job = _register_job(batch, job_name, section, image)
    job.declare_resource_group(ibis={'seg': '{root}.seg.gz', 'coef': '{root}.coef.gz', 'log': '{root}.log'})
    return job


def run_ibis_relatedness(
    pgen_path: str,
    pvar_path: str,
    psam_path: str,
    king_seg_path: str,
    remove_list_path: str,
    genetic_map_paths: dict[str, str],
    min_seg_cm: float,
    min_ibd2_cm: float,
    join_gap_cm: float,
    cutoffs_cm: list[float],
    output_pairs_path: str,
    output_sidecar_path: str,
    output_autosomes_seg_path: str,
    output_autosomes_coef_path: str,
    output_autosomes_log_path: str,
    output_x_seg_path: str,
    output_x_coef_path: str,
    output_x_log_path: str,
    job_name: str = 'ibis_relatedness',
) -> list[BashJob]:
    """
    Queue the plink2, prepare-maps, IBIS (autosomes and X) and summarise jobs.

    See the module docstring for what each job does and why. Each job's resources are read from
    ``[popgen_genotyping.ibis_relatedness.<job>]`` (``plink2``, ``prepare_maps``,
    ``ibis_autosomes``, ``ibis_x``, ``summarise``); the images are ``workflow.plink2_cm_map_image``,
    ``workflow.ibis_image`` and the driver image (for the two Python steps, which need only
    pandas).

    Args:
        pgen_path: Cloud path to the aggregate ``.pgen``.
        pvar_path: Cloud path to its ``.pvar``.
        psam_path: Cloud path to its ``.psam``.
        king_seg_path: Cloud path to the KING ``--ibdseg`` autosomal ``.seg`` for the same aggregate.
        remove_list_path: Cloud path to the plink2 ``--remove`` list (``#IID`` header, one IID per
            line) written by ``RelatednessExclusions``.
        genetic_map_paths: Cloud path of the genetic map of each chromosome 1-22 and X.
        min_seg_cm: IBIS ``-min_l``: shortest segment, and the shortest run kept (cM).
        min_ibd2_cm: IBIS ``-min_l2``: shortest IBD2 segment (cM).
        join_gap_cm: Segments of a pair on a chromosome with a gap below this (cM) are joined into one run.
        cutoffs_cm: Run lengths (cM) the kinship is re-estimated at; each at least ``min_seg_cm``.
        output_pairs_path: Cloud path for the gzipped pairs table.
        output_sidecar_path: Cloud path for the pairs table's sidecar JSON.
        output_autosomes_seg_path: Cloud path for the autosomal IBIS ``.seg.gz``.
        output_autosomes_coef_path: Cloud path for the autosomal IBIS ``.coef.gz``.
        output_autosomes_log_path: Cloud path for the autosomal IBIS log.
        output_x_seg_path: Cloud path for the X IBIS ``.seg.gz``.
        output_x_coef_path: Cloud path for the X IBIS ``.coef.gz``.
        output_x_log_path: Cloud path for the X IBIS log.
        job_name: Name prefix for the Hail Batch jobs.

    Returns:
        The queued jobs: ``[plink2, prepare_maps, ibis_autosomes, ibis_x, summarise]``.

    Raises:
        ValueError: If ``genetic_map_paths`` does not hold exactly chromosomes 1-22 and X, or the
            segment length floor, join gap or cutoffs cannot be used.
        FileNotFoundError: If a genetic map does not exist.
    """
    if set(genetic_map_paths) != set(MAP_CHROMOSOMES):
        raise ValueError(
            f'genetic_map_paths holds chromosomes {sorted(genetic_map_paths)}, expected {list(MAP_CHROMOSOMES)}',
        )
    check_parameters(min_seg_cm, join_gap_cm, cutoffs_cm)
    check_inputs_exist(list(genetic_map_paths.values()))

    b = get_batch()
    images: dict[str, str] = {
        'ibis': config_retrieve(['workflow', 'ibis_image']),
        'plink2': config_retrieve(['workflow', 'plink2_cm_map_image']),
        'python': get_driver_image(),
    }
    script_path = str(files('popgen_genotyping.scripts').joinpath('ibis_relatedness_summary.py'))
    summary_script = b.read_input(script_path)

    pfile = b.read_input_group(pgen=pgen_path, pvar=pvar_path, psam=psam_path)
    remove_list = b.read_input(remove_list_path)
    map_inputs = {chromosome: b.read_input(path) for chromosome, path in genetic_map_paths.items()}

    plink = _register_job(b, f'{job_name}_plink2', 'plink2', images['plink2'])
    plink.declare_resource_group(
        autosomes={'bed': '{root}.bed', 'bim': '{root}.bim', 'fam': '{root}.fam'},
        x={'bed': '{root}.bed', 'bim': '{root}.bim', 'fam': '{root}.fam'},
    )
    plink.command(plink_command(pfile, remove_list, map_inputs, plink.autosomes, plink.x))

    prepare = _register_job(b, f'{job_name}_prepare_maps', 'prepare_maps', images['python'])
    prepare.depends_on(plink)
    prepare.command(
        f"""
set -euo pipefail
python3 {summary_script} prepare-maps \\
    --auto-bim {plink.autosomes.bim} --x-bim {plink.x.bim} \\
    --out-auto-bim {prepare.autosomes_bim} --out-x-bim {prepare.x_bim} --out-stats {prepare.map_stats}
""",
    )

    autosomal_ibis = _register_ibis_job(b, f'{job_name}_ibis_autosomes', 'ibis_autosomes', images['ibis'])
    autosomal_ibis.depends_on(prepare)
    autosomal_ibis.command(
        ibis_command(
            plink.autosomes.bed,
            prepare.autosomes_bim,
            plink.autosomes.fam,
            autosomal_ibis.ibis,
            min_seg_cm,
            min_ibd2_cm,
            config_retrieve([*CONFIG_PATH, 'ibis_autosomes', 'cpu']),
        ),
    )

    x_ibis = _register_ibis_job(b, f'{job_name}_ibis_x', 'ibis_x', images['ibis'])
    x_ibis.depends_on(prepare)
    x_ibis.command(
        ibis_command(
            plink.x.bed,
            prepare.x_bim,
            plink.x.fam,
            x_ibis.ibis,
            min_seg_cm,
            min_ibd2_cm,
            config_retrieve([*CONFIG_PATH, 'ibis_x', 'cpu']),
            chromosome='X',
        ),
    )

    run_info = {
        'inputs': {
            'pgen': pgen_path,
            'pvar': pvar_path,
            'psam': psam_path,
            'king_seg': king_seg_path,
            'remove_list': remove_list_path,
            'genetic_maps': genetic_map_paths,
        },
        'images': images,
        'plink2_chromosomes': 'autosomes: --chr 1-22; X: --chr X --remove-nosex (PAR1 and PAR2 excluded)',
        'min_ibd2_cm': min_ibd2_cm,
        'ibis_command': (
            f'ibis <bed> <bim> <fam> -min_l {min_seg_cm:g} -min_l2 {min_ibd2_cm:g} -2 -printCoef -gzip -t <threads> '
            '[-chr X for the X run]'
        ),
    }
    king_seg = b.read_input(king_seg_path)
    summarise = _register_job(b, f'{job_name}_summarise', 'summarise', images['python'])
    summarise.depends_on(autosomal_ibis, x_ibis)
    summarise.command(
        f"""
set -euo pipefail
python3 {summary_script} summarise \\
    --ibis-auto-seg {autosomal_ibis.ibis.seg} --ibis-x-seg {x_ibis.ibis.seg} \\
    --ibis-auto-log {autosomal_ibis.ibis.log} \\
    --king-seg {king_seg} \\
    --fam {plink.autosomes.fam} --fam-x {plink.x.fam} --map-stats {prepare.map_stats} \\
    --min-seg-cm {min_seg_cm:g} --join-gap-cm {join_gap_cm:g} \\
    --cutoffs-cm {','.join(f'{cutoff:g}' for cutoff in cutoffs_cm)} \\
    --run-info {quote(json.dumps(run_info))} \\
    --out-table {summarise.table} --out-sidecar {summarise.sidecar}
""",
    )

    for resource, path in (
        (autosomal_ibis.ibis.seg, output_autosomes_seg_path),
        (autosomal_ibis.ibis.coef, output_autosomes_coef_path),
        (autosomal_ibis.ibis.log, output_autosomes_log_path),
        (x_ibis.ibis.seg, output_x_seg_path),
        (x_ibis.ibis.coef, output_x_coef_path),
        (x_ibis.ibis.log, output_x_log_path),
        (summarise.table, output_pairs_path),
        (summarise.sidecar, output_sidecar_path),
    ):
        b.write_output(resource, path)

    return [plink, prepare, autosomal_ibis, x_ibis, summarise]
