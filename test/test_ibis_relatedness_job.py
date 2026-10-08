"""Tests for the IBIS relatedness jobs: command text, job wiring, and plink2 output through prepare-maps."""

from __future__ import annotations

import gzip
import json
import random
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest
from cpg_utils.config import ConfigError
from cpg_utils.config import config_retrieve as real_config_retrieve

from popgen_genotyping.jobs import ibis_relatedness_job as job_module
from popgen_genotyping.jobs.ibis_relatedness_job import (
    GENETIC_MAP_LINK,
    MAP_CHROMOSOMES,
    check_inputs_exist,
    check_parameters,
    ibis_command,
    plink_command,
    resolve_genetic_map_paths,
    run_ibis_relatedness,
)
from popgen_genotyping.scripts import ibis_relatedness_summary as summary
from popgen_genotyping.stages import ExportCohortDatasets, IbisRelatedness, KingIbdseg, RelatednessExclusions

# -- Helpers ------------------------------------------------------------------


class _Resource:
    """Stands in for a Hail Batch resource: str() is its name, attributes are sub-resources."""

    def __init__(self, name: str) -> None:
        self.label = name

    def __getattr__(self, attr: str) -> _Resource:
        if attr.startswith('__'):
            raise AttributeError(attr)
        return _Resource(f'{self.label}.{attr}')

    def __str__(self) -> str:
        return self.label


class _Job(_Resource):
    """A job that records its commands, dependencies and declared resource groups."""

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.commands: list[str] = []
        self.dependencies: list[_Job] = []
        self.resource_groups: dict[str, Any] = {}

    def command(self, text: str) -> None:
        self.commands.append(text)

    def depends_on(self, *jobs: _Job) -> None:
        self.dependencies.extend(jobs)

    def declare_resource_group(self, **groups: Any) -> None:
        self.resource_groups.update(groups)


GENETIC_MAPS = {chromosome: f'gs://ref/maps/chr{chromosome}.txt.gz' for chromosome in MAP_CHROMOSOMES}

CONFIG: dict[str, Any] = {
    'workflow': {
        'driver_image': 'driver:latest',
        'ibis_image': 'ibis:1',
        'plink2_cm_map_image': 'plink2-cm-map:1',
    },
    'popgen_genotyping': {
        'ibis_relatedness': {
            'plink2': {'cpu': 4, 'memory': 'standard', 'storage': '20Gi'},
            'prepare_maps': {'cpu': 1, 'memory': 'standard', 'storage': '10Gi'},
            'ibis_autosomes': {'cpu': 16, 'memory': 'lowmem', 'storage': '20Gi'},
            'ibis_x': {'cpu': 4, 'memory': 'lowmem', 'storage': '20Gi'},
            'summarise': {'cpu': 2, 'memory': 'highmem', 'storage': '10Gi'},
        },
    },
}

OUTPUT_PATHS = {
    'output_pairs_path': 'gs://o/pairs.tsv.gz',
    'output_sidecar_path': 'gs://o/pairs.json',
    'output_autosomes_seg_path': 'gs://o/auto.seg.gz',
    'output_autosomes_coef_path': 'gs://o/auto.coef.gz',
    'output_autosomes_log_path': 'gs://o/auto.log',
    'output_x_seg_path': 'gs://o/x.seg.gz',
    'output_x_coef_path': 'gs://o/x.coef.gz',
    'output_x_log_path': 'gs://o/x.log',
}


def _run_args(**overrides: Any) -> dict[str, Any]:
    """Arguments for run_ibis_relatedness, with the config's parameter values unless overridden."""
    args: dict[str, Any] = {
        'pgen_path': 'gs://agg/COH1.pgen',
        'pvar_path': 'gs://agg/COH1.pvar',
        'psam_path': 'gs://agg/COH1.psam',
        'king_seg_path': 'gs://king/COH1_king.seg',
        'remove_list_path': 'gs://exclusions/COH1.remove_list.txt',
        'genetic_map_paths': GENETIC_MAPS,
        'min_seg_cm': 5.0,
        'min_ibd2_cm': 5.0,
        'join_gap_cm': 0.5,
        'cutoffs_cm': [5.0, 10.0, 15.0],
        **OUTPUT_PATHS,
        'job_name': 'IbisRelatedness_cohort',
    }
    args.update(overrides)
    return args


def _queue(config: dict[str, Any] = CONFIG, **overrides: Any) -> tuple[MagicMock, list[_Job], MagicMock]:
    """Queue run_ibis_relatedness with a mocked batch; return the batch, the jobs and the register_job mock."""
    mock_batch = MagicMock()
    # Echo each staged path back so commands can be asserted on file names.
    mock_batch.read_input.side_effect = lambda path: _Resource(f'IN[{path}]')
    mock_batch.read_input_group.side_effect = lambda **paths: _Resource(f'GROUP[{sorted(paths.values())[0]}]')

    registered: list[_Job] = []

    def register_job(**kwargs: Any) -> _Job:
        job = _Job(kwargs['job_name'])
        registered.append(job)
        return job

    def config_retrieve(key: list[str], default: Any = None) -> Any:
        del default
        return real_config_retrieve(key, config=config)

    with (
        patch.object(job_module, 'get_batch', return_value=mock_batch),
        patch.object(job_module, 'get_driver_image', return_value=config['workflow']['driver_image']),
        patch.object(job_module, 'config_retrieve', side_effect=config_retrieve),
        patch.object(job_module, 'check_inputs_exist') as mock_check,
        patch.object(job_module, 'register_job', side_effect=register_job) as mock_register,
    ):
        jobs = run_ibis_relatedness(**_run_args(**overrides))
    del mock_check
    assert jobs == registered
    return mock_batch, jobs, mock_register


# -- Tests: resolve_genetic_map_paths ------------------------------------------


class TestResolveGeneticMapPaths:
    """The per-chromosome maps come from the references config."""

    def test_one_reference_key_per_chromosome(self) -> None:
        """The template is filled with 1-22 and X and each key resolved with reference_path."""
        with patch.object(job_module, 'reference_path', side_effect=lambda key: f'gs://ref/{key}') as mock_reference:
            result = resolve_genetic_map_paths('genetic_maps_decode2019/sexavg_grch38_chr{chromosome}')

        assert list(result) == [*(str(c) for c in range(1, 23)), 'X']
        assert result['1'] == 'gs://ref/genetic_maps_decode2019/sexavg_grch38_chr1'
        assert result['X'] == 'gs://ref/genetic_maps_decode2019/sexavg_grch38_chrX'
        assert mock_reference.call_count == 23

    def test_resolves_against_the_references_config_layout(self) -> None:
        """A references Source with `files` is a section of per-file keys, which reference_path walks."""
        directory = 'gs://cpg-common-main/references/genetic_maps/decode2019'
        references = {
            'genetic_maps_decode2019': {
                **{
                    f'sexavg_grch38_chr{c}': f'{directory}/decode2019_sexavg_GRCh38_chr{c}.txt.gz'
                    for c in MAP_CHROMOSOMES
                },
                'sexavg_grch38_eagle': f'{directory}/decode2019_sexavg_GRCh38.eagle.txt.gz',
            },
        }
        with patch('cpg_utils.config.get_config', return_value={'references': references}):
            result = resolve_genetic_map_paths('genetic_maps_decode2019/sexavg_grch38_chr{chromosome}')

        assert result['1'] == f'{directory}/decode2019_sexavg_GRCh38_chr1.txt.gz'
        assert result['X'] == f'{directory}/decode2019_sexavg_GRCh38_chrX.txt.gz'
        assert len(set(result.values())) == 23

    def test_missing_reference_key_raises(self) -> None:
        """A key the references config lacks is a ConfigError, not a skipped chromosome."""
        references = {'genetic_maps_decode2019': {'sexavg_grch38_chr1': 'gs://ref/chr1.txt.gz'}}
        with (
            patch('cpg_utils.config.get_config', return_value={'references': references}),
            pytest.raises(ConfigError, match='sexavg_grch38_chr2'),
        ):
            resolve_genetic_map_paths('genetic_maps_decode2019/sexavg_grch38_chr{chromosome}')

    def test_template_without_placeholder_raises(self) -> None:
        """A key without {chromosome} would give every chromosome the same map."""
        with pytest.raises(ValueError, match=r'no \{chromosome\} placeholder'):
            resolve_genetic_map_paths('genetic_maps_decode2019/sexavg_grch38_chr1')


# -- Tests: input and parameter checks -----------------------------------------


class TestChecks:
    """Checks that fail before any job is queued."""

    def test_missing_inputs_are_all_listed(self) -> None:
        """Every missing path is named, not just the first."""

        def to_path(path: str) -> MagicMock:
            return MagicMock(**{'exists.return_value': path.endswith('here')})

        with (
            patch.object(job_module, 'to_path', side_effect=to_path),
            pytest.raises(FileNotFoundError) as excinfo,
        ):
            check_inputs_exist(['gs://a/here', 'gs://a/gone1', 'gs://a/gone2'])
        assert 'gs://a/gone1' in str(excinfo.value)
        assert 'gs://a/gone2' in str(excinfo.value)
        assert 'gs://a/here' not in str(excinfo.value)

    def test_existing_inputs_pass(self) -> None:
        """No error when every path exists."""
        with patch.object(job_module, 'to_path', return_value=MagicMock(**{'exists.return_value': True})):
            check_inputs_exist(['gs://a/b'])

    def test_valid_parameters_pass(self) -> None:
        """The configured parameters are accepted."""
        check_parameters(5.0, 0.5, [5.0, 10.0, 15.0])
        check_parameters(5.0, 0.0, [5.0])

    def test_cutoff_below_min_seg_cm_raises(self) -> None:
        """IBIS never reports a run shorter than -min_l, so such a cutoff is rejected."""
        with pytest.raises(ValueError, match='below min_seg_cm'):
            check_parameters(5.0, 0.5, [2.0, 10.0])

    def test_repeated_cutoff_raises(self) -> None:
        """A repeated cutoff would give two columns the same name."""
        with pytest.raises(ValueError, match='repeats a value'):
            check_parameters(5.0, 0.5, [10.0, 10.0])

    def test_no_cutoffs_raises(self) -> None:
        """At least one cutoff is needed for a kinship column."""
        with pytest.raises(ValueError, match='at least one cutoff'):
            check_parameters(5.0, 0.5, [])

    def test_negative_join_gap_raises(self) -> None:
        """A negative join gap is meaningless."""
        with pytest.raises(ValueError, match='join_gap_cm must not be negative'):
            check_parameters(5.0, -0.1, [5.0])


# -- Tests: shell commands -----------------------------------------------------


class TestIbisCommand:
    """The IBIS invocation."""

    def test_autosomal_command(self) -> None:
        """-min_l and -min_l2 are set, IBD2 and coefficients are written gzipped, threads match the job's cpus."""
        cmd = ibis_command(_Resource('BED'), _Resource('BIM'), _Resource('FAM'), _Resource('OUT'), 5.0, 5.0, threads=16)
        assert 'ibis BED BIM FAM -min_l 5 -min_l2 5 -2 -printCoef -gzip -t 16' in cmd
        assert '-chr' not in cmd
        assert '-f OUT' in cmd
        assert 'tee OUT.log' in cmd
        assert 'set -euo pipefail' in cmd
        assert 'OUT.seg OUT.coef' in cmd

    def test_x_command_restricts_to_chromosome_x(self) -> None:
        """The X run passes -chr X."""
        cmd = ibis_command(
            _Resource('BED'), _Resource('BIM'), _Resource('FAM'), _Resource('OUT'), 5.0, 5.0, 4, chromosome='X'
        )
        assert '-t 4' in cmd
        assert '-chr X -f OUT' in cmd

    def test_fractional_lengths_keep_their_decimals(self) -> None:
        """5.5 cM is passed as 5.5, not rounded."""
        cmd = ibis_command(_Resource('B'), _Resource('M'), _Resource('F'), _Resource('O'), 5.5, 2.5, 1)
        assert '-min_l 5.5 -min_l2 2.5' in cmd


def _plink_command_for(directory: Path) -> str:
    """The plink2 job's bash script, with local paths standing in for the staged inputs and outputs."""
    pfile = _Resource(str(directory / 'agg'))
    maps = {chromosome: str(directory / f'map{chromosome}.txt.gz') for chromosome in MAP_CHROMOSOMES}
    return plink_command(
        pfile,
        str(directory / 'remove.txt'),
        maps,
        _Resource(str(directory / 'autosomes')),
        _Resource(str(directory / 'x')),
    )


class TestPlinkCommand:
    """The plink2 job's bash script."""

    def test_is_valid_bash(self, tmp_path: Path) -> None:
        """The generated script parses."""
        subprocess.run(['bash', '-n', '-c', _plink_command_for(tmp_path)], check=True)  # noqa: S603, S607

    def test_autosomes_and_x_are_separate_cm_mapped_filesets(self, tmp_path: Path) -> None:
        """Autosomes then X (unknown sex removed), each with the cM map, chromosomes as 1-22 and X."""
        cmd = _plink_command_for(tmp_path)
        autosomes = f'plink2 --pfile {tmp_path / "agg"} --remove {tmp_path / "remove.txt"} --chr 1-22 \\\n'
        x = f'plink2 --pfile {tmp_path / "agg"} --remove {tmp_path / "remove.txt"} --remove-nosex --chr X \\\n'
        assert autosomes in cmd
        assert x in cmd
        assert cmd.count('--cm-map "$cm_map" --output-chr MT --make-bed') == 2

    def test_each_chromosome_map_is_linked_under_the_at_pattern_name(self, tmp_path: Path) -> None:
        """One link per chromosome, named so that plink2's '@' finds it."""
        cmd = _plink_command_for(tmp_path)
        for chromosome in MAP_CHROMOSOMES:
            link = GENETIC_MAP_LINK.format(chromosome=chromosome)
            assert f'ln -s {tmp_path / f"map{chromosome}.txt.gz"} maps/{link}' in cmd
        assert f"cm_map='maps/{GENETIC_MAP_LINK.format(chromosome='@')}'" in cmd

    def test_remove_list_header_and_sample_counts_are_checked(self, tmp_path: Path) -> None:
        """A remove list without the #IID header, or a wrong sample count after --remove, fails the job."""
        cmd = _plink_command_for(tmp_path)
        assert '!= "#IID"' in cmd
        assert 'expected $((n_input - n_excluded)) samples after --remove' in cmd
        assert 'implausible X sample count' in cmd

    def test_remove_list_header_check_rejects_a_list_without_header(self, tmp_path: Path) -> None:
        """Run just the guard: a list whose first line is not #IID exits 1 before plink2 is called."""
        cmd = _plink_command_for(tmp_path)
        guard = cmd.split('mkdir maps')[0]
        (tmp_path / 'remove.txt').write_text('CPG1\n')
        result = subprocess.run(['bash', '-c', guard], check=False, capture_output=True)  # noqa: S603, S607
        assert result.returncode == 1
        assert b'unexpected remove-list header' in result.stderr
        (tmp_path / 'remove.txt').write_text('#IID\nCPG1\n')
        subprocess.run(['bash', '-c', guard], check=True)  # noqa: S603, S607


# -- Tests: run_ibis_relatedness wiring ----------------------------------------


class TestRunIbisRelatedness:
    """Job order, images, resources, commands and outputs."""

    def test_five_jobs_in_dependency_order(self) -> None:
        """plink2, then prepare-maps, then IBIS autosomes and X in parallel, then summarise."""
        _batch, jobs, _register = _queue()
        plink, prepare, autosomal, x, summarise = jobs
        assert [job.label for job in jobs] == [
            'IbisRelatedness_cohort_plink2',
            'IbisRelatedness_cohort_prepare_maps',
            'IbisRelatedness_cohort_ibis_autosomes',
            'IbisRelatedness_cohort_ibis_x',
            'IbisRelatedness_cohort_summarise',
        ]
        assert plink.dependencies == []
        assert prepare.dependencies == [plink]
        assert autosomal.dependencies == [prepare]
        assert x.dependencies == [prepare]
        assert summarise.dependencies == [autosomal, x]

    def test_images_and_config_sections(self) -> None:
        """The cm-map plink2 and IBIS images come from workflow config; Python steps run in the driver image."""
        _batch, _jobs, mock_register = _queue()
        by_name = {call.kwargs['job_name']: call.kwargs for call in mock_register.call_args_list}
        sections = {
            'plink2': 'plink2-cm-map:1',
            'prepare_maps': 'driver:latest',
            'ibis_autosomes': 'ibis:1',
            'ibis_x': 'ibis:1',
            'summarise': 'driver:latest',
        }
        for section, image in sections.items():
            kwargs = by_name[f'IbisRelatedness_cohort_{section}']
            assert kwargs['image'] == image
            assert kwargs['config_path'] == ['popgen_genotyping', 'ibis_relatedness', section]

    def test_ibis_threads_equal_the_job_cpus(self) -> None:
        """-t is the configured cpu count of each IBIS job, and only the X job is restricted to X."""
        _batch, jobs, _register = _queue()
        autosomal, x = jobs[2], jobs[3]
        assert '-min_l 5 -min_l2 5 -2 -printCoef -gzip -t 16' in autosomal.commands[0]
        assert '-chr' not in autosomal.commands[0]
        assert '-t 4' in x.commands[0]
        assert '-chr X' in x.commands[0]

    def test_ibis_reads_the_corrected_bim_beside_the_plink_bed_and_fam(self) -> None:
        """IBIS gets plink2's .bed and .fam but prepare-maps' shifted .bim."""
        _batch, jobs, _register = _queue()
        plink, prepare, autosomal, x, _summarise = jobs
        assert (
            f'ibis {plink.label}.autosomes.bed {prepare.label}.autosomes_bim {plink.label}.autosomes.fam'
            in autosomal.commands[0]
        )
        assert f'ibis {plink.label}.x.bed {prepare.label}.x_bim {plink.label}.x.fam' in x.commands[0]

    def test_each_ibis_job_declares_seg_coef_and_log(self) -> None:
        """The three outputs of an IBIS run are one resource group."""
        _batch, jobs, _register = _queue()
        for job in (jobs[2], jobs[3]):
            assert job.resource_groups['ibis'] == {
                'seg': '{root}.seg.gz',
                'coef': '{root}.coef.gz',
                'log': '{root}.log',
            }

    def test_plink_job_gets_the_remove_list_and_maps(self) -> None:
        """The staged remove list and one map per chromosome reach the plink2 command."""
        _batch, jobs, _register = _queue()
        cmd = jobs[0].commands[0]
        assert '--remove IN[gs://exclusions/COH1.remove_list.txt]' in cmd
        for chromosome in MAP_CHROMOSOMES:
            assert f'ln -s IN[gs://ref/maps/chr{chromosome}.txt.gz] maps/' in cmd

    def test_summarise_command(self) -> None:
        """The summarise step gets the KING .seg, both fams, the parameters and a JSON run description."""
        _batch, jobs, _register = _queue()
        plink, prepare, autosomal, x, summarise = jobs
        cmd = summarise.commands[0]
        assert 'ibis_relatedness_summary.py' in cmd
        assert 'summarise' in cmd
        assert f'--ibis-auto-seg {autosomal.label}.ibis.seg --ibis-x-seg {x.label}.ibis.seg' in cmd
        assert f'--ibis-auto-log {autosomal.label}.ibis.log' in cmd
        assert '--king-seg IN[gs://king/COH1_king.seg]' in cmd
        assert f'--fam {plink.label}.autosomes.fam --fam-x {plink.label}.x.fam' in cmd
        assert f'--map-stats {prepare.label}.map_stats' in cmd
        assert '--min-seg-cm 5 --join-gap-cm 0.5' in cmd
        assert '--cutoffs-cm 5,10,15' in cmd
        run_info = json.loads(shlex.split(cmd.split('--run-info ')[1])[0])
        assert run_info['inputs']['king_seg'] == 'gs://king/COH1_king.seg'
        assert run_info['inputs']['remove_list'] == 'gs://exclusions/COH1.remove_list.txt'
        assert run_info['inputs']['genetic_maps'] == GENETIC_MAPS
        assert run_info['images'] == {'ibis': 'ibis:1', 'plink2': 'plink2-cm-map:1', 'python': 'driver:latest'}
        assert run_info['min_ibd2_cm'] == 5.0

    def test_prepare_maps_command(self) -> None:
        """prepare-maps reads both .bim files plink2 wrote and writes the corrected ones and the map stats."""
        _batch, jobs, _register = _queue()
        plink, prepare = jobs[0], jobs[1]
        cmd = prepare.commands[0]
        assert 'ibis_relatedness_summary.py' in cmd
        assert 'prepare-maps' in cmd
        assert f'--auto-bim {plink.label}.autosomes.bim --x-bim {plink.label}.x.bim' in cmd

    def test_all_eight_outputs_are_written_to_their_cloud_paths(self) -> None:
        """Both IBIS runs (seg, coef, log), the pairs table and its sidecar are written out."""
        batch, jobs, _register = _queue()
        _plink, _prepare, autosomal, x, summarise = jobs
        written = {(str(call.args[0]), call.args[1]) for call in batch.write_output.call_args_list}
        assert written == {
            (f'{autosomal.label}.ibis.seg', 'gs://o/auto.seg.gz'),
            (f'{autosomal.label}.ibis.coef', 'gs://o/auto.coef.gz'),
            (f'{autosomal.label}.ibis.log', 'gs://o/auto.log'),
            (f'{x.label}.ibis.seg', 'gs://o/x.seg.gz'),
            (f'{x.label}.ibis.coef', 'gs://o/x.coef.gz'),
            (f'{x.label}.ibis.log', 'gs://o/x.log'),
            (f'{summarise.label}.table', 'gs://o/pairs.tsv.gz'),
            (f'{summarise.label}.sidecar', 'gs://o/pairs.json'),
        }

    def test_missing_resource_key_raises_instead_of_using_a_default(self) -> None:
        """A job section without cpu/memory/storage raises: resources have no code defaults."""
        for missing in ('cpu', 'memory', 'storage'):
            config = json.loads(json.dumps(CONFIG))
            del config['popgen_genotyping']['ibis_relatedness']['ibis_x'][missing]
            with pytest.raises(ConfigError, match=missing):
                _queue(config=config)

    def test_missing_job_section_raises(self) -> None:
        """A missing job section raises."""
        config = json.loads(json.dumps(CONFIG))
        del config['popgen_genotyping']['ibis_relatedness']['summarise']
        with pytest.raises(ConfigError, match='summarise'):
            _queue(config=config)

    def test_invalid_parameters_raise_before_any_job_is_queued(self) -> None:
        """A bad cutoff stops the run at submission rather than in the summarise job."""
        with (
            patch.object(job_module, 'get_batch') as mock_get_batch,
            patch.object(job_module, 'register_job') as mock_register,
            pytest.raises(ValueError, match='below min_seg_cm'),
        ):
            run_ibis_relatedness(**_run_args(cutoffs_cm=[2.0]))
        mock_register.assert_not_called()
        mock_get_batch.assert_not_called()

    def test_missing_genetic_map_raises_before_any_job_is_queued(self) -> None:
        """A map that does not exist stops the run at submission."""
        with (
            patch.object(job_module, 'to_path', return_value=MagicMock(**{'exists.return_value': False})),
            patch.object(job_module, 'get_batch') as mock_get_batch,
            patch.object(job_module, 'register_job') as mock_register,
            pytest.raises(FileNotFoundError, match=r'gs://ref/maps/chrX\.txt\.gz'),
        ):
            run_ibis_relatedness(**_run_args())
        mock_register.assert_not_called()
        mock_get_batch.assert_not_called()

    def test_genetic_maps_must_cover_exactly_chromosomes_1_to_22_and_x(self) -> None:
        """A map set missing a chromosome, or with an extra one, is rejected."""
        short = {k: v for k, v in GENETIC_MAPS.items() if k != 'X'}
        with pytest.raises(ValueError, match='expected'):
            run_ibis_relatedness(**_run_args(genetic_map_paths=short))
        with pytest.raises(ValueError, match='expected'):
            run_ibis_relatedness(**_run_args(genetic_map_paths={**GENETIC_MAPS, 'Y': 'gs://ref/maps/chrY.txt.gz'}))


# -- Tests: plink2 output through prepare-maps ---------------------------------


def _plink2_reads_cm_maps() -> bool:
    """Whether the installed plink2 has --cm-map (the plink image's alpha build does not)."""
    if shutil.which('plink2') is None:
        return False
    result = subprocess.run(['plink2', '--help', 'cm-map'], check=False, capture_output=True, text=True)  # noqa: S607
    return '--cm-map' in result.stdout


def _write_aggregate_and_maps(directory: Path) -> None:
    """Write a 6-sample chr-prefixed PGEN (chr1-22 and chrX, one unknown-sex sample) and a map per chromosome.

    Each map starts after the first variants, so plink2 extrapolates them to negative cM, as the real
    deCODE maps do for the first markers of a chromosome.
    """
    rng = random.Random(1)  # noqa: S311 - seeded synthetic genotypes, not a security use
    samples = [f'S{i}' for i in range(6)]
    contigs = [f'chr{c}' for c in MAP_CHROMOSOMES]
    lines = ['##fileformat=VCFv4.2', *(f'##contig=<ID={contig}>' for contig in contigs)]
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    lines.append('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', *samples]))
    for contig in contigs:
        for i in range(8):
            # Outside the hg38 pseudoautosomal regions, which plink2 would insist on splitting.
            position = 10_000_000 + i * 1_000_000
            genotypes = [f'{rng.randint(0, 1)}/{rng.randint(0, 1)}' for _ in samples]
            lines.append('\t'.join([contig, str(position), f'{contig}_{i}', 'A', 'G', '.', '.', '.', 'GT', *genotypes]))
    (directory / 'agg.vcf').write_text('\n'.join(lines) + '\n')
    # S5 has unknown sex, so the X fileset has one sample fewer than the autosomal one.
    (directory / 'sex.txt').write_text('#IID\tSEX\n' + ''.join(f'S{i}\t{[1, 2, 1, 2, 1, 0][i]}\n' for i in range(6)))
    subprocess.run(  # noqa: S603
        [  # noqa: S607
            'plink2',
            '--vcf',
            str(directory / 'agg.vcf'),
            '--update-sex',
            str(directory / 'sex.txt'),
            '--output-chr',
            'chrM',
            '--make-pgen',
            '--out',
            str(directory / 'agg'),
        ],
        check=True,
        capture_output=True,
    )
    for chromosome in MAP_CHROMOSOMES:
        with gzip.open(directory / f'map{chromosome}.txt.gz', 'wt') as f:
            f.write('pposition rrate gposition\n')
            # 1 cM per Mb from 13 Mb, so the first three variants lie before the map.
            for megabase in range(13, 40):
                f.write(f'{megabase * 1_000_000} 1.0 {float(megabase - 13)}\n')


@pytest.mark.skipif(not _plink2_reads_cm_maps(), reason='needs a plink2 with --cm-map (2.0 beta 1 or later)')
class TestPlinkCommandWithPlink2:
    """Run the plink2 job's script with a real plink2 and check what prepare-maps needs from it."""

    def _run(self, tmp_path: Path, remove: list[str]) -> None:
        _write_aggregate_and_maps(tmp_path)
        (tmp_path / 'remove.txt').write_text('#IID\n' + ''.join(f'{iid}\n' for iid in remove))
        subprocess.run(  # noqa: S603
            ['bash', '-c', _plink_command_for(tmp_path)],  # noqa: S607
            check=True,
            capture_output=True,
            cwd=tmp_path,
        )

    def test_filesets_pass_prepare_maps(self, tmp_path: Path) -> None:
        """chr-prefixed input comes out as 1-22 and X with a cM map, which prepare-maps accepts and shifts."""
        self._run(tmp_path, remove=['S4'])

        autosomal = summary.read_bim(tmp_path / 'autosomes.bim')
        x_bim = summary.read_bim(tmp_path / 'x.bim')
        summary.check_chromosome_blocks(autosomal, summary.AUTOSOMES)
        summary.check_chromosome_blocks(x_bim, [summary.X_CHROMOSOME])
        # plink2 extrapolated the variants before each map's first row to negative cM.
        assert (autosomal['CM'] < 0).any()
        assert (x_bim['CM'] < 0).any()

        summary.main(
            [
                'prepare-maps',
                '--auto-bim',
                str(tmp_path / 'autosomes.bim'),
                '--x-bim',
                str(tmp_path / 'x.bim'),
                '--out-auto-bim',
                str(tmp_path / 'autosomes_fixed.bim'),
                '--out-x-bim',
                str(tmp_path / 'x_fixed.bim'),
                '--out-stats',
                str(tmp_path / 'stats.json'),
            ],
        )
        stats = json.loads((tmp_path / 'stats.json').read_text())
        # The shift keeps every distance between markers: 5 variants per chromosome span 5 Mb = 5 cM at
        # 1 cM per Mb, once past the map start.
        assert stats['autosomal_map_length_cm'] == pytest.approx(22 * 7.0, rel=1e-3)
        assert stats['x_map_length_cm'] == pytest.approx(7.0, rel=1e-3)
        assert stats['n_autosomal_variants'] == 22 * 8
        assert stats['n_x_variants'] == 8

    def test_removed_samples_and_unknown_sex_are_dropped(self, tmp_path: Path) -> None:
        """The remove list applies to both filesets, and the X fileset also drops the unknown-sex sample."""
        self._run(tmp_path, remove=['S4'])
        autosomal_samples = [line.split()[1] for line in (tmp_path / 'autosomes.fam').read_text().splitlines()]
        x_samples = [line.split()[1] for line in (tmp_path / 'x.fam').read_text().splitlines()]
        assert autosomal_samples == ['S0', 'S1', 'S2', 'S3', 'S5']
        assert x_samples == ['S0', 'S1', 'S2', 'S3']

    def test_header_only_remove_list_removes_nothing(self, tmp_path: Path) -> None:
        """An empty exclusion (header-only list) leaves every sample in the autosomal fileset."""
        self._run(tmp_path, remove=[])
        assert len((tmp_path / 'autosomes.fam').read_text().splitlines()) == 6

    def test_count_check_fails_when_a_remove_list_sample_is_not_in_the_aggregate(self, tmp_path: Path) -> None:
        """A removed ID that is not in the aggregate makes the sample count differ, which fails the job."""
        with pytest.raises(subprocess.CalledProcessError):
            self._run(tmp_path, remove=['NOT_A_SAMPLE'])


# -- Tests: IbisRelatedness stage ---------------------------------------------


class TestIbisRelatednessExpectedOutputs:
    """Path layout that Metamist `analysis_keys=['pairs']` depends on."""

    def test_path_layout_and_naming(self) -> None:
        """Eight outputs under the standard prefix, cohort-keyed so successive aggregates never collide."""
        prefix = Path('/cohort/ibis')
        mock_cohort = MagicMock()
        mock_cohort.id = 'COH123'
        mock_self = MagicMock()
        mock_self.name = 'IbisRelatedness'

        with patch('popgen_genotyping.stages.get_output_prefix', return_value=prefix) as mock_prefix:
            result = IbisRelatedness.expected_outputs(mock_self, mock_cohort)

        assert result == {
            'pairs': prefix / 'COH123_ibis_relatedness_pairs.tsv.gz',
            'sidecar': prefix / 'COH123_ibis_relatedness_pairs.json',
            'autosomes_seg': prefix / 'COH123_ibis_autosomes.seg.gz',
            'autosomes_coef': prefix / 'COH123_ibis_autosomes.coef.gz',
            'autosomes_log': prefix / 'COH123_ibis_autosomes.log',
            'x_seg': prefix / 'COH123_ibis_x.seg.gz',
            'x_coef': prefix / 'COH123_ibis_x.coef.gz',
            'x_log': prefix / 'COH123_ibis_x.log',
        }
        mock_prefix.assert_called_once_with(dataset=mock_cohort.dataset, stage_name='IbisRelatedness')


class TestIbisRelatednessQueueJobs:
    """Wire-up between upstream outputs, config, the references and the job factory."""

    def test_passes_inputs_config_and_resolved_maps_to_job(self) -> None:
        """Aggregate PGEN, remove list, KING .seg, parameters from config and the resolved maps reach the job."""
        mock_cohort = MagicMock()
        mock_cohort.name = 'my_cohort'

        def as_path(target: object, stage: Any, key: str) -> Path:
            del target
            return {
                (ExportCohortDatasets, 'pgen'): Path('/agg/COH1.pgen'),
                (ExportCohortDatasets, 'pvar'): Path('/agg/COH1.pvar'),
                (ExportCohortDatasets, 'psam'): Path('/agg/COH1.psam'),
                (RelatednessExclusions, 'remove_list'): Path('/exclusions/COH1.remove_list.txt'),
                (KingIbdseg, 'seg'): Path('/king/COH1_king.seg'),
            }[(stage, key)]

        mock_inputs = MagicMock()
        mock_inputs.as_path.side_effect = as_path

        expected_outputs = {
            key: Path(f'/out/{key}')
            for key in (
                'pairs',
                'sidecar',
                'autosomes_seg',
                'autosomes_coef',
                'autosomes_log',
                'x_seg',
                'x_coef',
                'x_log',
            )
        }
        mock_self = MagicMock()
        mock_self.expected_outputs.return_value = expected_outputs

        # TOML integers: the stage must turn them into floats for the job.
        values = {
            'genetic_map_reference_key': 'genetic_maps_decode2019/sexavg_grch38_chr{chromosome}',
            'min_seg_cm': 5,
            'min_ibd2_cm': 5,
            'join_gap_cm': 0.5,
            'cutoffs_cm': [5, 10, 15],
        }

        def config_retrieve(path: list[str]) -> Any:
            # No default is passed: a missing key must raise, so each lookup is a one-argument call.
            assert path[:2] == ['popgen_genotyping', 'ibis_relatedness']
            return values[path[2]]

        with (
            patch('popgen_genotyping.stages.config_retrieve', side_effect=config_retrieve),
            patch('popgen_genotyping.stages.resolve_genetic_map_paths', return_value=GENETIC_MAPS) as mock_maps,
            patch('popgen_genotyping.stages.run_ibis_relatedness') as mock_run,
        ):
            IbisRelatedness.queue_jobs(mock_self, mock_cohort, mock_inputs)

        mock_maps.assert_called_once_with('genetic_maps_decode2019/sexavg_grch38_chr{chromosome}')
        mock_run.assert_called_once_with(
            pgen_path='/agg/COH1.pgen',
            pvar_path='/agg/COH1.pvar',
            psam_path='/agg/COH1.psam',
            king_seg_path='/king/COH1_king.seg',
            remove_list_path='/exclusions/COH1.remove_list.txt',
            genetic_map_paths=GENETIC_MAPS,
            min_seg_cm=5.0,
            min_ibd2_cm=5.0,
            join_gap_cm=0.5,
            cutoffs_cm=[5.0, 10.0, 15.0],
            output_pairs_path='/out/pairs',
            output_sidecar_path='/out/sidecar',
            output_autosomes_seg_path='/out/autosomes_seg',
            output_autosomes_coef_path='/out/autosomes_coef',
            output_autosomes_log_path='/out/autosomes_log',
            output_x_seg_path='/out/x_seg',
            output_x_coef_path='/out/x_coef',
            output_x_log_path='/out/x_log',
            job_name='IbisRelatedness_my_cohort',
        )
        assert all(isinstance(cutoff, float) for cutoff in mock_run.call_args.kwargs['cutoffs_cm'])
        mock_self.make_outputs.assert_called_once_with(mock_cohort, data=expected_outputs, jobs=mock_run.return_value)
