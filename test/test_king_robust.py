"""Tests for KingRobust: job command, stage wiring, and a real plink2 run when plink2 is installed."""

from __future__ import annotations

import gzip
import random
import shutil
import subprocess
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest

from popgen_genotyping.jobs.king_robust_job import run_king_robust
from popgen_genotyping.scripts.merge_qc import KING_ROBUST_RELATIVES_COL, process_king_robust
from popgen_genotyping.stages import ExportCohortDatasets, KingRobust, RelatednessExclusions

# -- Helpers ------------------------------------------------------------------


class _ResourceGroup:
    """Stands in for a Hail Batch resource group: str() is the root, attributes are the files."""

    def __init__(self, root: Path) -> None:
        self.root = str(root)
        self.kin0 = f'{root}.kin0'
        self.kin0_gz = f'{root}.kin0.gz'
        self.log = f'{root}.log'

    def __str__(self) -> str:
        return self.root


def _queue_job(
    tmp_path: Path, pfile_prefix: str = 'PFILE', remove_list: str = 'REMOVE'
) -> tuple[MagicMock, MagicMock, MagicMock, _ResourceGroup]:
    """Queue run_king_robust with mocked Hail Batch.

    Args:
        tmp_path: Root for the simulated job output files.
        pfile_prefix: What the staged pgen/pvar/psam input group renders as in the command.
        remove_list: What the staged remove list renders as in the command.

    Returns:
        The mocked batch, the job, the register_job mock, and the job's output resource group.
    """
    mock_batch = MagicMock()
    mock_batch.read_input_group.return_value = pfile_prefix
    mock_batch.read_input.return_value = remove_list
    outputs = _ResourceGroup(tmp_path / 'king_robust')
    mock_job = MagicMock()
    mock_job.king_robust = outputs

    with (
        patch('popgen_genotyping.jobs.king_robust_job.get_batch', return_value=mock_batch),
        patch('popgen_genotyping.jobs.king_robust_job.config_retrieve', return_value='plink-image:1.0'),
        patch('popgen_genotyping.jobs.king_robust_job.register_job', return_value=mock_job) as mock_register,
    ):
        returned = run_king_robust(
            pgen_path='gs://x/agg.pgen',
            pvar_path='gs://x/agg.pvar',
            psam_path='gs://x/agg.psam',
            remove_list_path='gs://x/remove.txt',
            kinship_min=0.0442,
            output_kin0_path='gs://o/out.kin0.gz',
            output_log_path='gs://o/out.log',
            job_name='KingRobust_cohort',
        )
    assert returned is mock_job
    return mock_batch, mock_job, mock_register, outputs


# -- Tests: run_king_robust ---------------------------------------------------


class TestRunKingRobust:
    """The plink2 command and its wiring."""

    def test_inputs_are_staged_from_the_given_cloud_paths(self, tmp_path: Path) -> None:
        """pgen/pvar/psam are staged as one input group and the remove list as one file."""
        batch, _job, _register, _outputs = _queue_job(tmp_path)
        batch.read_input_group.assert_called_once_with(
            pgen='gs://x/agg.pgen', pvar='gs://x/agg.pvar', psam='gs://x/agg.psam'
        )
        batch.read_input.assert_called_once_with('gs://x/remove.txt')

    def test_job_uses_the_plink_image_and_its_config_section(self, tmp_path: Path) -> None:
        """The existing workflow plink image runs it, with resources from [popgen_genotyping.king_robust]."""
        _batch, _job, mock_register, _outputs = _queue_job(tmp_path)
        kwargs = mock_register.call_args.kwargs
        assert kwargs['config_path'] == ['popgen_genotyping', 'king_robust']
        assert kwargs['image'] == 'plink-image:1.0'
        assert kwargs['job_name'] == 'KingRobust_cohort'

    def test_command_flags(self, tmp_path: Path) -> None:
        """All variants, autosomes only, the remove list, the kinship floor, and a fixed column set."""
        _batch, job, _register, outputs = _queue_job(tmp_path, pfile_prefix='PFILE', remove_list='REMOVE')
        cmd: str = job.command.call_args[0][0]
        assert 'plink2 --pfile PFILE' in cmd
        assert '--remove REMOVE' in cmd
        assert '--autosome' in cmd
        assert '--make-king-table cols=id,nsnp,hethet,ibs0,kinship' in cmd
        assert '--king-table-filter 0.0442' in cmd
        assert f'--out {outputs.root}' in cmd
        # No SNP-QC filtering: nothing restricts the variants beyond --autosome.
        for flag in ('--extract', '--exclude', '--maf', '--geno', '--hwe'):
            assert flag not in cmd

    def test_job_fails_when_plink2_wrote_nothing(self, tmp_path: Path) -> None:
        """A missing or empty .kin0 or log fails the job rather than writing out an empty table."""
        _batch, job, _register, outputs = _queue_job(tmp_path)
        cmd: str = job.command.call_args[0][0]
        # Replace plink2 with a no-op: no .kin0 and no log exist, so the check must exit non-zero.
        script = cmd.split('plink2 --pfile')[0] + cmd[cmd.index('for output in') :]
        result = subprocess.run(['bash', '-c', script], check=False, capture_output=True)  # noqa: S603, S607
        assert result.returncode == 1
        assert b'plink2 wrote no' in result.stderr
        assert not Path(outputs.kin0_gz).exists()

    def test_every_output_file_is_named_in_the_command(self, tmp_path: Path) -> None:
        """Hail Batch only writes out resources the command mentions, so each written file must appear."""
        _batch, job, _register, outputs = _queue_job(tmp_path)
        cmd: str = job.command.call_args[0][0]
        assert outputs.kin0_gz in cmd
        assert outputs.log in cmd

    def test_table_is_gzipped_for_python_to_read(self, tmp_path: Path) -> None:
        """plink2 writes the plain table; the job gzips it (zstd would need an extra Python package)."""
        _batch, job, _register, outputs = _queue_job(tmp_path)
        cmd: str = job.command.call_args[0][0]
        assert f'gzip -c {outputs.kin0} > {outputs.kin0_gz}' in cmd
        assert ' zs' not in cmd

    def test_gzipped_table_and_log_are_written_out(self, tmp_path: Path) -> None:
        """The gzipped table and the log go to their cloud paths; the plain intermediate does not."""
        batch, _job, _register, outputs = _queue_job(tmp_path)
        written = {call.args for call in batch.write_output.call_args_list}
        assert written == {(outputs.kin0_gz, 'gs://o/out.kin0.gz'), (outputs.log, 'gs://o/out.log')}


# -- Tests: the command against a real plink2 ----------------------------------

N_SAMPLES = 20
# Enough variants that unrelated pairs sit well below the 0.0442 floor: the kinship estimate's noise
# shrinks with the variant count.
N_VARIANTS_PER_CONTIG = 2000


def _degrees(relatives_cell: str) -> dict[str, str]:
    """Relative ID -> degree from a ``REL_ID:KINSHIP:DEGREE;...`` cell."""
    return {rel_id: degree for rel_id, _kinship, degree in (item.split(':') for item in relatives_cell.split(';'))}


def _write_synthetic_pgen(directory: Path) -> str:
    """Write a 20-sample chr-prefixed PGEN fileset with a duplicate pair and a parent-offspring-like pair.

    S1 duplicates S0 and S3 shares one haplotype with S2. Every other sample is random, so unrelated
    pairs sit near kinship 0 and are filtered out by the 0.0442 floor. Variants are on chr1 and chr2,
    plus 200 on chrX that the autosomes-only command must leave out.

    Args:
        directory: Directory to write the VCF and the PGEN fileset into.

    Returns:
        The PGEN fileset prefix.
    """
    rng = random.Random(1)  # noqa: S311 - seeded synthetic genotypes, not a security use
    samples = [f'S{i}' for i in range(N_SAMPLES)]
    contigs = [('chr1', N_VARIANTS_PER_CONTIG), ('chr2', N_VARIANTS_PER_CONTIG), ('chrX', 200)]
    lines = ['##fileformat=VCFv4.2', *(f'##contig=<ID={contig}>' for contig, _ in contigs)]
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    lines.append('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', *samples]))
    for contig, n_variants in contigs:
        # chrX sits well clear of the hg38 pseudoautosomal regions, which plink2 would insist on splitting.
        offset = 10_000_000 if contig == 'chrX' else 0
        for i in range(1, n_variants + 1):
            freq = rng.uniform(0.1, 0.9)
            haplotypes = {s: (int(rng.random() < freq), int(rng.random() < freq)) for s in samples}
            haplotypes['S1'] = haplotypes['S0']
            haplotypes['S3'] = (haplotypes['S2'][0], int(rng.random() < freq))
            genotypes = [f'{haplotypes[s][0]}/{haplotypes[s][1]}' for s in samples]
            lines.append(
                '\t'.join([contig, str(offset + i * 1000), f'{contig}_{i}', 'A', 'G', '.', '.', '.', 'GT', *genotypes])
            )
    vcf = directory / 'synthetic.vcf'
    vcf.write_text('\n'.join(lines) + '\n')
    sex = directory / 'sex.txt'
    sex.write_text('#IID\tSEX\n' + ''.join(f'S{i}\t{2 if i % 2 else 1}\n' for i in range(N_SAMPLES)))
    prefix = directory / 'agg'
    subprocess.run(  # noqa: S603
        [  # noqa: S607
            'plink2',
            '--vcf',
            str(vcf),
            '--update-sex',
            str(sex),
            '--output-chr',
            'chrM',
            '--make-pgen',
            '--out',
            str(prefix),
        ],
        check=True,
        capture_output=True,
    )
    return str(prefix)


@pytest.mark.skipif(shutil.which('plink2') is None, reason='plink2 is not installed')
class TestKingRobustWithPlink2:
    """Run the job's command with a real plink2 and parse its output as QcReport does."""

    def _run(self, tmp_path: Path, remove: list[str]) -> Path:
        """Execute the job's bash command on a synthetic fileset; return the gzipped .kin0 path."""
        pfile_prefix = _write_synthetic_pgen(tmp_path)
        remove_list = tmp_path / 'remove.txt'
        remove_list.write_text('#IID\n' + ''.join(f'{iid}\n' for iid in remove))
        _batch, job, _register, outputs = _queue_job(tmp_path, pfile_prefix=pfile_prefix, remove_list=str(remove_list))
        subprocess.run(['bash', '-c', job.command.call_args[0][0]], check=True, capture_output=True)  # noqa: S603, S607
        return Path(outputs.kin0_gz)

    def test_duplicate_and_close_pair_are_reported_with_degrees(self, tmp_path: Path) -> None:
        """The duplicate pair is Dup/MZ, the shared-haplotype pair 1st degree, and unrelated pairs are filtered."""
        kin0_gz = self._run(tmp_path, remove=['S19'])

        with gzip.open(kin0_gz, 'rt') as f:
            header = f.readline().rstrip('\n').split('\t')
        assert header == ['#IID1', 'IID2', 'NSNP', 'HETHET', 'IBS0', 'KINSHIP']

        relatives = process_king_robust(str(kin0_gz)).set_index('IID')[KING_ROBUST_RELATIVES_COL]
        assert _degrees(relatives['S0']) == {'S1': 'Dup/MZ'}
        assert _degrees(relatives['S1']) == {'S0': 'Dup/MZ'}
        assert _degrees(relatives['S2']) == {'S3': '1st'}
        assert _degrees(relatives['S3']) == {'S2': '1st'}
        assert set(relatives.index) == {'S0', 'S1', 'S2', 'S3'}

    def test_only_autosomal_variants_are_used(self, tmp_path: Path) -> None:
        """NSNP counts the chr1 and chr2 variants, not the chrX ones."""
        kin0_gz = self._run(tmp_path, remove=[])
        with gzip.open(kin0_gz, 'rt') as f:
            header = f.readline().rstrip('\n').split('\t')
            rows = [dict(zip(header, line.split('\t'), strict=True)) for line in f]
        assert rows
        assert {int(row['NSNP']) for row in rows} == {2 * N_VARIANTS_PER_CONTIG}

    def test_removed_samples_are_absent(self, tmp_path: Path) -> None:
        """A sample in the remove list never appears in the table, even as part of a duplicate pair."""
        kin0_gz = self._run(tmp_path, remove=['S1'])
        relatives = process_king_robust(str(kin0_gz))
        assert 'S1' not in set(relatives['IID'])
        assert 'S0' not in set(relatives['IID'])  # its only relative was S1
        assert 'S1:' not in ';'.join(relatives[KING_ROBUST_RELATIVES_COL])

    def test_header_only_table_when_no_pair_reaches_the_floor(self, tmp_path: Path) -> None:
        """Removing one member of each related pair leaves a header-only table, which still parses."""
        kin0_gz = self._run(tmp_path, remove=['S0', 'S2'])
        relatives = process_king_robust(str(kin0_gz))
        assert relatives.empty


# -- Tests: KingRobust.expected_outputs ----------------------------------------


class TestKingRobustExpectedOutputs:
    """Path layout that Metamist `analysis_keys=['kin0']` depends on."""

    def test_path_layout_and_naming(self) -> None:
        """A gzipped table and a log under the standard prefix, cohort-keyed."""
        prefix = Path('/cohort/king_robust')
        mock_cohort = MagicMock()
        mock_cohort.id = 'COH123'
        mock_self = MagicMock()
        mock_self.name = 'KingRobust'

        with patch('popgen_genotyping.stages.get_output_prefix', return_value=prefix) as mock_prefix:
            result = KingRobust.expected_outputs(mock_self, mock_cohort)

        assert result == {
            'kin0': prefix / 'COH123_king_robust.kin0.gz',
            'log': prefix / 'COH123_king_robust.log',
        }
        mock_prefix.assert_called_once_with(dataset=mock_cohort.dataset, stage_name='KingRobust')


# -- Tests: KingRobust.queue_jobs ---------------------------------------------


class TestKingRobustQueueJobs:
    """Wire-up between ExportCohortDatasets/RelatednessExclusions inputs, config, and the job factory."""

    def test_passes_pgen_remove_list_and_kinship_min_to_job(self) -> None:
        """The aggregate pgen/pvar/psam, the remove list and the configured floor reach run_king_robust."""
        mock_cohort = MagicMock()
        mock_cohort.name = 'my_cohort'

        def as_path(target: object, stage: Any, key: str) -> Path:
            del target
            return {
                (ExportCohortDatasets, 'pgen'): Path('/agg/COH1.pgen'),
                (ExportCohortDatasets, 'pvar'): Path('/agg/COH1.pvar'),
                (ExportCohortDatasets, 'psam'): Path('/agg/COH1.psam'),
                (RelatednessExclusions, 'remove_list'): Path('/exclusions/COH1.remove_list.txt'),
            }[(stage, key)]

        mock_inputs = MagicMock()
        mock_inputs.as_path.side_effect = as_path

        expected_outputs = {'kin0': Path('/out/COH1_king_robust.kin0.gz'), 'log': Path('/out/COH1_king_robust.log')}
        mock_self = MagicMock()
        mock_self.expected_outputs.return_value = expected_outputs

        def config_retrieve(path: list[str]) -> float:
            # No default is passed: a missing key must raise, so the lookup is a one-argument call.
            assert path == ['popgen_genotyping', 'king_robust', 'kinship_min']
            return 0.0442

        with (
            patch('popgen_genotyping.stages.config_retrieve', side_effect=config_retrieve),
            patch('popgen_genotyping.stages.run_king_robust') as mock_run,
        ):
            KingRobust.queue_jobs(mock_self, mock_cohort, mock_inputs)

        mock_run.assert_called_once_with(
            pgen_path='/agg/COH1.pgen',
            pvar_path='/agg/COH1.pvar',
            psam_path='/agg/COH1.psam',
            remove_list_path='/exclusions/COH1.remove_list.txt',
            kinship_min=0.0442,
            output_kin0_path='/out/COH1_king_robust.kin0.gz',
            output_log_path='/out/COH1_king_robust.log',
            job_name='KingRobust_my_cohort',
        )
        mock_self.make_outputs.assert_called_once_with(mock_cohort, data=expected_outputs, jobs=[mock_run.return_value])
