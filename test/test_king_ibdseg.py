"""Tests for KingIbdseg: placeholder bash backfill + stage expected_outputs."""

from __future__ import annotations

import gzip
import re
import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

from popgen_genotyping.jobs.king_ibdseg_job import (
    _AUTOSOME_SEG_HEADER,
    _X_SEG_HEADER,
    run_king_ibdseg,
)
from popgen_genotyping.stages import KingIbdseg, MergeCohortPlink, Plink2Qc

# -- Helpers ------------------------------------------------------------------


def _capture_run_king_ibdseg_commands(tmp_path: Path) -> tuple[str, str, str, dict[str, Path]]:
    """Invoke run_king_ibdseg with mocked Hail Batch; return the three queued commands.

    run_king_ibdseg queues three jobs (exclusions, a plink2 chromosome recode +
    --remove, then KING). The five KING resource-group outputs are wired to
    paths under tmp_path so that executing the captured KING command writes to
    a real local directory.

    Args:
        tmp_path (Path): pytest tmp_path fixture; root for the simulated KING outputs.

    Returns:
        tuple[str, str, str, dict[str, Path]]: the exclusions bash string, the
            recode bash string, the KING bash string, and the KING
            resource-group key -> local file path mapping.
    """
    outputs: dict[str, Path] = {
        'seg': tmp_path / 'k.seg',
        'segments_gz': tmp_path / 'k.segments.gz',
        'seg_x': tmp_path / 'kX.seg',
        'segments_gz_x': tmp_path / 'kX.segments.gz',
        'log': tmp_path / 'k.log',
    }

    mock_batch = MagicMock()
    # read_input is used both to stage the script and to read the bafregress/smiss
    # inputs; echoing the path back lets the captured command assert on file names.
    mock_batch.read_input.side_effect = lambda path: path

    exclusions_job = MagicMock()
    exclusions_job.excluded_tsv = str(tmp_path / 'k.excluded_samples.tsv')
    exclusions_job.remove_list = str(tmp_path / 'remove.txt')
    recode_job = MagicMock()
    king_job = MagicMock()
    for key, path in outputs.items():
        setattr(king_job.king_outputs, key, str(path))

    with (
        patch('popgen_genotyping.jobs.king_ibdseg_job.get_batch', return_value=mock_batch),
        patch('popgen_genotyping.jobs.king_ibdseg_job.get_driver_image', return_value='driver:latest'),
        patch(
            'popgen_genotyping.jobs.king_ibdseg_job.register_job',
            side_effect=[exclusions_job, recode_job, king_job],
        ),
        patch('popgen_genotyping.jobs.king_ibdseg_job.config_retrieve', return_value='img:1.0'),
    ):
        run_king_ibdseg(
            bed_path='gs://x/in.bed',
            bim_path='gs://x/in.bim',
            fam_path='gs://x/in.fam',
            output_seg_path='gs://o/out.seg',
            output_segments_path='gs://o/out.segments.gz',
            output_seg_x_path='gs://o/outX.seg',
            output_segments_x_path='gs://o/outX.segments.gz',
            output_log_path='gs://o/out.log',
            bafregress_paths=['gs://baf/COHP1.BAFRegress.txt'],
            smiss_path='gs://qc/cohort.smiss',
            contamination_max=0.02,
            fmiss_max=0.02,
            output_excluded_path='gs://o/out.excluded_samples.tsv',
        )

    exclusions_cmd: str = exclusions_job.command.call_args[0][0]
    recode_cmd: str = recode_job.command.call_args[0][0]
    king_cmd: str = king_job.command.call_args[0][0]
    return exclusions_cmd, recode_cmd, king_cmd, outputs


def _strip_king_invocation(cmd: str, replacement: str = ':') -> str:
    """Replace the `king ... | tee <log>` line with replacement shell code.

    The pytest environment has no KING binary, and we want full control over
    which of the four KING output files exist on disk when the backfill
    `if [ ! -s ... ]` guards run. The replacement is dropped in verbatim, so
    it can either be a no-op (`:`) or shell code that synthesises KING output.

    Args:
        cmd (str): The full bash command captured from run_king_ibdseg.
        replacement (str): Shell code to substitute for the KING invocation.

    Returns:
        str: The bash command with the KING invocation replaced.
    """
    new, n = re.subn(
        r'king\s+-b\b.*?\|\s*tee\s+\S+',
        replacement,
        cmd,
        count=1,
        flags=re.DOTALL,
    )
    assert n == 1, 'expected exactly one KING invocation to substitute'
    return new


# -- Tests: run_king_ibdseg placeholder backfill ------------------------------


class TestKingIbdsegPlaceholderBackfill:
    """Schema contract for the four header-only placeholder outputs."""

    def test_round_trips_when_king_omits_all_outputs(self, tmp_path: Path) -> None:
        """All four placeholders are produced with the documented schemas."""
        _exclusions_cmd, _recode_cmd, cmd, outputs = _capture_run_king_ibdseg_commands(tmp_path)

        # Replace `king | tee` with a no-op so none of the KING outputs exist.
        script = _strip_king_invocation(cmd, replacement=':')
        subprocess.run(['bash', '-c', script], check=True)  # noqa: S603, S607

        # Autosomal .seg: bash writes the module constant verbatim, and the
        # constant carries the schema downstream parsers depend on.
        seg_text = outputs['seg'].read_text()
        assert seg_text == _AUTOSOME_SEG_HEADER + '\n'
        assert seg_text.rstrip('\n').split('\t') == [
            'FID1',
            'ID1',
            'FID2',
            'ID2',
            'IBD1Seg',
            'IBD2Seg',
            'PropIBD',
            'InfType',
        ]

        # X-chr .seg: same shape, X schema (no InfType; has Sex/MaxIBD).
        seg_x_text = outputs['seg_x'].read_text()
        assert seg_x_text == _X_SEG_HEADER + '\n'
        assert seg_x_text.rstrip('\n').split('\t') == [
            'FID1',
            'ID1',
            'FID2',
            'ID2',
            'Sex1',
            'Sex2',
            'MaxIBD1',
            'MaxIBD2',
            'IBD1Seg',
            'IBD2Seg',
            'PropIBD',
        ]

        # .segments.gz placeholders: valid gzip with empty payload.
        with gzip.open(outputs['segments_gz'], 'rb') as fh:
            assert fh.read() == b''
        with gzip.open(outputs['segments_gz_x'], 'rb') as fh:
            assert fh.read() == b''

    def test_king_failure_fails_the_job_instead_of_backfilling(self, tmp_path: Path) -> None:
        """A non-zero KING exit propagates through `| tee` and no placeholders are written.

        KING exits 1 on a fatal error but the pipeline's status is tee's, so without
        pipefail the job would succeed with header-only relatedness outputs.
        """
        _exclusions_cmd, _recode_cmd, cmd, outputs = _capture_run_king_ibdseg_commands(tmp_path)

        script = _strip_king_invocation(cmd, replacement=f'false | tee {outputs["log"]}')
        result = subprocess.run(['bash', '-c', script], check=False)  # noqa: S603, S607

        assert result.returncode != 0
        assert not outputs['seg'].exists()
        assert not outputs['seg_x'].exists()

    def test_preserves_real_king_output(self, tmp_path: Path) -> None:
        """When KING writes real output the `[ -s ]` guards leave it untouched."""
        _exclusions_cmd, _recode_cmd, cmd, outputs = _capture_run_king_ibdseg_commands(tmp_path)

        seg_path = outputs['seg']
        seg_x_path = outputs['seg_x']
        segments_gz_path = outputs['segments_gz']
        segments_gz_x_path = outputs['segments_gz_x']
        fake_king = (
            f'printf "real-autosome-content\\n" > {seg_path} && '
            f'printf "real-autosome-segments\\n" | gzip -c > {segments_gz_path} && '
            f'printf "real-x-content\\n" > {seg_x_path} && '
            f'printf "real-x-segments\\n" | gzip -c > {segments_gz_x_path}'
        )
        script = _strip_king_invocation(cmd, replacement=fake_king)
        subprocess.run(['bash', '-c', script], check=True)  # noqa: S603, S607

        assert outputs['seg'].read_text() == 'real-autosome-content\n'
        assert outputs['seg_x'].read_text() == 'real-x-content\n'
        with gzip.open(outputs['segments_gz'], 'rb') as fh:
            assert fh.read() == b'real-autosome-segments\n'
        with gzip.open(outputs['segments_gz_x'], 'rb') as fh:
            assert fh.read() == b'real-x-segments\n'

    def test_recodes_chromosomes_to_numeric_before_king(self, tmp_path: Path) -> None:
        """A plink2 recode to numeric chromosome codes precedes KING.

        KING only parses numeric codes (23=X, 24=Y, 26=MT); the merged fileset is
        chr-prefixed, so the recode job must run ``plink2 --output-chr 26`` and
        KING must consume that recoded fileset, not the raw input.
        """
        _exclusions_cmd, recode_cmd, king_cmd, _outputs = _capture_run_king_ibdseg_commands(tmp_path)
        assert 'plink2 --bfile' in recode_cmd
        assert '--output-chr 26' in recode_cmd
        # KING reads the recode job's resource group, never the chr-prefixed input.
        assert 'king -b' in king_cmd


# -- Tests: relatedness-exclusion job wiring -----------------------------------


class TestKingIbdsegExclusionJob:
    """The exclusions job invocation and the recode job's --remove + count assert."""

    def test_exclusions_job_invokes_script_with_all_inputs(self, tmp_path: Path) -> None:
        """The exclusions job runs relatedness_exclusions.py with every threshold and input."""
        exclusions_cmd, _recode_cmd, _king_cmd, _outputs = _capture_run_king_ibdseg_commands(tmp_path)
        assert 'relatedness_exclusions.py' in exclusions_cmd
        assert '--bafregress' in exclusions_cmd
        assert '--smiss' in exclusions_cmd
        assert '--contamination-max 0.02' in exclusions_cmd
        assert '--fmiss-max 0.02' in exclusions_cmd
        assert '--output-tsv' in exclusions_cmd
        assert '--output-remove-list' in exclusions_cmd

    def test_recode_always_passes_remove(self, tmp_path: Path) -> None:
        """The recode job always passes --remove; a header-only list removes nothing."""
        _exclusions_cmd, recode_cmd, _king_cmd, _outputs = _capture_run_king_ibdseg_commands(tmp_path)
        assert '--remove' in recode_cmd

    def test_recode_count_assert_present(self, tmp_path: Path) -> None:
        """The recode command asserts recoded count == input count - (remove-list rows)."""
        _exclusions_cmd, recode_cmd, _king_cmd, _outputs = _capture_run_king_ibdseg_commands(tmp_path)
        assert 'removed_count=$(($(wc -l < ' in recode_cmd
        assert 'expected_count=$((input_count - removed_count))' in recode_cmd
        assert 'if [ "$recoded_count" -ne "$expected_count" ]; then' in recode_cmd
        assert 'exit 1' in recode_cmd


# -- Tests: KingIbdseg.expected_outputs ---------------------------------------


class TestKingIbdsegExpectedOutputs:
    """Path layout that Metamist `analysis_keys=['seg', 'seg_x']` depends on."""

    def test_path_layout_and_naming(self) -> None:
        """Six outputs under the standard prefix, cohort-keyed, no dot before X."""
        prefix = Path('/cohort/king_ibdseg')
        mock_cohort = MagicMock()
        mock_cohort.id = 'COH123'
        mock_self = MagicMock()
        mock_self.name = 'KingIbdseg'

        with patch('popgen_genotyping.stages.get_output_prefix', return_value=prefix) as mock_prefix:
            result = KingIbdseg.expected_outputs(mock_self, mock_cohort)

        assert result == {
            'seg': prefix / 'COH123_king.seg',
            'segments': prefix / 'COH123_king.segments.gz',
            'seg_x': prefix / 'COH123_kingX.seg',
            'segments_x': prefix / 'COH123_kingX.segments.gz',
            'log': prefix / 'COH123_king.log',
            'excluded': prefix / 'COH123_king.excluded_samples.tsv',
        }
        mock_prefix.assert_called_once_with(
            dataset=mock_cohort.dataset,
            stage_name='KingIbdseg',
        )


# -- Tests: KingIbdseg.queue_jobs ---------------------------------------------


class TestKingIbdsegQueueJobs:
    """Wire-up between MergeCohortPlink/Plink2Qc inputs, expected_outputs, and the job factory."""

    def test_passes_merged_plink_and_outputs_to_job(self) -> None:
        """Bed/bim/fam and .smiss map through; BafRegress paths are resolved and deduped;
        thresholds come from config; everything reaches run_king_ibdseg in its named slot."""
        mock_cohort = MagicMock()
        mock_cohort.name = 'my_cohort'
        mock_cohort.get_sequencing_group_ids.return_value = ['CPG1', 'CPG2', 'CPG3']

        merged_plink: dict[str, Path] = {
            'bed': Path('/merged/cohort.bed'),
            'bim': Path('/merged/cohort.bim'),
            'fam': Path('/merged/cohort.fam'),
        }

        def as_dict(target: object, stage: object) -> dict[str, Path]:
            del target, stage
            return merged_plink

        def as_path(target: object, stage: object, key: str) -> Path:
            del target
            assert stage is Plink2Qc
            assert key == 'smiss'
            return Path('/qc/cohort.smiss')

        mock_inputs = MagicMock()
        mock_inputs.as_dict.side_effect = as_dict
        mock_inputs.as_path.side_effect = as_path

        expected_outputs: dict[str, Path] = {
            'seg': Path('/out/20260115_king.seg'),
            'segments': Path('/out/20260115_king.segments.gz'),
            'seg_x': Path('/out/20260115_kingX.seg'),
            'segments_x': Path('/out/20260115_kingX.segments.gz'),
            'log': Path('/out/20260115_king.log'),
            'excluded': Path('/out/20260115_king.excluded_samples.tsv'),
        }
        mock_self = MagicMock()
        mock_self.expected_outputs.return_value = expected_outputs

        bafregress_map = {
            'CPG1': 'gs://baf/COHP1.BAFRegress.txt',
            'CPG2': 'gs://baf/COHP1.BAFRegress.txt',
            'CPG3': 'gs://baf/COHP2.BAFRegress.txt',
        }

        thresholds: dict[str, float] = {'contamination_max': 0.02, 'fmiss_max': 0.03}

        def config_retrieve(path: list[str]) -> float:
            assert path[:2] == ['popgen_genotyping', 'king_ibdseg']
            return thresholds[path[2]]

        with (
            patch('popgen_genotyping.stages.resolve_bafregress_map', return_value=bafregress_map) as mock_baf,
            patch('popgen_genotyping.stages.config_retrieve', side_effect=config_retrieve),
            patch('popgen_genotyping.stages.run_king_ibdseg') as mock_run,
        ):
            KingIbdseg.queue_jobs(mock_self, mock_cohort, mock_inputs)

        # Inputs are pulled from the right upstream stages.
        mock_inputs.as_dict.assert_called_once_with(target=mock_cohort, stage=MergeCohortPlink)

        # BafRegress is resolved for the full cohort membership and deduped before passing on.
        mock_baf.assert_called_once_with(sg_ids=['CPG1', 'CPG2', 'CPG3'])

        # The job factory receives plink inputs, the merged-cohort .smiss, both thresholds,
        # deduped BafRegress paths, and the excluded-samples output path in their named slots
        # -- this is the seam where cross-wiring would manifest.
        mock_run.assert_called_once_with(
            bed_path='/merged/cohort.bed',
            bim_path='/merged/cohort.bim',
            fam_path='/merged/cohort.fam',
            output_seg_path='/out/20260115_king.seg',
            output_segments_path='/out/20260115_king.segments.gz',
            output_seg_x_path='/out/20260115_kingX.seg',
            output_segments_x_path='/out/20260115_kingX.segments.gz',
            output_log_path='/out/20260115_king.log',
            bafregress_paths=['gs://baf/COHP1.BAFRegress.txt', 'gs://baf/COHP2.BAFRegress.txt'],
            smiss_path='/qc/cohort.smiss',
            contamination_max=0.02,
            fmiss_max=0.03,
            output_excluded_path='/out/20260115_king.excluded_samples.tsv',
            job_name='KingIbdseg_my_cohort',
        )

        # The same outputs dict flows through to make_outputs, and the job list
        # returned by run_king_ibdseg (exclusions + recode + KING) is passed through verbatim.
        mock_self.make_outputs.assert_called_once_with(
            mock_cohort,
            data=expected_outputs,
            jobs=mock_run.return_value,
        )
