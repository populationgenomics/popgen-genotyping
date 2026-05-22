"""Tests for KingIbdseg: placeholder bash backfill + stage expected_outputs."""

from __future__ import annotations

import gzip
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock, patch

from popgen_genotyping.jobs.king_ibdseg_job import (
    _AUTOSOME_SEG_HEADER,
    _X_SEG_HEADER,
    run_king_ibdseg,
)
from popgen_genotyping.stages import KingIbdseg

# -- Helpers ------------------------------------------------------------------


def _capture_run_king_ibdseg_command(tmp_path: Path) -> tuple[str, dict[str, Path]]:
    """Invoke run_king_ibdseg with mocked Hail Batch; return (queued bash, tmp output paths).

    The five KING resource-group outputs are wired to paths under tmp_path so
    that executing the captured command writes to a real local directory.

    Args:
        tmp_path (Path): pytest tmp_path fixture; root for the simulated KING outputs.

    Returns:
        tuple[str, dict[str, Path]]: The bash string the production code queued,
            plus the resource-group key -> local file path mapping.
    """
    outputs: dict[str, Path] = {
        'seg': tmp_path / 'k.seg',
        'segments_gz': tmp_path / 'k.segments.gz',
        'seg_x': tmp_path / 'kX.seg',
        'segments_gz_x': tmp_path / 'kX.segments.gz',
        'log': tmp_path / 'k.log',
    }

    mock_batch = MagicMock()
    plink_input = MagicMock()
    plink_input.bed = str(tmp_path / 'in.bed')
    plink_input.bim = str(tmp_path / 'in.bim')
    plink_input.fam = str(tmp_path / 'in.fam')
    mock_batch.read_input_group.return_value = plink_input

    mock_job = MagicMock()
    for key, path in outputs.items():
        setattr(mock_job.king_outputs, key, str(path))

    with (
        patch('popgen_genotyping.jobs.king_ibdseg_job.get_batch', return_value=mock_batch),
        patch('popgen_genotyping.jobs.king_ibdseg_job.register_job', return_value=mock_job),
        patch('popgen_genotyping.jobs.king_ibdseg_job.config_retrieve', return_value='king-image:1.0'),
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
        )

    cmd: str = mock_job.command.call_args[0][0]
    return cmd, outputs


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
        cmd, outputs = _capture_run_king_ibdseg_command(tmp_path)

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

    def test_preserves_real_king_output(self, tmp_path: Path) -> None:
        """When KING writes real output the `[ -s ]` guards leave it untouched."""
        cmd, outputs = _capture_run_king_ibdseg_command(tmp_path)

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


# -- Tests: KingIbdseg.expected_outputs ---------------------------------------


class TestKingIbdsegExpectedOutputs:
    """Path layout that Metamist `analysis_keys=['seg', 'seg_x']` depends on."""

    def test_path_layout_and_naming(self) -> None:
        """Five outputs under the standard prefix, date-stamped, no dot before X."""
        prefix = Path('/cohort/king_ibdseg')
        mock_multicohort = MagicMock()
        mock_self = MagicMock()
        mock_self.name = 'KingIbdseg'
        fixed_now = datetime(2026, 1, 15, tzinfo=timezone.utc)

        with (
            patch('popgen_genotyping.stages.get_output_prefix', return_value=prefix) as mock_prefix,
            patch('popgen_genotyping.stages.datetime') as mock_datetime,
        ):
            mock_datetime.now.return_value = fixed_now
            result = KingIbdseg.expected_outputs(mock_self, mock_multicohort)

        assert result == {
            'seg': prefix / '20260115_king.seg',
            'segments': prefix / '20260115_king.segments.gz',
            'seg_x': prefix / '20260115_kingX.seg',
            'segments_x': prefix / '20260115_kingX.segments.gz',
            'log': prefix / '20260115_king.log',
        }
        mock_prefix.assert_called_once_with(
            dataset=mock_multicohort.analysis_dataset,
            stage_name='KingIbdseg',
        )
