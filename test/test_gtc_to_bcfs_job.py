"""Tests for GtcToBcfs: the sample_mapping inversion check and run_gtc_to_bcfs's reheadered-BCF assertion."""

from __future__ import annotations

import os
import stat
import subprocess
from typing import TYPE_CHECKING
from unittest.mock import MagicMock, patch

import pytest

from popgen_genotyping.jobs.gtc_to_bcfs_job import run_gtc_to_bcfs
from popgen_genotyping.stages import GtcToBcfs

if TYPE_CHECKING:
    from pathlib import Path

# barcode_position -> SG ID, as built by the GtcToBcfs stage.
_SAMPLE_MAPPING = {
    '210297820108_R01C01': 'SG3',
    '210297820108_R01C02': 'SG1',
    '210297820108_R02C01': 'SG2',
}

# A stand-in for bcftools: `query -l` prints the sample names in $FAKE_BCF_SAMPLES, every other
# subcommand drains stdin and succeeds. Lets the captured command run without real GTC/BCF data.
_FAKE_BCFTOOLS = """#!/usr/bin/env bash
if [ "$1" = "query" ]; then
    cat "$FAKE_BCF_SAMPLES"
else
    cat > /dev/null
fi
"""


def _capture_command(sample_mapping: dict[str, str], tmp_path: Path) -> str:
    """Invoke run_gtc_to_bcfs against a mocked Batch and return the queued bash command.

    The output resource-group paths point under tmp_path so the command can be executed locally.

    Args:
        sample_mapping (dict[str, str]): Mapping of barcode_position to SG ID.
        tmp_path (Path): pytest tmp_path fixture; root for the simulated job outputs.

    Returns:
        str: The bash string passed to j.command().
    """
    mock_batch = MagicMock()
    # Echo the path back so the captured command names real, inspectable files.
    mock_batch.read_input.side_effect = lambda path: path
    mock_batch.read_input_group.return_value.base = 'ref.fasta'
    mock_job = MagicMock()
    mock_job.heavy_bcf.bcf = str(tmp_path / 'heavy.bcf')
    mock_job.light_bcf.bcf = str(tmp_path / 'light.bcf')
    mock_job.metadata_tsv.tsv = str(tmp_path / 'metadata.tsv')

    with (
        patch('popgen_genotyping.jobs.gtc_to_bcfs_job.get_batch', return_value=mock_batch),
        patch('popgen_genotyping.jobs.gtc_to_bcfs_job.register_job', return_value=mock_job),
        patch('popgen_genotyping.jobs.gtc_to_bcfs_job.config_retrieve', return_value='bcftools-image:1.0'),
    ):
        run_gtc_to_bcfs(
            gtc_paths=['gs://x/a.gtc'],
            sample_mapping=sample_mapping,
            output_heavy_bcf_path='gs://o/heavy.bcf',
            output_light_bcf_path='gs://o/light.bcf',
            output_metadata_path='gs://o/meta.tsv',
            bpm_manifest_path='gs://x/manifest.bpm',
            egt_cluster_path='gs://x/clusters.egt',
            fasta_ref_path='gs://x/ref.fasta',
        )

    command: str = mock_job.command.call_args[0][0]
    return command


def _run_command(command: str, bcf_samples: list[str], tmp_path: Path) -> subprocess.CompletedProcess[str]:
    """Execute the captured bash command with a fake bcftools whose BCF holds bcf_samples.

    Args:
        command (str): The bash string captured from run_gtc_to_bcfs.
        bcf_samples (list[str]): Sample names the fake `bcftools query -l` reports for the heavy BCF.
        tmp_path (Path): pytest tmp_path fixture; the command's working directory.

    Returns:
        subprocess.CompletedProcess[str]: The finished process, with stderr captured.
    """
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    fake_bcftools = bin_dir / 'bcftools'
    fake_bcftools.write_text(_FAKE_BCFTOOLS)
    fake_bcftools.chmod(fake_bcftools.stat().st_mode | stat.S_IEXEC)

    samples_file = tmp_path / 'fake_bcf_samples.txt'
    samples_file.write_text(''.join(f'{name}\n' for name in bcf_samples))
    # The command ends by moving the metadata file that `bcftools +gtc2vcf --extra` would write.
    (tmp_path / 'metadata_raw.tsv').write_text('')

    env = {
        **os.environ,
        'PATH': f'{bin_dir}{os.pathsep}{os.environ["PATH"]}',
        'FAKE_BCF_SAMPLES': str(samples_file),
        'BATCH_TMPDIR': str(tmp_path),
    }
    return subprocess.run(  # noqa: S603
        ['bash', '-c', command],  # noqa: S607
        cwd=tmp_path,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )


def test_command_asserts_bcf_samples_equal_expected_sg_ids(tmp_path: Path) -> None:
    """The command writes the sorted expected SG IDs and diffs them against the heavy BCF's samples."""
    command = _capture_command(_SAMPLE_MAPPING, tmp_path)

    # The expected set is the mapping's values (SG IDs), not its keys (barcodes), one per line.
    assert 'cat <<EOF > expected_samples.txt\nSG1\nSG2\nSG3\nEOF\n' in command
    assert f'bcftools query -l {tmp_path / "heavy.bcf"}' in command
    assert 'diff expected_samples.txt actual_samples.txt' in command
    assert 'exit 1' in command

    # The check runs on the heavy BCF immediately after it is written, before the light BCF is derived.
    reheader = command.index('bcftools reheader')
    check = command.index('bcftools query -l')
    annotate = command.index('bcftools annotate')
    assert reheader < check < annotate


def test_command_passes_when_bcf_samples_match(tmp_path: Path) -> None:
    """A BCF holding exactly the expected SG IDs (in any order) lets the command finish."""
    command = _capture_command(_SAMPLE_MAPPING, tmp_path)

    result = _run_command(command, bcf_samples=['SG2', 'SG3', 'SG1'], tmp_path=tmp_path)

    assert result.returncode == 0, result.stderr


def test_command_fails_and_names_unmatched_barcodes(tmp_path: Path) -> None:
    """Names reheader could not map (still raw barcodes) fail the command and appear in the log."""
    command = _capture_command(_SAMPLE_MAPPING, tmp_path)

    # The COH16076 failure mode: no barcode matched, so every sample kept its raw name.
    raw_barcodes = list(_SAMPLE_MAPPING)
    result = _run_command(command, bcf_samples=raw_barcodes, tmp_path=tmp_path)

    assert result.returncode != 0
    assert 'not exactly the expected SG IDs' in result.stderr
    assert 'expected 3, found 3, 6 differing lines' in result.stderr
    for barcode in raw_barcodes:
        assert f'> {barcode}' in result.stderr
    for sg_id in _SAMPLE_MAPPING.values():
        assert f'< {sg_id}' in result.stderr
    # The light BCF is never derived from a mis-named heavy BCF.
    assert 'bcftools annotate' not in result.stderr


def test_command_fails_when_a_sample_is_missing(tmp_path: Path) -> None:
    """A BCF with fewer samples than the cohort fails, even if every sample present is correct."""
    command = _capture_command(_SAMPLE_MAPPING, tmp_path)

    result = _run_command(command, bcf_samples=['SG1', 'SG2'], tmp_path=tmp_path)

    assert result.returncode != 0
    assert 'expected 3, found 2, 1 differing lines' in result.stderr
    assert '< SG3' in result.stderr


def test_command_fails_when_a_sample_is_duplicated_in_the_bcf(tmp_path: Path) -> None:
    """Comparing sorted lists (not sets) means a repeated sample name is caught."""
    command = _capture_command(_SAMPLE_MAPPING, tmp_path)

    result = _run_command(command, bcf_samples=['SG1', 'SG2', 'SG3', 'SG3'], tmp_path=tmp_path)

    assert result.returncode != 0
    assert 'expected 3, found 4, 1 differing lines' in result.stderr


def test_empty_sample_mapping_raises(tmp_path: Path) -> None:
    """With no expected SG IDs the assertion would be vacuous, so fail before queueing anything."""
    with pytest.raises(ValueError, match='sample_mapping is empty'):
        _capture_command({}, tmp_path)


def test_duplicate_sg_ids_in_sample_mapping_raise(tmp_path: Path) -> None:
    """Two barcodes mapping to one SG ID would give the BCF duplicate sample names: fail, naming the SG ID."""
    mapping = {'210297820108_R01C01': 'SG1', '210297820108_R01C02': 'SG1', '210297820108_R02C01': 'SG2'}

    with pytest.raises(ValueError, match=r"same SG ID, duplicated: \['SG1'\]"):
        _capture_command(mapping, tmp_path)


def _queue_gtc_to_bcfs(mapping_data: dict[str, dict[str, str]]) -> MagicMock:
    """Run GtcToBcfs.queue_jobs with the Metamist mapping and the job function mocked.

    Args:
        mapping_data (dict[str, dict[str, str]]): SG ID -> {'gtc': path, 'old_name': barcode_pos}, as
            returned by resolve_cohort_gtc_mapping.

    Returns:
        MagicMock: The mocked run_gtc_to_bcfs, to inspect the sample_mapping it was called with.
    """
    mock_cohort = MagicMock()
    mock_cohort.id = 'COH999'
    mock_self = MagicMock()
    mock_self.expected_outputs.return_value = {
        'heavy_bcf': 'gs://o/heavy.bcf',
        'light_bcf': 'gs://o/light.bcf',
        'metadata_tsv': 'gs://o/meta.tsv',
    }

    with (
        patch('popgen_genotyping.stages.config_retrieve', return_value='gs://x/ref'),
        patch('popgen_genotyping.stages.resolve_cohort_gtc_mapping', return_value=mapping_data),
        patch('popgen_genotyping.stages.run_gtc_to_bcfs') as mock_run,
    ):
        GtcToBcfs.queue_jobs(mock_self, mock_cohort, MagicMock())
    return mock_run


def test_queue_jobs_inverts_mapping_to_barcode_position_keyed_sg_ids() -> None:
    """With unique names, run_gtc_to_bcfs receives barcode_position -> SG ID for every SG."""
    mapping_data = {sg_id: {'gtc': f'gs://x/{sg_id}.gtc', 'old_name': name} for name, sg_id in _SAMPLE_MAPPING.items()}

    mock_run = _queue_gtc_to_bcfs(mapping_data)

    assert mock_run.call_args.kwargs['sample_mapping'] == _SAMPLE_MAPPING
    assert len(mock_run.call_args.kwargs['gtc_paths']) == len(_SAMPLE_MAPPING)


def test_queue_jobs_raises_when_sgs_share_a_barcode_position() -> None:
    """Two SGs with the same barcode_position must fail, not collapse to one dict key.

    A manifest barcode corrupted to scientific notation gives many chips one name, so the plain
    inversion silently kept one SG per name and the expected-sample set stopped matching the cohort.
    """
    mangled = '2.10298E+11_R01C01'
    mapping_data = {
        'SG1': {'gtc': 'gs://x/1.gtc', 'old_name': mangled},
        'SG2': {'gtc': 'gs://x/2.gtc', 'old_name': mangled},
        'SG3': {'gtc': 'gs://x/3.gtc', 'old_name': '210297820108_R02C01'},
    }

    with pytest.raises(
        ValueError, match=r"1 barcode_position name\(s\) shared.*2\.10298E\+11_R01C01: \['SG1', 'SG2'\]"
    ):
        _queue_gtc_to_bcfs(mapping_data)
