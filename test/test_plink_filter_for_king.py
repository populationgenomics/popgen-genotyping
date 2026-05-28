"""Tests for the BAFRegress contamination filter sub-job in KingIbdseg.

The job is pure awk + sort + plink: we strip the plink invocation (no plink
binary in the pytest env) and assert against the `remove_samples.tsv` that
the awk pipeline produces. KING 2.3.2 silently ignores `--remove` (see
[[reference-king-232-no-remove-flag]]), so getting this remove-list right is
the only thing standing between a contaminated sample and the IBD graph.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

from popgen_genotyping.jobs.plink_filter_for_king_job import (
    BAFREGRESS_ESTIMATE_COL,
    BAFREGRESS_THRESHOLD,
    run_plink_filter_for_king,
)

# -- Helpers ------------------------------------------------------------------


def _capture_filter_command(
    tmp_path: Path,
    fam_iids: list[tuple[str, str]],
    bafregress_files: dict[str, list[tuple[str, str]]],
) -> tuple[str, Path]:
    """Invoke run_plink_filter_for_king with mocked Hail Batch.

    Writes a fake .fam (FID/IID pairs) and one BAFRegress TSV per entry in
    ``bafregress_files`` to ``tmp_path``, then captures the bash command the
    production code queued.

    Args:
        tmp_path (Path): pytest tmp_path fixture; fixture files are written
            under this directory and the captured bash is expected to be
            executed from here.
        fam_iids (list[tuple[str, str]]): (FID, IID) pairs to write into the
            .fam fixture in the standard 6-column PLINK format.
        bafregress_files (dict[str, list[tuple[str, str]]]): filename ->
            list of (sample_id, baf_regress_string) rows. The estimate is
            written verbatim so callers can drive non-numeric branches.

    Returns:
        tuple[str, Path]: The captured bash and the path at which
            ``remove_samples.tsv`` will be written when the bash is executed
            with ``cwd=tmp_path``.
    """
    fam_path: Path = tmp_path / 'in.fam'
    fam_path.write_text(
        ''.join(f'{fid}\t{iid}\t0\t0\t0\t-9\n' for fid, iid in fam_iids),
    )

    baf_input_paths: list[str] = []
    for filename, rows in bafregress_files.items():
        baf_path: Path = tmp_path / filename
        lines: list[str] = [f'sample_id\t{BAFREGRESS_ESTIMATE_COL}\tNhom\n']
        lines.extend(f'{sid}\t{est}\t1000\n' for sid, est in rows)
        baf_path.write_text(''.join(lines))
        baf_input_paths.append(f'gs://fake-bucket/{filename}')

    mock_batch = MagicMock()
    plink_resource = MagicMock()
    plink_resource.fam = str(fam_path)
    mock_batch.read_input_group.return_value = plink_resource
    # `read_input` redirects each gs:// path onto the local fixture file with
    # the same basename, so the awk loop reads from tmp_path.
    mock_batch.read_input = MagicMock(side_effect=lambda p: str(tmp_path / Path(p).name))

    mock_job = MagicMock()

    with (
        patch('popgen_genotyping.jobs.plink_filter_for_king_job.get_batch', return_value=mock_batch),
        patch('popgen_genotyping.jobs.plink_filter_for_king_job.register_job', return_value=mock_job),
        patch('popgen_genotyping.jobs.plink_filter_for_king_job.config_retrieve', return_value='plink-image:1.0'),
    ):
        run_plink_filter_for_king(
            bed_path='gs://x/in.bed',
            bim_path='gs://x/in.bim',
            fam_path='gs://x/in.fam',
            bafregress_paths=baf_input_paths,
        )

    cmd: str = mock_job.command.call_args[0][0]
    return cmd, tmp_path / 'remove_samples.tsv'


def _strip_plink_invocation(cmd: str, replacement: str = ':') -> str:
    """Replace the trailing ``plink --bfile ...`` invocation with ``replacement``.

    The pytest env has no plink binary and we only need to inspect the
    remove-list construction; the plink ``--make-bed`` is a no-op for these
    tests. The plink call is the last command in the queued bash, so we
    consume from ``plink --bfile`` to end-of-string — that avoids depending
    on the formatted form of the input/output ResourceGroup references
    (Hail Batch's deferred-substitution placeholders contain spaces and
    angle brackets when fixture mocks fall back to their default repr).

    Args:
        cmd (str): The full bash command captured from run_plink_filter_for_king.
        replacement (str): Shell code to substitute for the plink invocation.

    Returns:
        str: The bash command with the plink invocation replaced.
    """
    new, n = re.subn(
        r'plink\s+--bfile\b.*\Z',
        replacement,
        cmd,
        count=1,
        flags=re.DOTALL,
    )
    assert n == 1, 'expected exactly one plink invocation to substitute'
    return new


def _run_filter_bash(tmp_path: Path, cmd: str) -> subprocess.CompletedProcess[bytes]:
    """Execute the captured bash (sans plink) under tmp_path."""
    script: str = _strip_plink_invocation(cmd)
    return subprocess.run(  # noqa: S603
        ['bash', '-c', script],  # noqa: S607
        cwd=tmp_path,
        capture_output=True,
        check=False,
    )


def _read_remove_iids(remove_path: Path) -> set[str]:
    """Return the IIDs (column 2) listed in ``remove_samples.tsv``."""
    if not remove_path.exists():
        return set()
    return {line.split('\t')[1] for line in remove_path.read_text().splitlines() if line}


# -- Tests --------------------------------------------------------------------


class TestFilterRemoveList:
    """Awk + shell remove-list construction is the contract of this job."""

    def test_above_threshold_is_removed(self, tmp_path: Path) -> None:
        """Samples whose baf_regress > BAFREGRESS_THRESHOLD are dropped."""
        cmd, remove_path = _capture_filter_command(
            tmp_path,
            fam_iids=[('F1', 'S1'), ('F2', 'S2')],
            bafregress_files={'baf.tsv': [('S1', '0.05'), ('S2', '0.01')]},
        )
        result = _run_filter_bash(tmp_path, cmd)
        assert result.returncode == 0, result.stderr
        assert _read_remove_iids(remove_path) == {'S1'}

    def test_at_or_below_threshold_is_kept(self, tmp_path: Path) -> None:
        """The 3% boundary is inclusive; clean low values are kept."""
        cmd, remove_path = _capture_filter_command(
            tmp_path,
            fam_iids=[('F', 'S_exact'), ('F', 'S_below'), ('F', 'S_zero')],
            bafregress_files={
                'baf.tsv': [
                    ('S_exact', str(BAFREGRESS_THRESHOLD)),
                    ('S_below', '0.029'),
                    ('S_zero', '0.0'),
                ],
            },
        )
        result = _run_filter_bash(tmp_path, cmd)
        assert result.returncode == 0, result.stderr
        assert _read_remove_iids(remove_path) == set()

    def test_nan_is_removed(self, tmp_path: Path) -> None:
        """BAFRegress failure mode `-nan` is non-numeric and drops the sample."""
        cmd, remove_path = _capture_filter_command(
            tmp_path,
            fam_iids=[('F', 'S_nan'), ('F', 'S_clean')],
            bafregress_files={
                'baf.tsv': [('S_nan', '-nan'), ('S_clean', '0.01')],
            },
        )
        result = _run_filter_bash(tmp_path, cmd)
        assert result.returncode == 0, result.stderr
        assert _read_remove_iids(remove_path) == {'S_nan'}

    def test_non_numeric_garbage_is_removed(self, tmp_path: Path) -> None:
        """Non-matching string values (e.g. `NA`) are treated as failures."""
        cmd, remove_path = _capture_filter_command(
            tmp_path,
            fam_iids=[('F', 'S_na'), ('F', 'S_ok')],
            bafregress_files={'baf.tsv': [('S_na', 'NA'), ('S_ok', '0.005')]},
        )
        result = _run_filter_bash(tmp_path, cmd)
        assert result.returncode == 0, result.stderr
        assert _read_remove_iids(remove_path) == {'S_na'}

    def test_sample_absent_from_bafregress_is_removed(self, tmp_path: Path) -> None:
        """.fam IIDs missing from every BAFRegress output are dropped."""
        cmd, remove_path = _capture_filter_command(
            tmp_path,
            fam_iids=[('F', 'S_in'), ('F', 'S_missing')],
            bafregress_files={'baf.tsv': [('S_in', '0.01')]},
        )
        result = _run_filter_bash(tmp_path, cmd)
        assert result.returncode == 0, result.stderr
        assert _read_remove_iids(remove_path) == {'S_missing'}

    def test_union_across_bafregress_files(self, tmp_path: Path) -> None:
        """A sample passing in any per-cohort file is kept (union semantics)."""
        cmd, remove_path = _capture_filter_command(
            tmp_path,
            fam_iids=[('F', 'S_in_A'), ('F', 'S_in_B'), ('F', 'S_in_neither')],
            bafregress_files={
                'bafA.tsv': [('S_in_A', '0.01')],
                'bafB.tsv': [('S_in_B', '0.02')],
            },
        )
        result = _run_filter_bash(tmp_path, cmd)
        assert result.returncode == 0, result.stderr
        assert _read_remove_iids(remove_path) == {'S_in_neither'}

    def test_fid_iid_emitted_in_plink_remove_format(self, tmp_path: Path) -> None:
        """remove_samples.tsv rows are `FID\\tIID` per `plink --remove`."""
        cmd, remove_path = _capture_filter_command(
            tmp_path,
            fam_iids=[('FAM_A', 'S1')],
            bafregress_files={'baf.tsv': [('S1', '0.99')]},
        )
        result = _run_filter_bash(tmp_path, cmd)
        assert result.returncode == 0, result.stderr
        assert remove_path.read_text() == 'FAM_A\tS1\n'

    def test_missing_required_column_errors(self, tmp_path: Path) -> None:
        """Missing `sample_id` or `baf_regress` header trips the awk ERROR branch."""
        fam_path: Path = tmp_path / 'in.fam'
        fam_path.write_text('F\tS1\t0\t0\t0\t-9\n')
        bad_baf: Path = tmp_path / 'bad.tsv'
        bad_baf.write_text('sample\testimate\nS1\t0.01\n')

        mock_batch = MagicMock()
        plink_resource = MagicMock()
        plink_resource.fam = str(fam_path)
        mock_batch.read_input_group.return_value = plink_resource
        mock_batch.read_input = MagicMock(side_effect=lambda _p: str(bad_baf))
        mock_job = MagicMock()

        with (
            patch('popgen_genotyping.jobs.plink_filter_for_king_job.get_batch', return_value=mock_batch),
            patch('popgen_genotyping.jobs.plink_filter_for_king_job.register_job', return_value=mock_job),
            patch('popgen_genotyping.jobs.plink_filter_for_king_job.config_retrieve', return_value='plink-image:1.0'),
        ):
            run_plink_filter_for_king(
                bed_path='gs://x/in.bed',
                bim_path='gs://x/in.bim',
                fam_path='gs://x/in.fam',
                bafregress_paths=['gs://bad.tsv'],
            )
        captured: str = mock_job.command.call_args[0][0]
        result = _run_filter_bash(tmp_path, captured)
        assert result.returncode != 0
        assert b'ERROR' in result.stderr
