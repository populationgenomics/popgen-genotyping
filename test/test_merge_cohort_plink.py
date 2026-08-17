"""Tests for run_merge_plink: allele-order preservation and the final super-cohort --keep trim."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from popgen_genotyping.jobs.merge_cohort_plink_job import run_merge_plink

_DEFAULT_COHORT_PATHS = [{'bed': 'gs://x/c1.bed', 'bim': 'gs://x/c1.bim', 'fam': 'gs://x/c1.fam'}]


def _run_merge_plink(
    *,
    cohort_plink_paths: list[dict[str, str]] | None = None,
    previous_aggregate_resource: MagicMock | None = None,
    samples_to_remove: list[str] | None = None,
    keep_samples: list[str] | None = None,
) -> tuple[list[str], MagicMock]:
    """Invoke run_merge_plink against a mocked Batch and capture the queued commands.

    Args:
        cohort_plink_paths: New per-plate filesets to merge (defaults to a single cohort).
        previous_aggregate_resource: Optional rolling-aggregate resource group.
        samples_to_remove: Optional withdrawn SGs to drop from the previous aggregate.
        keep_samples: Super-cohort membership to trim to (defaults to a two-SG cohort; pass an
            explicit list — including [] — to exercise the trim/guard behaviour).

    Returns:
        tuple[list[str], MagicMock]: (bash strings passed to j.command() in queue order, to_path mock).
    """
    mock_batch = MagicMock()
    mock_job = MagicMock()
    mock_batch.read_input_group.return_value = MagicMock()

    with (
        patch('popgen_genotyping.jobs.merge_cohort_plink_job.get_batch', return_value=mock_batch),
        patch('popgen_genotyping.jobs.merge_cohort_plink_job.register_job', return_value=mock_job),
        patch('popgen_genotyping.jobs.merge_cohort_plink_job.config_retrieve', return_value='plink-image:1.0'),
        patch('popgen_genotyping.jobs.merge_cohort_plink_job.to_path') as mock_to_path,
    ):
        run_merge_plink(
            cohort_plink_paths=cohort_plink_paths if cohort_plink_paths is not None else _DEFAULT_COHORT_PATHS,
            output_prefix='gs://o/out',
            keep_samples=keep_samples if keep_samples is not None else ['SG1', 'SG2'],
            previous_aggregate_resource=previous_aggregate_resource,
            samples_to_remove=samples_to_remove,
        )

    commands = [call.args[0] for call in mock_job.command.call_args_list]
    return commands, mock_to_path


def test_rolling_aggregate_path_preserves_allele_order_throughout() -> None:
    """Every PLINK 1.9 --make-bed step in the rolling-aggregate path keeps allele order.

    Without --keep-allele-order, PLINK 1.9 resets A1 to the minor allele by
    frequency. The merged cohort then has mixed A1 semantics, and downstream
    reference-panel merges silently flip allele counts at common variants.

    Guards both:
      - the --remove filter step on the previous aggregate, and
      - the final --merge-list step.
    """
    commands, _ = _run_merge_plink(previous_aggregate_resource=MagicMock(), samples_to_remove=['SG_X'])

    # First queued command: the --remove filter on the previous aggregate.
    assert '--remove' in commands[0], 'precondition: first command should be the remove step'
    assert '--keep-allele-order' in commands[0], (
        '--remove --make-bed must keep A1=ALT/A2=REF orientation; without --keep-allele-order '
        'PLINK 1.9 swaps A1 to the minor allele.'
    )

    # Second queued command: the final --merge-list step.
    assert '--merge-list' in commands[1], 'precondition: second command should be the merge-list step'
    assert '--keep-allele-order' in commands[1]


def test_empty_keep_samples_raises() -> None:
    """An empty keep list is a caller/config bug and must fail fast.

    Trimming to the super cohort is mandatory, so an empty membership would otherwise write the
    untrimmed merge as the aggregate: the exact outcome the trim exists to prevent.
    """
    with pytest.raises(ValueError, match='keep_samples is empty'):
        _run_merge_plink(keep_samples=[])


def test_keep_samples_appends_final_keep_trim() -> None:
    """keep_samples queues a final --keep pass to super-cohort membership, allele order preserved."""
    keep = ['SG1', 'SG2']
    commands, _ = _run_merge_plink(keep_samples=keep)

    # The merge runs first, the trim last.
    assert len(commands) == 2
    assert '--keep ' not in commands[0], 'merge step should not filter samples'
    trim = commands[-1]
    assert '--keep ' in trim, 'final step must be the super-cohort --keep trim'
    assert '--keep-allele-order' in trim
    assert '--output-chr chrM' in trim
    # Membership assert: the trimmed .fam must contain exactly len(keep_samples) rows, else
    # plink --keep silently dropped a claimed SG (super ⊆ merged goes unverified otherwise).
    assert 'wc -l' in trim and '-ne 2' in trim, 'final step must assert the kept-sample count'

    # --keep silently ignores IDs absent from the input, so the pass yields the intersection of the
    # merged set and keep_samples, not equality. The in-job count is what turns a shortfall into a
    # failure rather than a quietly smaller aggregate registered against the super cohort.
    assert f'-ne {len(keep)}' in trim, 'retained count must be compared against the full keep-list length'
    assert 'exit 1' in trim, 'a count mismatch must fail the job, not merely log'


def test_keep_samples_writes_fid_iid_list() -> None:
    """The --keep list is written with the FID=0 / IID convention used by the --remove list."""
    _, mock_to_path = _run_merge_plink(keep_samples=['SG1', 'SG2'])

    mock_to_path.assert_called_once_with('gs://o/out_samples_to_keep.txt')
    mock_to_path.return_value.write_text.assert_called_once_with('0\tSG1\n0\tSG2')


def test_duplicate_keep_samples_do_not_inflate_expected_count() -> None:
    """A repeated ID must not raise the expected count: --keep matches on set membership.

    plink emits one .fam row per distinct sample, so counting the raw list would fail the in-job
    check on a duplicate that changed nothing about the output.
    """
    commands, mock_to_path = _run_merge_plink(keep_samples=['SG2', 'SG1', 'SG2'])

    mock_to_path.return_value.write_text.assert_called_once_with('0\tSG1\n0\tSG2')
    assert '-ne 2' in commands[-1], 'expected count must be the number of distinct IDs'


def test_keep_samples_composes_with_rolling_aggregate() -> None:
    """Rolling aggregate + keep_samples: remove, then merge, then the super-cohort trim (in order)."""
    commands, _ = _run_merge_plink(
        previous_aggregate_resource=MagicMock(),
        samples_to_remove=['SG_X'],
        keep_samples=['SG1', 'SG2'],
    )

    assert len(commands) == 3
    assert '--remove' in commands[0]
    assert '--merge-list' in commands[1]
    assert '--keep ' in commands[2]
    assert '--keep-allele-order' in commands[2]
