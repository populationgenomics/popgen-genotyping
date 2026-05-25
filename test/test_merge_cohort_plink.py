"""Tests for run_merge_plink: allele-order preservation on the rolling-aggregate path."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

from popgen_genotyping.jobs.merge_cohort_plink_job import run_merge_plink


def _capture_merge_plink_commands_with_rolling_aggregate() -> list[str]:
    """Invoke run_merge_plink with the rolling-aggregate branch active.

    Returns:
        list[str]: Bash strings passed to j.command(), in queue order.
    """
    mock_batch = MagicMock()
    mock_job = MagicMock()
    prev_resource = MagicMock()
    mock_batch.read_input_group.return_value = MagicMock()

    with (
        patch('popgen_genotyping.jobs.merge_cohort_plink_job.get_batch', return_value=mock_batch),
        patch('popgen_genotyping.jobs.merge_cohort_plink_job.register_job', return_value=mock_job),
        patch('popgen_genotyping.jobs.merge_cohort_plink_job.config_retrieve', return_value='plink-image:1.0'),
        patch('popgen_genotyping.jobs.merge_cohort_plink_job.to_path'),
    ):
        run_merge_plink(
            cohort_plink_paths=[{'bed': 'gs://x/c1.bed', 'bim': 'gs://x/c1.bim', 'fam': 'gs://x/c1.fam'}],
            output_prefix='gs://o/out',
            previous_aggregate_resource=prev_resource,
            samples_to_remove=['SG_X'],
        )

    return [call.args[0] for call in mock_job.command.call_args_list]


def test_rolling_aggregate_path_preserves_allele_order_throughout() -> None:
    """Every PLINK 1.9 --make-bed step in the rolling-aggregate path keeps allele order.

    Without --keep-allele-order, PLINK 1.9 resets A1 to the minor allele by
    frequency. The merged cohort then has mixed A1 semantics, and downstream
    reference-panel merges silently flip allele counts at common variants.

    Guards both:
      - the --remove filter step on the previous aggregate, and
      - the final --merge-list step.
    """
    commands: list[str] = _capture_merge_plink_commands_with_rolling_aggregate()

    # First queued command: the --remove filter on the previous aggregate.
    assert '--remove' in commands[0], 'precondition: first command should be the remove step'
    assert '--keep-allele-order' in commands[0], (
        '--remove --make-bed must keep A1=ALT/A2=REF orientation; without --keep-allele-order '
        'PLINK 1.9 swaps A1 to the minor allele.'
    )

    # Second queued command: the final --merge-list step.
    assert '--merge-list' in commands[1], 'precondition: second command should be the merge-list step'
    assert '--keep-allele-order' in commands[1]
