"""Tests for run_merge_with_reference_panel and MergeWithReferencePanel.queue_jobs."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock, patch

from popgen_genotyping.jobs.merge_with_reference_panel_job import run_merge_with_reference_panel
from popgen_genotyping.stages import MergeCohortPlink, MergeWithReferencePanel

# -- Helpers ------------------------------------------------------------------


def _capture_merge_with_reference_panel_command() -> str:
    """Invoke run_merge_with_reference_panel with mocked Hail Batch; return the queued bash.

    Returns:
        str: The single bash command string the production code queued.
    """
    mock_batch = MagicMock()
    mock_job = MagicMock()
    mock_batch.read_input_group.return_value = MagicMock()

    with (
        patch('popgen_genotyping.jobs.merge_with_reference_panel_job.get_batch', return_value=mock_batch),
        patch('popgen_genotyping.jobs.merge_with_reference_panel_job.register_job', return_value=mock_job),
        patch(
            'popgen_genotyping.jobs.merge_with_reference_panel_job.config_retrieve',
            return_value='plink-image:1.0',
        ),
    ):
        run_merge_with_reference_panel(
            cohort_plink_paths={
                'bed': 'gs://x/cohort.bed',
                'bim': 'gs://x/cohort.bim',
                'fam': 'gs://x/cohort.fam',
            },
            reference_panel_paths={
                'bed': 'gs://ref/ref.bed',
                'bim': 'gs://ref/ref.bim',
                'fam': 'gs://ref/ref.fam',
            },
            fasta_ref_path='gs://ref/genome.fasta',
            expected_contig_style='with_chr',
            expected_variant_id_pattern='^chr[0-9XYM]+:[0-9]+:[ACGT]+:[ACGT]+$',
            output_pgen_prefix='gs://o/merged',
            output_log_path='gs://o/merged.log',
            output_stats_path='gs://o/merged_stats.tsv',
        )

    return mock_job.command.call_args[0][0]


# -- Tests: run_merge_with_reference_panel command shape ----------------------


class TestRunMergeWithReferencePanelCommand:
    """Structural contract of the bash command queued by the job factory."""

    def test_normalize_step_uses_fasta_anchored_ref(self) -> None:
        """NormalizeCohort must set REF authoritatively from the FASTA.

        Without `--ref-from-fa force`, the cohort's BIM A1/A2 orientation
        depends on whatever the upstream produced (subject to the rolling-
        aggregate allele-swap bug). Anchoring against the FASTA is what
        makes this stage robust regardless of #28's merge status.
        """
        cmd: str = _capture_merge_with_reference_panel_command()
        normalize_section = cmd.split('NormalizeCohort')[1].split('ValidateAgainstExpectations')[0]
        assert '--fa' in normalize_section
        assert '--ref-from-fa force' in normalize_section
        assert "--set-all-var-ids '@:#:$r:$a'" in normalize_section
        assert '--output-chr chrM' in normalize_section
        assert '--max-alleles 2' in normalize_section

    def test_validate_step_asserts_both_sides(self) -> None:
        """Validation must check contig style + variant ID pattern on cohort AND reference."""
        cmd: str = _capture_merge_with_reference_panel_command()
        validate_section = cmd.split('ValidateAgainstExpectations')[1].split('Merge (first attempt)')[0]
        # contig-style assertions: one per side
        assert validate_section.count('contig style') == 2
        assert 'normalized_cohort.bim' in validate_section
        # variant-ID-pattern assertions: one per side
        assert validate_section.count('variant IDs') == 2
        # Both assertions exit non-zero on mismatch
        assert validate_section.count('exit 1') >= 4

    def test_all_plink_19_steps_preserve_allele_order(self) -> None:
        """Every PLINK 1.9 step (both bmerges, both retry excludes) keeps allele order.

        Four invocations expected in the queued bash:
          1. first --bmerge attempt
          2. cohort --exclude first_merge.missnp (retry pre-step)
          3. reference --exclude first_merge.missnp (retry pre-step)
          4. final --bmerge on cleaned filesets

        plink2 calls (NormalizeCohort, ConvertToPlink2) do not need this flag
        because plink2's --make-bed preserves REF/ALT orientation by default.
        """
        cmd: str = _capture_merge_with_reference_panel_command()
        assert cmd.count('--bmerge') == 2
        assert cmd.count('--exclude first_merge.missnp') == 2
        assert cmd.count('--keep-allele-order') == 4

    def test_missnp_retry_targets_both_sides_with_exclude(self) -> None:
        """Retry path applies --exclude to both cohort and reference (not just one)."""
        cmd: str = _capture_merge_with_reference_panel_command()
        retry_section = cmd.split('Drop .missnp on both sides, retry')[1].split('ConvertToPlink2')[0]
        assert retry_section.count('--exclude first_merge.missnp') == 2
        # cohort_clean and reference_clean are both produced before the retry merge
        assert '--out cohort_clean' in retry_section
        assert '--out reference_clean' in retry_section
        assert '--bmerge reference_clean' in retry_section

    def test_converts_final_output_to_plink2(self) -> None:
        """Final output is PLINK2 PGEN, written by plink2 --make-pgen."""
        cmd: str = _capture_merge_with_reference_panel_command()
        assert 'plink2 --bfile final_merge' in cmd
        assert '--make-pgen' in cmd


# -- Tests: MergeWithReferencePanel.expected_outputs --------------------------


class TestMergeWithReferencePanelExpectedOutputs:
    """Output-path layout that Metamist analysis_keys=['pgen'] depends on."""

    def test_path_layout_includes_panel_id_and_datestamp(self) -> None:
        """Five outputs under a date-stamped, panel-id-suffixed basename."""
        prefix = Path('/cohort/merge_ref')
        mock_multicohort = MagicMock()
        mock_self = MagicMock()
        mock_self.name = 'MergeWithReferencePanel'

        fixed_now = datetime(2026, 1, 15, tzinfo=timezone.utc)

        with (
            patch('popgen_genotyping.stages.get_output_prefix', return_value=prefix),
            patch('popgen_genotyping.stages.config_retrieve', return_value='hgdp_1kg_v3.1.2'),
            patch('popgen_genotyping.stages.datetime') as mock_datetime,
        ):
            mock_datetime.now.return_value = fixed_now
            result = MergeWithReferencePanel.expected_outputs(mock_self, mock_multicohort)

        base = '20260115_cohort_plus_hgdp_1kg_v3.1.2'
        assert result == {
            'pgen': prefix / f'{base}.pgen',
            'pvar': prefix / f'{base}.pvar',
            'psam': prefix / f'{base}.psam',
            'log': prefix / f'{base}.log',
            'stats': prefix / f'{base}_variant_intersection_stats.tsv',
        }


# -- Tests: MergeWithReferencePanel.queue_jobs --------------------------------


class TestMergeWithReferencePanelQueueJobs:
    """Wire-up between MergeCohortPlink inputs, config, and the job factory."""

    def test_passes_inputs_and_config_to_job(self) -> None:
        """Cohort bed/bim/fam, reference paths, FASTA, and assertions flow through."""
        mock_multicohort = MagicMock()
        mock_multicohort.name = 'my_multicohort'

        cohort_plink: dict[str, Path] = {
            'bed': Path('/merged/cohort.bed'),
            'bim': Path('/merged/cohort.bim'),
            'fam': Path('/merged/cohort.fam'),
        }
        mock_inputs = MagicMock()
        mock_inputs.as_dict.return_value = cohort_plink

        expected_outputs: dict[str, Path] = {
            'pgen': Path('/out/x.pgen'),
            'pvar': Path('/out/x.pvar'),
            'psam': Path('/out/x.psam'),
            'log': Path('/out/x.log'),
            'stats': Path('/out/x_stats.tsv'),
        }
        mock_self = MagicMock()
        mock_self.expected_outputs.return_value = expected_outputs

        ref_cfg = {
            'bed_path': 'gs://r/ref.bed',
            'bim_path': 'gs://r/ref.bim',
            'fam_path': 'gs://r/ref.fam',
            'panel_id': 'hgdp_1kg_v3.1.2',
            'expected_contig_style': 'with_chr',
            'expected_variant_id_pattern': '^chr[0-9XYM]+:[0-9]+:[ACGT]+:[ACGT]+$',
        }

        def fake_config_retrieve(path: list[str], default: object = None) -> object:  # noqa: ARG001
            """Return the reference-panel sub-dict or the FASTA path by config key."""
            if path == ['popgen_genotyping', 'references', 'reference_panel']:
                return ref_cfg
            if path == ['popgen_genotyping', 'references', 'fasta_ref_path']:
                return 'gs://r/genome.fasta'
            return None

        with (
            patch('popgen_genotyping.stages.run_merge_with_reference_panel') as mock_run,
            patch('popgen_genotyping.stages.config_retrieve', side_effect=fake_config_retrieve),
        ):
            MergeWithReferencePanel.queue_jobs(mock_self, mock_multicohort, mock_inputs)

        # Upstream wired from MergeCohortPlink
        mock_inputs.as_dict.assert_called_once_with(target=mock_multicohort, stage=MergeCohortPlink)

        # Job factory called with cohort PLINK, reference paths, FASTA, expectations
        mock_run.assert_called_once_with(
            cohort_plink_paths={
                'bed': '/merged/cohort.bed',
                'bim': '/merged/cohort.bim',
                'fam': '/merged/cohort.fam',
            },
            reference_panel_paths={
                'bed': 'gs://r/ref.bed',
                'bim': 'gs://r/ref.bim',
                'fam': 'gs://r/ref.fam',
            },
            fasta_ref_path='gs://r/genome.fasta',
            expected_contig_style='with_chr',
            expected_variant_id_pattern='^chr[0-9XYM]+:[0-9]+:[ACGT]+:[ACGT]+$',
            output_pgen_prefix='/out/x',
            output_log_path='/out/x.log',
            output_stats_path='/out/x_stats.tsv',
            job_name='MergeWithReferencePanel_my_multicohort',
        )

        # Outputs flow to make_outputs verbatim
        mock_self.make_outputs.assert_called_once_with(
            mock_multicohort,
            data=expected_outputs,
            jobs=[mock_run.return_value],
        )
