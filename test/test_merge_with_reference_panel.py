"""Tests for run_merge_with_reference_panel and MergeWithReferencePanel.queue_jobs."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock, patch

from popgen_genotyping.jobs.merge_with_reference_panel_job import run_merge_with_reference_panel
from popgen_genotyping.stages import MergeWithReferencePanel

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
            cohort_pgen_paths={
                'pgen': 'gs://x/cohort.pgen',
                'pvar': 'gs://x/cohort.pvar',
                'psam': 'gs://x/cohort.psam',
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

    def test_convert_step_round_trips_cohort_pgen_to_plink1(self) -> None:
        """ConvertCohortToPlink1 produces a `cohort_plink1` BED/BIM/FAM via `plink2 --pfile`.

        The normalize, validate, intersect and merge steps run in PLINK 1.9
        and consume `cohort_plink1`.
        """
        cmd: str = _capture_merge_with_reference_panel_command()
        convert_section = cmd.split('ConvertCohortToPlink1')[1].split('NormalizeCohort')[0]
        assert '--pfile' in convert_section
        assert '--make-bed' in convert_section
        assert '--out cohort_plink1' in convert_section
        normalize_section = cmd.split('NormalizeCohort')[1].split('ValidateAgainstExpectations')[0]
        assert '--bfile cohort_plink1' in normalize_section

    def test_normalize_step_uses_fasta_anchored_ref(self) -> None:
        """NormalizeCohort sets REF authoritatively from the FASTA.

        Anchoring REF against the FASTA decouples this stage from any
        upstream A1/A2 orientation in the cohort BIM.
        """
        cmd: str = _capture_merge_with_reference_panel_command()
        normalize_section = cmd.split('NormalizeCohort')[1].split('ValidateAgainstExpectations')[0]
        assert '--fa' in normalize_section
        assert '--ref-from-fa force' in normalize_section
        assert "--set-all-var-ids '@:#:$r:$a'" in normalize_section
        assert '--output-chr chrM' in normalize_section

    def test_normalize_step_restricts_to_biallelic_acgt_snps(self) -> None:
        """NormalizeCohort drops indels, non-ACGT alleles, and multiallelics.

        The reference panel is bi-allelic SNPs only; aligning the cohort
        to the same surface keeps the merge boundary clean.
        """
        cmd: str = _capture_merge_with_reference_panel_command()
        normalize_section = cmd.split('NormalizeCohort')[1].split('ValidateAgainstExpectations')[0]
        assert '--snps-only just-acgt' in normalize_section
        assert '--max-alleles 2' in normalize_section

    def test_normalize_step_drops_all_variants_at_duplicate_positions(self) -> None:
        """Every variant at any position with >1 record is excluded (keeping none).

        `--rm-dup` only matches identical variant IDs, so e.g. (A,C) and
        (A,T) at the same coordinate both survive `--rm-dup` but would
        break `plink --bmerge`. The awk pass keys off chr+pos and emits
        every ID at any duplicated position.
        """
        cmd: str = _capture_merge_with_reference_panel_command()
        normalize_section = cmd.split('NormalizeCohort')[1].split('ValidateAgainstExpectations')[0]
        assert '--out normalized_cohort_pre_dedup' in normalize_section
        assert 'duplicate_position_var_ids.txt' in normalize_section
        assert 'count[$1' in normalize_section
        assert '$4' in normalize_section
        assert 'count[$1' in normalize_section and '> 1' in normalize_section
        assert '--out normalized_cohort' in normalize_section

    def test_normalize_step_drops_strand_ambiguous_sites(self) -> None:
        """{A,T} and {C,G} variants are excluded before the merge.

        FASTA-anchoring cannot resolve strand for these sites, so an
        external panel could disagree on REF/ALT orientation; dropping
        them is the safest pre-merge posture.
        """
        cmd: str = _capture_merge_with_reference_panel_command()
        normalize_section = cmd.split('NormalizeCohort')[1].split('ValidateAgainstExpectations')[0]
        assert 'strand_ambiguous_var_ids.txt' in normalize_section
        assert '$5=="A" && $6=="T"' in normalize_section
        assert '$5=="T" && $6=="A"' in normalize_section
        assert '$5=="C" && $6=="G"' in normalize_section
        assert '$5=="G" && $6=="C"' in normalize_section

    def test_normalize_step_combines_exclude_lists(self) -> None:
        """Duplicate-position and strand-ambiguous IDs are unioned for a single --exclude."""
        cmd: str = _capture_merge_with_reference_panel_command()
        normalize_section = cmd.split('NormalizeCohort')[1].split('ValidateAgainstExpectations')[0]
        assert 'cat duplicate_position_var_ids.txt strand_ambiguous_var_ids.txt' in normalize_section
        assert 'sort -u > normalize_exclude.txt' in normalize_section
        assert '--exclude normalize_exclude.txt' in normalize_section

    def test_validate_step_asserts_both_sides(self) -> None:
        """Validation must check contig style + variant ID pattern on cohort AND reference."""
        cmd: str = _capture_merge_with_reference_panel_command()
        validate_section = cmd.split('ValidateAgainstExpectations')[1].split('ComputeIntersect')[0]
        # contig-style assertions: one per side
        assert validate_section.count('contig style') == 2
        assert 'normalized_cohort.bim' in validate_section
        # variant-ID-pattern assertions: one per side
        assert validate_section.count('variant IDs') == 2
        # Both assertions exit non-zero on mismatch
        assert validate_section.count('exit 1') >= 4

    def test_intersect_step_builds_common_ids_via_awk_id_hash(self) -> None:
        """Common-ID list is computed by awk hash on the variant-ID column ($2) of both BIMs.

        Both BIMs use the FASTA-anchored `chr:pos:ref:alt` ID scheme, so an
        ID match implies an allele match — pre-filtering to the intersect
        keeps `plink --bmerge` conflict-free without a retry path.
        """
        cmd: str = _capture_merge_with_reference_panel_command()
        intersect_section = cmd.split('ComputeIntersect')[1].split('# ---- Merge')[0]
        assert 'NR==FNR' in intersect_section
        assert 'ref[$2]' in intersect_section
        assert '$2 in ref' in intersect_section
        assert 'common.ids' in intersect_section

    def test_intersect_step_fails_fast_on_empty_common_ids(self) -> None:
        """A fully-disjoint panel exits non-zero before bmerge with an explicit message."""
        cmd: str = _capture_merge_with_reference_panel_command()
        intersect_section = cmd.split('ComputeIntersect')[1].split('# ---- Merge')[0]
        assert '[ ! -s common.ids ]' in intersect_section
        assert 'no variants in common' in intersect_section
        assert 'exit 1' in intersect_section

    def test_intersect_step_extracts_both_sides_before_merge(self) -> None:
        """Both cohort and reference are pre-filtered with `--extract common.ids`."""
        cmd: str = _capture_merge_with_reference_panel_command()
        intersect_section = cmd.split('ComputeIntersect')[1].split('ConvertToPlink2')[0]
        assert intersect_section.count('--extract common.ids') == 2
        assert '--out cohort_intersect' in intersect_section
        assert '--out reference_intersect' in intersect_section
        assert '--bmerge reference_intersect' in intersect_section

    def test_all_plink_19_steps_preserve_allele_order(self) -> None:
        """Every PLINK 1.9 step (both extracts, single bmerge) passes `--keep-allele-order`.

        plink2 calls do not need this flag — `plink2 --make-bed` preserves
        REF/ALT orientation by default.
        """
        cmd: str = _capture_merge_with_reference_panel_command()
        assert cmd.count('--bmerge') == 1
        assert cmd.count('--extract common.ids') == 2
        assert cmd.count('--keep-allele-order') == 3

    def test_stats_tsv_emits_intersect_and_side_only_counts(self) -> None:
        """Stats TSV reports the intersect-derived locus accounting."""
        cmd: str = _capture_merge_with_reference_panel_command()
        stats_section = cmd.split('# ---- Stats')[1].split('# ---- Capture log')[0]
        assert 'intersect_variants' in stats_section
        assert 'cohort_only_variants' in stats_section
        assert 'reference_only_variants' in stats_section
        assert 'cohort_variants_dropped_duplicate_position' in stats_section
        assert 'cohort_variants_dropped_strand_ambiguous' in stats_section

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
    """Wire-up between config-driven inputs and the job factory.

    The stage has no `required_stages`; both the cohort PGEN/PVAR/PSAM and the
    reference panel BED/BIM/FAM are pointed at by config.
    """

    def test_passes_inputs_and_config_to_job(self) -> None:
        """Cohort pgen/pvar/psam, reference paths, FASTA, and assertions flow through."""
        mock_multicohort = MagicMock()
        mock_multicohort.name = 'my_multicohort'

        # `inputs` is unused by this stage — it reads everything from config.
        mock_inputs = MagicMock()

        expected_outputs: dict[str, Path] = {
            'pgen': Path('/out/x.pgen'),
            'pvar': Path('/out/x.pvar'),
            'psam': Path('/out/x.psam'),
            'log': Path('/out/x.log'),
            'stats': Path('/out/x_stats.tsv'),
        }
        mock_self = MagicMock()
        mock_self.expected_outputs.return_value = expected_outputs

        cohort_cfg = {
            'pgen_path': 'gs://c/cohort.pgen',
            'pvar_path': 'gs://c/cohort.pvar',
            'psam_path': 'gs://c/cohort.psam',
        }
        ref_cfg = {
            'bed_path': 'gs://r/ref.bed',
            'bim_path': 'gs://r/ref.bim',
            'fam_path': 'gs://r/ref.fam',
            'panel_id': 'hgdp_1kg_v3.1.2',
            'expected_contig_style': 'with_chr',
            'expected_variant_id_pattern': '^chr[0-9XYM]+:[0-9]+:[ACGT]+:[ACGT]+$',
        }

        def fake_config_retrieve(path: list[str], default: object = None) -> object:  # noqa: ARG001
            """Return the cohort_aggregate / reference_panel sub-dict, or the FASTA path."""
            if path == ['popgen_genotyping', 'references', 'cohort_aggregate']:
                return cohort_cfg
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

        # Stage does not consult upstream cpg-flow inputs.
        mock_inputs.as_dict.assert_not_called()

        # Job factory called with cohort PGEN, reference paths, FASTA, expectations
        mock_run.assert_called_once_with(
            cohort_pgen_paths={
                'pgen': 'gs://c/cohort.pgen',
                'pvar': 'gs://c/cohort.pvar',
                'psam': 'gs://c/cohort.psam',
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
