"""Tests for scripts/aggregate_cluster_stats.py."""

from __future__ import annotations

import io
from typing import TYPE_CHECKING

import pytest

from popgen_genotyping.scripts.aggregate_cluster_stats import (
    HEADER_COLS,
    build_subset_indices,
    load_sex_map,
    moments,
    run,
)

if TYPE_CHECKING:
    from pathlib import Path


def _row_to_dict(row: str) -> dict[str, str]:
    """Parse a single TSV row into a column-keyed dict."""
    return dict(zip(HEADER_COLS, row.split('\t'), strict=True))


def _build_indices(
    sample_order: list[str],
    sex_map: dict[str, str],
) -> dict[str, list[int]]:
    """Shorthand for tests that don't care about file IO."""
    return build_subset_indices(sample_order, sex_map)


# -- moments ------------------------------------------------------------------


class TestMoments:
    """Tests for the (mean, variance) reducer."""

    def test_empty_cell(self) -> None:
        """Zero observations -> both NA."""
        cell = {'n_x': 0.0, 'sum_x': 0.0, 'sumsq_x': 0.0}
        assert moments(cell, 'n_x', 'sum_x', 'sumsq_x') == ('NA', 'NA')

    def test_single_observation(self) -> None:
        """One observation -> mean reported, variance NA."""
        cell = {'n_x': 1.0, 'sum_x': 5.0, 'sumsq_x': 25.0}
        assert moments(cell, 'n_x', 'sum_x', 'sumsq_x') == ('5', 'NA')

    def test_two_observations(self) -> None:
        """Two observations -> sample variance with n-1 denominator."""
        # values 1.0, 3.0 -> mean 2, sumsq 10, var = (10 - 2*4)/1 = 2
        cell = {'n_x': 2.0, 'sum_x': 4.0, 'sumsq_x': 10.0}
        assert moments(cell, 'n_x', 'sum_x', 'sumsq_x') == ('2', '2')

    def test_three_observations(self) -> None:
        """1, 2, 3 -> mean 2, var 1."""
        cell = {'n_x': 3.0, 'sum_x': 6.0, 'sumsq_x': 14.0}
        assert moments(cell, 'n_x', 'sum_x', 'sumsq_x') == ('2', '1')

    def test_zero_variance_when_all_equal(self) -> None:
        """All observations equal -> variance is exactly 0."""
        # Four observations of 0.25: sum=1.0, sumsq=0.25.
        cell = {'n_x': 4.0, 'sum_x': 1.0, 'sumsq_x': 0.25}
        mean, var = moments(cell, 'n_x', 'sum_x', 'sumsq_x')
        assert mean == '0.25'
        assert var == '0'

    def test_clips_tiny_negative_variance(self) -> None:
        """Roundoff making (sumsq - n*mean^2) slightly negative is clipped to 0."""
        # Construct a cell where the analytic variance would be slightly negative.
        # The clipping branch returns '0' regardless of the magnitude.
        cell = {'n_x': 4.0, 'sum_x': 4.0, 'sumsq_x': 1.0}
        _mean, var = moments(cell, 'n_x', 'sum_x', 'sumsq_x')
        assert var == '0'


# -- load_sex_map -------------------------------------------------------------


class TestLoadSexMap:
    """Tests for parsing the sample-id to sex-code TSV."""

    def test_basic(self, tmp_path: Path) -> None:
        """Codes 1 and 2 map to M and F; others are dropped."""
        path = tmp_path / 'sex.tsv'
        path.write_text('S1\t1\nS2\t2\nS3\t0\nS4\t-9\n')
        assert load_sex_map(path) == {'S1': 'M', 'S2': 'F'}

    def test_skips_blank_and_comment_lines(self, tmp_path: Path) -> None:
        """Empty lines and `#`-prefixed lines are ignored."""
        path = tmp_path / 'sex.tsv'
        path.write_text('# header\n\nS1\t1\n\n#S2\t2\nS3\t2\n')
        assert load_sex_map(path) == {'S1': 'M', 'S3': 'F'}

    def test_missing_sex_column(self, tmp_path: Path) -> None:
        """A row with no sex column is treated as unknown."""
        path = tmp_path / 'sex.tsv'
        path.write_text('S1\nS2\t1\n')
        assert load_sex_map(path) == {'S2': 'M'}


# -- build_subset_indices -----------------------------------------------------


class TestBuildSubsetIndices:
    """Tests for the per-subset column-index lists."""

    def test_basic_partition(self) -> None:
        """`all` keeps order; `female`/`male` filter and preserve order."""
        sample_order = ['S1', 'S2', 'S3', 'S4', 'S5']
        sex_map = {'S1': 'M', 'S2': 'F', 'S3': 'F', 'S4': 'M', 'S5': 'F'}
        idx = build_subset_indices(sample_order, sex_map)
        assert idx['all'] == [0, 1, 2, 3, 4]
        assert idx['female'] == [1, 2, 4]
        assert idx['male'] == [0, 3]

    def test_unknown_sex_only_in_all(self) -> None:
        """Samples with no sex entry are present in `all` but not stratified."""
        sample_order = ['S1', 'S2', 'S3']
        sex_map = {'S2': 'F'}
        idx = build_subset_indices(sample_order, sex_map)
        assert idx['all'] == [0, 1, 2]
        assert idx['female'] == [1]
        assert idx['male'] == []


# -- run (the streaming pass) -------------------------------------------------


@pytest.fixture
def five_sample_indices() -> dict[str, list[int]]:
    """5 samples with mixed sexes: M F F M F."""
    sample_order = ['S1', 'S2', 'S3', 'S4', 'S5']
    sex_map = {'S1': 'M', 'S2': 'F', 'S3': 'F', 'S4': 'M', 'S5': 'F'}
    return _build_indices(sample_order, sex_map)


class TestRunAutosomal:
    """Autosomal variants emit a single row with sex_subset='all'."""

    def test_basic_counts_and_moments(self, five_sample_indices: dict[str, list[int]]) -> None:
        """Hand-computed cluster counts and per-cluster (mean, var) for THETA and R."""
        # S1 0/0 0.1 2.0  -> AA
        # S2 0/0 0.2 2.5  -> AA
        # S3 0/1 0.5 2.2  -> AB
        # S4 1/1 0.9 1.8  -> BB
        # S5 ./. .   .    -> nocall
        line = '\t'.join(
            [
                'chr1',
                '100',
                'A',
                'T',
                'rs1',
                '0/0,0.1,2.0',
                '0/0,0.2,2.5',
                '0/1,0.5,2.2',
                '1/1,0.9,1.8',
                './.,.,.',
            ],
        )

        buf = io.StringIO()
        run([line], five_sample_indices, buf)
        lines = buf.getvalue().rstrip('\n').split('\n')

        assert lines[0] == '\t'.join(HEADER_COLS)
        assert len(lines) == 2, 'autosomal variant must emit exactly one data row'

        row = _row_to_dict(lines[1])
        assert row['chrom'] == 'chr1'
        assert row['sex_subset'] == 'all'
        assert row['n_AA'] == '2'
        assert row['n_AB'] == '1'
        assert row['n_BB'] == '1'
        assert row['n_nocall'] == '1'

        # AA: mean=(0.1+0.2)/2=0.15; var=((0.01+0.04)-2*0.15^2)/1=0.005
        assert row['mean_THETA_AA'] == '0.15'
        assert row['var_THETA_AA'] == '0.005'
        # AA R: mean=(2.0+2.5)/2=2.25; var=((4.0+6.25)-2*2.25^2)/1=0.125
        assert row['mean_R_AA'] == '2.25'
        assert row['var_R_AA'] == '0.125'

        # Single-sample clusters: mean reported, variance NA.
        assert row['mean_THETA_AB'] == '0.5'
        assert row['var_THETA_AB'] == 'NA'
        assert row['mean_R_AB'] == '2.2'
        assert row['mean_THETA_BB'] == '0.9'
        assert row['var_THETA_BB'] == 'NA'

    def test_empty_cluster(self, five_sample_indices: dict[str, list[int]]) -> None:
        """A cluster with zero samples reports n=0 and NA moments."""
        # All five samples call AA; AB and BB are empty.
        line = '\t'.join(
            [
                'chr1',
                '200',
                'A',
                'T',
                'rs2',
                '0/0,0.10,2.0',
                '0/0,0.11,2.0',
                '0/0,0.12,2.0',
                '0/0,0.13,2.0',
                '0/0,0.14,2.0',
            ],
        )
        buf = io.StringIO()
        run([line], five_sample_indices, buf)
        row = _row_to_dict(buf.getvalue().rstrip('\n').split('\n')[1])

        assert row['n_AA'] == '5'
        assert row['n_AB'] == '0'
        assert row['n_BB'] == '0'
        assert row['n_nocall'] == '0'
        assert row['mean_THETA_AB'] == 'NA'
        assert row['var_THETA_AB'] == 'NA'
        assert row['mean_R_BB'] == 'NA'

    def test_missing_format_value_keeps_count(self, five_sample_indices: dict[str, list[int]]) -> None:
        """A called sample with `.` for THETA is counted in n but skipped from moments."""
        # S1 has THETA=. (rare in practice but possible); R is fine.
        line = '\t'.join(
            [
                'chr1',
                '300',
                'A',
                'T',
                'rs3',
                '0/0,.,2.0',
                '0/0,0.2,2.5',
                '0/0,0.3,2.6',
                './.,.,.',
                './.,.,.',
            ],
        )
        buf = io.StringIO()
        run([line], five_sample_indices, buf)
        row = _row_to_dict(buf.getvalue().rstrip('\n').split('\n')[1])

        assert row['n_AA'] == '3'  # cluster count includes the THETA-missing sample
        # THETA mean averages only the two parseable samples.
        assert row['mean_THETA_AA'] == '0.25'
        # R was present for all three.
        assert row['mean_R_AA'] == format((2.0 + 2.5 + 2.6) / 3, '.6g')


class TestRunSexStratified:
    """chrX/Y/MT variants emit three rows: all, female, male."""

    def test_chrx_partition(self, five_sample_indices: dict[str, list[int]]) -> None:
        """Female-only and male-only rows reflect the sex partition."""
        # Same five-sample setup as autosomal test, on chrX.
        line = '\t'.join(
            [
                'chrX',
                '500',
                'G',
                'C',
                'rsX1',
                '0/0,0.1,2.0',  # S1 male AA
                '0/0,0.2,2.5',  # S2 female AA
                '0/1,0.5,2.2',  # S3 female AB
                '1/1,0.9,1.8',  # S4 male BB
                './.,.,.',  # S5 female nocall
            ],
        )

        buf = io.StringIO()
        run([line], five_sample_indices, buf)
        data_rows = [_row_to_dict(r) for r in buf.getvalue().rstrip('\n').split('\n')[1:]]

        assert len(data_rows) == 3
        subsets = [r['sex_subset'] for r in data_rows]
        assert subsets == ['all', 'female', 'male']

        all_row, female_row, male_row = data_rows

        # 'all' matches the autosomal expectations.
        assert all_row['n_AA'] == '2'
        assert all_row['n_AB'] == '1'
        assert all_row['n_BB'] == '1'
        assert all_row['n_nocall'] == '1'

        # 'female': S2 (AA), S3 (AB), S5 (nocall)
        assert female_row['n_AA'] == '1'
        assert female_row['n_AB'] == '1'
        assert female_row['n_BB'] == '0'
        assert female_row['n_nocall'] == '1'
        assert female_row['mean_THETA_AA'] == '0.2'
        assert female_row['var_THETA_AA'] == 'NA'  # only one sample
        assert female_row['mean_THETA_AB'] == '0.5'

        # 'male': S1 (AA), S4 (BB)
        assert male_row['n_AA'] == '1'
        assert male_row['n_AB'] == '0'
        assert male_row['n_BB'] == '1'
        assert male_row['n_nocall'] == '0'
        assert male_row['mean_THETA_AA'] == '0.1'
        assert male_row['mean_THETA_BB'] == '0.9'
        assert male_row['mean_THETA_AB'] == 'NA'

    @pytest.mark.parametrize('chrom', ['chrY', 'chrM', 'MT', 'Y'])
    def test_other_sex_strat_chroms_also_split(
        self,
        five_sample_indices: dict[str, list[int]],
        chrom: str,
    ) -> None:
        """chrY/chrM(T) and unprefixed Y all trigger three-row emission."""
        line = '\t'.join(
            [chrom, '10', 'A', 'T', 'v', '0/0,0.1,2.0', '0/0,0.2,2.0', '0/0,0.3,2.0', '1/1,0.9,1.8', './.,.,.'],
        )
        buf = io.StringIO()
        run([line], five_sample_indices, buf)
        subsets = [_row_to_dict(r)['sex_subset'] for r in buf.getvalue().rstrip('\n').split('\n')[1:]]
        assert subsets == ['all', 'female', 'male']


class TestRunHeader:
    """The header is emitted exactly once and matches the declared column tuple."""

    def test_header_only_on_empty_input(self, five_sample_indices: dict[str, list[int]]) -> None:
        """Empty input still emits the header row."""
        buf = io.StringIO()
        run([], five_sample_indices, buf)
        assert buf.getvalue() == '\t'.join(HEADER_COLS) + '\n'
