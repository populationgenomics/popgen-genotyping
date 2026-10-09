"""Unit tests for the IBIS relatedness table functions, on synthetic segments only."""

import gzip
import json
from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from popgen_genotyping.scripts import ibis_relatedness_summary as summary

# A short map keeps the kinship arithmetic readable: 100 cM of IBD1 is 100 / 4 / 400 = 0.0625.
MAP_LENGTH_CM = 400.0
CUTOFFS = [5.0, 10.0, 15.0]
MIN_SEG_CM = 5.0
JOIN_GAP_CM = 0.5
# Gap `segments` leaves between segments of one pair and chromosome: more than JOIN_GAP_CM, so they stay apart.
SEPARATION_CM = 10.0


def segments(rows: list[tuple[str, str, str, str, float]]) -> pd.DataFrame:
    """IBIS segments as (ID1, ID2, CHROM, IBD_TYPE, LENGTH_CM) rows, each its own run.

    Segments of one pair on one chromosome are laid end to end SEPARATION_CM apart, so none joins.
    """
    next_start: dict[tuple[str, str, str], float] = {}
    placed = []
    for id1, id2, chromosome, ibd_type, length in rows:
        key = (min(id1, id2), max(id1, id2), chromosome)
        start = next_start.get(key, 0.0)
        placed.append((id1, id2, chromosome, ibd_type, start, start + length, length))
        next_start[key] = start + length + SEPARATION_CM
    return pd.DataFrame(placed, columns=list(summary.IBIS_SEG_USED_COLUMNS))


def pieces(rows: list[tuple[str, str, str, str, float, float]]) -> pd.DataFrame:
    """IBIS segments as (ID1, ID2, CHROM, IBD_TYPE, CM_START, CM_END) rows, for placing them against each other."""
    return pd.DataFrame(
        [
            (id1, id2, chromosome, ibd_type, start, end, end - start)
            for id1, id2, chromosome, ibd_type, start, end in rows
        ],
        columns=list(summary.IBIS_SEG_USED_COLUMNS),
    )


def king(rows: list[tuple[str, str, float, float, str]]) -> pd.DataFrame:
    """KING .seg rows as (ID1, ID2, IBD1Seg, IBD2Seg, InfType)."""
    return pd.DataFrame(rows, columns=['ID1', 'ID2', 'IBD1Seg', 'IBD2Seg', 'InfType'])


NO_KING = king([])
NO_X = segments([])


def build(
    autosomal: pd.DataFrame,
    x: pd.DataFrame = NO_X,
    king_seg: pd.DataFrame = NO_KING,
    cutoffs: list[float] = CUTOFFS,
    min_seg_cm: float = MIN_SEG_CM,
    join_gap_cm: float = JOIN_GAP_CM,
) -> pd.DataFrame:
    return summary.build_relatedness_table(autosomal, x, king_seg, MAP_LENGTH_CM, cutoffs, min_seg_cm, join_gap_cm)


def row(table: pd.DataFrame, id1: str, id2: str) -> pd.Series:
    selected = table[(table['ID1'] == id1) & (table['ID2'] == id2)]
    assert len(selected) == 1
    return selected.iloc[0]


# -- kinship and degree per cutoff -----------------------------------------------------------


def test_one_long_segment_versus_many_short_ones():
    """One 100 cM segment is related at every cutoff; twenty 6 cM segments drop out at >=10 cM."""
    long_pair = [('A', 'B', '1', 'IBD1', 100.0)]
    short_pair = [('C', 'D', str(1 + i % 22), 'IBD1', 6.0) for i in range(20)]
    table = build(segments(long_pair + short_pair))

    long_row = row(table, 'A', 'B')
    for label in ('5', '10', '15'):
        assert long_row[f'IBIS_KINSHIP_GE{label}CM'] == pytest.approx(0.0625)
        assert long_row[f'IBIS_DEGREE_GE{label}CM'] == '3rd'

    short_row = row(table, 'C', 'D')
    assert short_row['IBIS_KINSHIP_GE5CM'] == pytest.approx(120 / 4 / MAP_LENGTH_CM)
    assert short_row['IBIS_DEGREE_GE5CM'] == '3rd'
    assert short_row['IBIS_KINSHIP_GE10CM'] == 0.0
    assert pd.isna(short_row['IBIS_DEGREE_GE10CM'])
    assert short_row['IBIS_KINSHIP_GE15CM'] == 0.0
    assert pd.isna(short_row['IBIS_DEGREE_GE15CM'])


def test_segment_counts_totals_and_longest():
    table = build(
        segments(
            [
                ('A', 'B', '1', 'IBD1', 50.0),
                ('A', 'B', '2', 'IBD1', 9.0),
                ('A', 'B', '3', 'IBD2', 7.0),
                ('A', 'B', '4', 'IBD1', 5.5),
            ],
        ),
    )
    pair = row(table, 'A', 'B')
    assert pair['IBIS_N_SEG'] == 4
    assert pair['IBIS_TOTAL_CM'] == pytest.approx(71.5)
    assert pair['IBIS_LONGEST_CM'] == pytest.approx(50.0)
    assert pair['IBIS_N_SEG_GE7CM'] == 3
    assert pair['IBIS_N_SEG_GE10CM'] == 1


def test_ibd2_counts_double_ibd1_in_kinship():
    table = build(
        segments(
            [
                ('A', 'B', '1', 'IBD1', 80.0),
                ('C', 'D', '1', 'IBD2', 80.0),
                ('E', 'F', '1', 'IBD1', 40.0),
                ('E', 'F', '2', 'IBD2', 40.0),
            ],
        ),
    )
    assert row(table, 'A', 'B')['IBIS_KINSHIP_GE5CM'] == pytest.approx(80 / 4 / MAP_LENGTH_CM)
    assert row(table, 'C', 'D')['IBIS_KINSHIP_GE5CM'] == pytest.approx(80 / 2 / MAP_LENGTH_CM)
    assert row(table, 'E', 'F')['IBIS_KINSHIP_GE5CM'] == pytest.approx((40 / 4 + 40 / 2) / MAP_LENGTH_CM)


def test_degrees_follow_king_thresholds():
    """IBD2 over 360 cM of a 400 cM map is kinship 0.45 (Dup/MZ); 100 cM of IBD2 is 0.125 (2nd)."""
    table = build(
        segments(
            [
                ('A', 'B', '1', 'IBD2', 360.0),
                ('C', 'D', '1', 'IBD1', 140.0),
                ('E', 'F', '1', 'IBD2', 100.0),
                ('G', 'H', '1', 'IBD1', 40.0),
            ],
        ),
    )
    assert row(table, 'A', 'B')['IBIS_DEGREE_GE5CM'] == 'Dup/MZ'
    assert row(table, 'C', 'D')['IBIS_KINSHIP_GE5CM'] == pytest.approx(0.0875)
    assert row(table, 'C', 'D')['IBIS_DEGREE_GE5CM'] == '3rd'
    assert row(table, 'E', 'F')['IBIS_DEGREE_GE5CM'] == '2nd'
    assert pd.isna(row(table, 'G', 'H')['IBIS_DEGREE_GE5CM'])  # 40 / 4 / 400 = 0.025, unrelated


def test_kinship_to_degree_boundaries():
    kinship = pd.Series([0.4, 0.354, 0.2, 0.177, 0.1, 0.0884, 0.05, 0.0442, 0.04, 0.0, float('nan')])
    expected = ['Dup/MZ', '1st', '1st', '1st', '2nd', '2nd', '3rd', '3rd', None, None, None]
    assert summary.kinship_to_degree(kinship).tolist() == [value if value else pd.NA for value in expected]


def test_king_kinship_is_ibd1_over_four_plus_ibd2_over_two():
    kinship = summary.king_kinship(pd.Series([0.5, 1.0, 0.0]), pd.Series([0.25, 0.0, 1.0]))
    assert kinship.tolist() == pytest.approx([0.25, 0.25, 0.5])


# -- KING and IBIS pairs combined ----------------------------------------------------------


def test_king_only_pair_has_na_ibis_columns():
    table = build(
        segments([('A', 'B', '1', 'IBD1', 100.0)]),
        king_seg=king([('C', 'D', 0.5, 0.25, 'FS')]),
    )
    king_only = row(table, 'C', 'D')
    assert king_only['KING_INFTYPE'] == 'FS'
    assert king_only['KING_KINSHIP'] == pytest.approx(0.25)
    ibis_columns = [column for column in table.columns if column.startswith(('IBIS_', 'X_'))]
    assert king_only[ibis_columns].isna().all()


def test_ibis_only_pair_has_na_king_columns():
    table = build(segments([('A', 'B', '1', 'IBD1', 100.0)]))
    ibis_only = row(table, 'A', 'B')
    assert pd.isna(ibis_only['KING_INFTYPE'])
    assert pd.isna(ibis_only['KING_KINSHIP'])


def test_pair_in_both_is_one_row():
    table = build(
        segments([('A', 'B', '1', 'IBD1', 100.0)]),
        king_seg=king([('A', 'B', 0.3, 0.0, '3rd')]),
    )
    assert len(table) == 1
    pair = table.iloc[0]
    assert pair['KING_INFTYPE'] == '3rd'
    assert pair['IBIS_N_SEG'] == 1


def test_pair_order_is_normalised_across_tools():
    """IBIS lists pairs by sample index and KING in its own order; both must land on ID1 < ID2."""
    table = build(
        segments([('B', 'A', '1', 'IBD1', 100.0), ('A', 'B', '2', 'IBD1', 50.0)]),
        king_seg=king([('B', 'A', 0.3, 0.0, '3rd')]),
    )
    assert len(table) == 1
    pair = table.iloc[0]
    assert (pair['ID1'], pair['ID2']) == ('A', 'B')
    assert pair['IBIS_N_SEG'] == 2
    assert pair['KING_INFTYPE'] == '3rd'


def test_duplicate_king_pair_is_rejected():
    with pytest.raises(ValueError, match='more than once'):
        build(NO_X.iloc[0:0], king_seg=king([('A', 'B', 0.3, 0.0, '3rd'), ('B', 'A', 0.3, 0.0, '3rd')]))


def test_unknown_king_inftype_is_rejected(tmp_path: Path):
    path = tmp_path / 'king.seg'
    path.write_text('FID1\tID1\tFID2\tID2\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType\n0\tA\t0\tB\t0.1\t0\t0.1\tCousin\n')
    with pytest.raises(ValueError, match='InfType'):
        summary.read_king_seg(path)


def test_table_is_sorted_and_has_the_documented_columns():
    table = build(
        segments([('C', 'D', '1', 'IBD1', 100.0), ('A', 'B', '1', 'IBD1', 100.0)]),
        king_seg=king([('A', 'E', 0.3, 0.0, '3rd')]),
    )
    assert list(table.columns) == summary.ibis_columns(CUTOFFS)
    assert table[['ID1', 'ID2']].to_numpy().tolist() == [['A', 'B'], ['A', 'E'], ['C', 'D']]
    assert 'IBIS_KINSHIP_GE5CM' in table.columns
    assert 'IBIS_DEGREE_GE15CM' in table.columns


# -- X chromosome ---------------------------------------------------------------------------


def test_x_segments_stay_out_of_autosomal_kinship():
    autosomal = segments([('A', 'B', '1', 'IBD1', 100.0)])
    x = segments([('A', 'B', 'X', 'IBD1', 80.0), ('A', 'B', 'X', 'IBD2', 30.0)])
    with_x = row(build(autosomal, x), 'A', 'B')
    without_x = row(build(autosomal), 'A', 'B')

    assert with_x['IBIS_KINSHIP_GE5CM'] == without_x['IBIS_KINSHIP_GE5CM']
    assert with_x['IBIS_N_SEG'] == 1
    assert with_x['IBIS_TOTAL_CM'] == pytest.approx(100.0)
    assert with_x['X_N_SEG'] == 2
    assert with_x['X_TOTAL_CM'] == pytest.approx(110.0)
    assert with_x['X_LONGEST_CM'] == pytest.approx(80.0)
    assert pd.isna(without_x['X_N_SEG'])


def test_pair_with_only_x_segments_is_not_listed():
    table = build(segments([('A', 'B', '1', 'IBD1', 100.0)]), segments([('C', 'D', 'X', 'IBD1', 40.0)]))
    assert table[['ID1', 'ID2']].to_numpy().tolist() == [['A', 'B']]


# -- segment length floor and cutoffs ----------------------------------------------------------


def test_segments_below_min_seg_cm_are_ignored():
    """IBIS keeps IBD2 segments down to 2 cM and IBD1 pieces beside them; none count below the floor."""
    table = build(
        segments(
            [
                ('A', 'B', '1', 'IBD1', 60.0),
                ('A', 'B', '2', 'IBD2', 3.0),
                ('A', 'B', '3', 'IBD1', 1.5),
                ('C', 'D', '1', 'IBD2', 4.0),
            ],
        ),
    )
    assert table[['ID1', 'ID2']].to_numpy().tolist() == [['A', 'B']]
    assert row(table, 'A', 'B')['IBIS_N_SEG'] == 1
    assert row(table, 'A', 'B')['IBIS_TOTAL_CM'] == pytest.approx(60.0)


# -- joining abutting segments into runs ------------------------------------------------------


def test_abutting_ibd1_ibd2_ibd1_pieces_join_into_one_run():
    """Pieces of 6, 5 and 7 cM fail the 10 and 15 cM cutoffs one by one; as an 18 cM run they pass both."""
    table = build(
        pieces(
            [
                ('A', 'B', '1', 'IBD1', 11.0, 18.0),  # out of order on purpose
                ('A', 'B', '1', 'IBD1', 0.0, 6.0),
                ('A', 'B', '1', 'IBD2', 6.0, 11.0),
            ],
        ),
    )
    pair = row(table, 'A', 'B')
    assert pair['IBIS_N_SEG'] == 1
    assert pair['IBIS_TOTAL_CM'] == pytest.approx(18.0)
    assert pair['IBIS_LONGEST_CM'] == pytest.approx(18.0)
    assert pair['IBIS_N_SEG_GE7CM'] == 1
    assert pair['IBIS_N_SEG_GE10CM'] == 1
    expected_kinship = ((6 + 7) / 4 + 5 / 2) / MAP_LENGTH_CM  # the weights are those of the pieces, not the run
    for label in ('5', '10', '15'):
        assert pair[f'IBIS_KINSHIP_GE{label}CM'] == pytest.approx(expected_kinship)


def test_a_run_passes_a_cutoff_as_a_whole_or_not_at_all():
    """A 12 cM run is below the 15 cM cutoff with all of its kinship, and above 10 cM with all of it."""
    table = build(
        pieces(
            [
                ('A', 'B', '1', 'IBD1', 0.0, 8.0),
                ('A', 'B', '1', 'IBD2', 8.0, 12.0),
            ],
        ),
    )
    pair = row(table, 'A', 'B')
    assert pair['IBIS_KINSHIP_GE10CM'] == pytest.approx((8 / 4 + 4 / 2) / MAP_LENGTH_CM)
    assert pair['IBIS_KINSHIP_GE15CM'] == 0.0
    assert pair['IBIS_N_SEG_GE10CM'] == 1


def test_segments_join_only_when_the_gap_is_below_join_gap_cm():
    table = build(
        pieces(
            [
                ('A', 'B', '1', 'IBD1', 0.0, 20.0),
                ('A', 'B', '1', 'IBD1', 20.5, 40.0),  # gap 0.5 = join_gap_cm: not joined
                ('C', 'D', '1', 'IBD1', 0.0, 20.0),
                ('C', 'D', '1', 'IBD1', 20.25, 40.0),  # gap 0.25: joined
            ],
        ),
    )
    apart = row(table, 'A', 'B')
    assert apart['IBIS_N_SEG'] == 2
    assert apart['IBIS_TOTAL_CM'] == pytest.approx(39.5)
    assert apart['IBIS_LONGEST_CM'] == pytest.approx(20.0)

    joined = row(table, 'C', 'D')
    assert joined['IBIS_N_SEG'] == 1
    assert joined['IBIS_TOTAL_CM'] == pytest.approx(40.0)  # the span, so it includes the 0.25 cM gap
    assert joined['IBIS_KINSHIP_GE5CM'] == pytest.approx((20 + 19.75) / 4 / MAP_LENGTH_CM)  # but the gap has no kinship


def test_join_gap_cm_of_zero_joins_nothing():
    table = build(
        pieces([('A', 'B', '1', 'IBD1', 0.0, 20.0), ('A', 'B', '1', 'IBD2', 20.0, 30.0)]),
        join_gap_cm=0.0,
    )
    assert row(table, 'A', 'B')['IBIS_N_SEG'] == 2


def test_negative_join_gap_cm_is_rejected():
    with pytest.raises(ValueError, match='join_gap_cm'):
        build(segments([('A', 'B', '1', 'IBD1', 60.0)]), join_gap_cm=-0.1)


def test_runs_are_per_pair_and_per_chromosome():
    """Segments that abut in cM but sit on different chromosomes, or belong to different pairs, never join."""
    table = build(
        pieces(
            [
                ('A', 'B', '1', 'IBD1', 0.0, 10.0),
                ('A', 'B', '2', 'IBD1', 10.0, 20.0),
                ('B', 'A', '3', 'IBD1', 0.0, 10.0),  # same pair as A-B, listed the other way round
                ('A', 'C', '1', 'IBD1', 10.0, 20.0),
            ],
        ),
    )
    assert row(table, 'A', 'B')['IBIS_N_SEG'] == 3
    assert row(table, 'A', 'C')['IBIS_N_SEG'] == 1


def test_pair_listed_both_ways_joins_into_one_run():
    table = build(
        pieces([('A', 'B', '1', 'IBD1', 0.0, 10.0), ('B', 'A', '1', 'IBD2', 10.0, 20.0)]),
    )
    pair = row(table, 'A', 'B')
    assert pair['IBIS_N_SEG'] == 1
    assert pair['IBIS_TOTAL_CM'] == pytest.approx(20.0)


def test_short_pieces_inside_a_long_run_are_not_dropped_by_the_floor():
    """A 2 cM IBD2 piece and a 1 cM IBD1 piece are below min_seg_cm, but part of an 83 cM run they still count."""
    table = build(
        pieces(
            [
                ('A', 'B', '1', 'IBD1', 0.0, 40.0),
                ('A', 'B', '1', 'IBD2', 40.0, 42.0),
                ('A', 'B', '1', 'IBD1', 42.0, 82.0),
                ('A', 'B', '1', 'IBD1', 82.0, 83.0),
            ],
        ),
    )
    pair = row(table, 'A', 'B')
    assert pair['IBIS_N_SEG'] == 1
    assert pair['IBIS_TOTAL_CM'] == pytest.approx(83.0)
    assert pair['IBIS_KINSHIP_GE15CM'] == pytest.approx(((40 + 40 + 1) / 4 + 2 / 2) / MAP_LENGTH_CM)


def test_run_shorter_than_min_seg_cm_is_dropped_but_a_longer_one_of_short_pieces_is_kept():
    table = build(
        pieces(
            [
                ('A', 'B', '1', 'IBD1', 0.0, 3.0),  # a 4 cM run: below the 5 cM floor
                ('A', 'B', '1', 'IBD2', 3.0, 4.0),
                ('C', 'D', '1', 'IBD1', 0.0, 3.0),  # a 6 cM run of two 3 cM pieces: kept
                ('C', 'D', '1', 'IBD2', 3.0, 6.0),
                ('E', 'F', '1', 'IBD1', 0.0, 60.0),
            ],
        ),
    )
    assert table[['ID1', 'ID2']].to_numpy().tolist() == [['C', 'D'], ['E', 'F']]
    kept = row(table, 'C', 'D')
    assert kept['IBIS_N_SEG'] == 1
    assert kept['IBIS_TOTAL_CM'] == pytest.approx(6.0)


def test_overlapping_pieces_are_rejected():
    overlapping = pieces([('A', 'B', '1', 'IBD1', 0.0, 10.0), ('A', 'B', '1', 'IBD2', 9.0, 12.0)])
    with pytest.raises(ValueError, match=r'start before the previous segment .* chromosome 1, overlapping by 1\.0'):
        build(overlapping)


def test_overlapping_x_pieces_join_one_run():
    # IBIS called one X stretch of two males as both IBD1 and IBD2; the gap to the next piece is
    # measured from the longer one's end.
    x = pieces(
        [
            ('A', 'B', 'X', 'IBD1', 106.0, 111.6),
            ('A', 'B', 'X', 'IBD2', 106.0, 111.8),
            ('A', 'B', 'X', 'IBD1', 111.7, 113.0),
        ],
    )
    table = build(segments([('A', 'B', '1', 'IBD1', 6.0)]), x)
    joined = row(table, 'A', 'B')
    assert joined['X_N_SEG'] == 1
    assert joined['X_TOTAL_CM'] == pytest.approx(7.0)


def test_join_into_runs_returns_one_row_per_run():
    runs = summary.join_into_runs(
        pieces(
            [
                ('B', 'A', '1', 'IBD2', 10.0, 14.0),
                ('B', 'A', '1', 'IBD1', 0.0, 10.0),
                ('B', 'A', '1', 'IBD1', 20.0, 30.0),
                ('B', 'A', '2', 'IBD1', 5.0, 8.0),
            ],
        ),
        JOIN_GAP_CM,
    )
    assert list(runs.columns) == list(summary.RUN_COLUMNS)
    assert runs[['ID1', 'ID2', 'CHROM']].to_numpy().tolist() == [['A', 'B', '1'], ['A', 'B', '1'], ['A', 'B', '2']]
    assert runs['CM_START'].tolist() == [0.0, 20.0, 5.0]
    assert runs['CM_END'].tolist() == [14.0, 30.0, 8.0]
    assert runs['LENGTH_CM'].tolist() == pytest.approx([14.0, 10.0, 3.0])
    assert runs['KINSHIP_CM'].tolist() == pytest.approx([10 / 4 + 4 / 2, 10 / 4, 3 / 4])


def test_join_into_runs_of_no_segments_is_empty():
    runs = summary.join_into_runs(NO_X, JOIN_GAP_CM)
    assert list(runs.columns) == list(summary.RUN_COLUMNS)
    assert len(runs) == 0


def test_x_pieces_are_joined_into_runs_too():
    autosomal = segments([('A', 'B', '1', 'IBD1', 100.0)])
    x = pieces(
        [
            ('A', 'B', 'X', 'IBD1', 0.0, 4.0),
            ('A', 'B', 'X', 'IBD2', 4.0, 6.0),  # the first two join into a 6 cM run
            ('A', 'B', 'X', 'IBD1', 30.0, 50.0),
            ('C', 'D', 'X', 'IBD1', 0.0, 3.0),  # a 3 cM run: below the floor
        ],
    )
    table = build(autosomal, x)
    assert len(table) == 1  # C-D has no autosomal or KING row, and an X-only pair is not listed
    pair = row(table, 'A', 'B')
    assert pair['X_N_SEG'] == 2
    assert pair['X_TOTAL_CM'] == pytest.approx(26.0)
    assert pair['X_LONGEST_CM'] == pytest.approx(20.0)


def test_cutoff_below_min_seg_cm_is_rejected():
    with pytest.raises(ValueError, match='below min_seg_cm'):
        build(segments([('A', 'B', '1', 'IBD1', 60.0)]), cutoffs=[2.0, 10.0])


def test_repeated_cutoff_is_rejected():
    with pytest.raises(ValueError, match='repeat'):
        build(segments([('A', 'B', '1', 'IBD1', 60.0)]), cutoffs=[10.0, 10.0])


def test_non_positive_map_length_is_rejected():
    with pytest.raises(ValueError, match='map length'):
        summary.build_relatedness_table(NO_X, NO_X, NO_KING, 0.0, CUTOFFS, MIN_SEG_CM, JOIN_GAP_CM)


def test_fractional_cutoff_label_and_columns():
    table = build(segments([('A', 'B', '1', 'IBD1', 60.0)]), cutoffs=[7.5])
    assert row(table, 'A', 'B')['IBIS_KINSHIP_GE7.5CM'] == pytest.approx(60 / 4 / MAP_LENGTH_CM)


def test_no_segments_gives_king_only_rows():
    table = build(NO_X, king_seg=king([('A', 'B', 0.5, 0.25, 'FS')]))
    assert len(table) == 1
    assert table.iloc[0]['KING_INFTYPE'] == 'FS'
    assert pd.isna(table.iloc[0]['IBIS_N_SEG'])


def test_table_counts():
    table = build(
        segments([('A', 'B', '1', 'IBD1', 100.0)]),
        segments([('A', 'B', 'X', 'IBD1', 40.0)]),
        king_seg=king([('C', 'D', 0.5, 0.25, 'FS'), ('A', 'B', 0.3, 0.0, '3rd')]),
    )
    counts = summary.table_counts(table, CUTOFFS)
    assert counts['n_pairs'] == 2
    assert counts['n_king_pairs'] == 2
    assert counts['n_ibis_pairs'] == 1
    assert counts['n_king_pairs_without_ibis_segments'] == 1
    assert counts['n_pairs_with_x_segments'] == 1
    assert counts['n_pairs_by_ibis_degree']['GE5CM'] == {'3rd': 1}
    assert counts['n_pairs_by_ibis_degree']['GE15CM'] == {'3rd': 1}


# -- IBIS IDs: the FID prefix -----------------------------------------------------------------


def test_fid_prefix_is_stripped():
    ids = pd.Series(['0:SAMPLE_A', '0:SAMPLE_B'])
    assert summary.strip_fid_prefix(ids, '0').tolist() == ['SAMPLE_A', 'SAMPLE_B']


def test_id_without_the_expected_fid_prefix_is_rejected():
    with pytest.raises(ValueError, match='do not start with "0:"'):
        summary.strip_fid_prefix(pd.Series(['0:SAMPLE_A', 'SAMPLE_B']), '0')
    with pytest.raises(ValueError, match='do not start with "0:"'):
        summary.strip_fid_prefix(pd.Series(['1:SAMPLE_A']), '0')


def test_ids_not_in_the_fam_are_rejected():
    summary.check_ids_known(pd.Series(['SAMPLE_A', 'SAMPLE_B']), {'SAMPLE_A', 'SAMPLE_B', 'SAMPLE_C'}, 'KING')
    with pytest.raises(ValueError, match=r'1 KING sample IDs are not in the \.fam'):
        summary.check_ids_known(pd.Series(['SAMPLE_A', 'SAMPLE_Z']), {'SAMPLE_A', 'SAMPLE_B'}, 'KING')


def test_single_fid(tmp_path: Path):
    fam = tmp_path / 'a.fam'
    fam.write_text('0 SAMPLE_A 0 0 1 -9\n0 SAMPLE_B 0 0 2 -9\n')
    assert summary.single_fid(summary.read_fam(fam)) == '0'
    fam.write_text('0 SAMPLE_A 0 0 1 -9\nF2 SAMPLE_B 0 0 2 -9\n')
    with pytest.raises(ValueError, match='one FID'):
        summary.single_fid(summary.read_fam(fam))


# -- reading IBIS files -------------------------------------------------------------------------

IBIS_LINE = (
    '0:SAMPLE_A\t0:SAMPLE_B\t1\t1510000\t201500000\tIBD2\t0.002537\t195.556442\t195.553909\t20000\t0\t0.000000\n'
)


def test_read_ibis_segments_plain_and_gzip(tmp_path: Path):
    plain = tmp_path / 'ibis.seg'
    plain.write_text(IBIS_LINE + IBIS_LINE.replace('IBD2', 'IBD1'))
    zipped = tmp_path / 'ibis.seg.gz'
    with gzip.open(zipped, 'wt') as handle:
        handle.write(IBIS_LINE)

    both = summary.read_ibis_segments(plain)
    assert both['IBD_TYPE'].tolist() == ['IBD2', 'IBD1']
    assert both['CM_START'].tolist() == pytest.approx([0.002537, 0.002537])
    assert both['CM_END'].tolist() == pytest.approx([195.556442, 195.556442])
    assert both['LENGTH_CM'].tolist() == pytest.approx([195.553909, 195.553909])
    assert both['ID1'].tolist() == ['0:SAMPLE_A', '0:SAMPLE_A']
    assert both['CHROM'].tolist() == ['1', '1']
    assert len(summary.read_ibis_segments(zipped)) == 1


def test_read_ibis_segments_empty_file(tmp_path: Path):
    path = tmp_path / 'empty.seg'
    path.write_text('')
    empty = summary.read_ibis_segments(path)
    assert list(empty.columns) == list(summary.IBIS_SEG_USED_COLUMNS)
    assert len(empty) == 0


def test_read_ibis_segments_rejects_wrong_layout(tmp_path: Path):
    path = tmp_path / 'bad.seg'
    path.write_text('0:SAMPLE_A\t0:SAMPLE_B\t1\t1510000\n')
    with pytest.raises(ValueError, match='columns'):
        summary.read_ibis_segments(path)
    path.write_text(IBIS_LINE.replace('IBD2', 'HBD'))
    with pytest.raises(ValueError, match='segment type'):
        summary.read_ibis_segments(path)


def test_ibis_map_length_log_check():
    log = 'Defining Windows... done.\nTotal Genetic Length in use:3391.234375\nOrganizing genotype data\n'
    ibis_length = summary.parse_ibis_map_length(log)
    assert ibis_length == pytest.approx(3391.234375)
    summary.check_map_length_matches_ibis(3391.2345, ibis_length)
    with pytest.raises(ValueError, match='differs from the length IBIS used'):
        summary.check_map_length_matches_ibis(3394.29, ibis_length)  # a map that differs by ~3 cM
    with pytest.raises(ValueError, match='no "Total Genetic Length'):
        summary.parse_ibis_map_length('nothing here')


# -- .bim cM shift and checks ---------------------------------------------------------------------


def bim(rows: list[tuple[str, float]]) -> pd.DataFrame:
    """A .bim from (CHROM, CM) rows, with positions increasing within each chromosome."""
    return pd.DataFrame(
        {
            'CHROM': [chromosome for chromosome, _ in rows],
            'ID': [f'v{i}' for i in range(len(rows))],
            'CM': [cm for _, cm in rows],
            'POS': [1000 * (i + 1) for i in range(len(rows))],
            'A1': 'A',
            'A2': 'G',
        },
    )


def test_negative_start_is_shifted_up_per_chromosome_and_distances_are_kept():
    original = bim([('1', -0.5), ('1', -0.1), ('1', 2.0), ('2', 0.25), ('2', 3.0), ('3', -0.2), ('3', 4.0)])
    shifted = summary.shift_to_non_negative_cm(original)
    assert shifted['CM'].tolist() == pytest.approx([0.0, 0.4, 2.5, 0.25, 3.0, 0.0, 4.2])
    assert original['CM'].iloc[0] == -0.5  # the input is unchanged
    for chromosome in ('1', '2', '3'):
        before = original[original['CHROM'] == chromosome]['CM'].diff().dropna()
        after = shifted[shifted['CHROM'] == chromosome]['CM'].diff().dropna()
        assert after.tolist() == pytest.approx(before.tolist())
    assert summary.cm_shift_by_chromosome(original).to_dict() == pytest.approx({'1': 0.5, '2': 0.0, '3': 0.2})


def test_clamping_would_leave_ibis_without_a_map_but_the_shift_does_not():
    """The real chromosome 1 starts with two variants before the map's first row, both below 0 cM."""
    start = bim([('1', -0.026), ('1', -0.020), ('1', 0.5), ('1', 90.0)])
    clamped = start.assign(CM=start['CM'].clip(lower=0.0))
    with pytest.raises(ValueError, match='no genetic map'):
        summary.check_ibis_reads_map(clamped)
    shifted, length = summary.finalise_bim(start, ['1'])
    assert shifted['CM'].iloc[:2].tolist() == pytest.approx([0.0, 0.006])
    assert length == pytest.approx(90.026)  # last minus first, exactly as unshifted


def test_ibis_map_check_uses_float32_like_ibis():
    summary.check_ibis_reads_map(bim([('1', 0.0), ('1', 0.002)]))
    with pytest.raises(ValueError, match='no genetic map'):
        summary.check_ibis_reads_map(bim([('1', 0.0), ('1', 1e-50)]))  # underflows to 0 in float32


def test_cm_decrease_within_a_chromosome_is_rejected():
    summary.check_cm_non_decreasing(
        bim([('1', 0.0), ('1', 1.0), ('1', 1.0), ('2', 0.0), ('2', 3.0)])
    )  # a plateau is fine
    with pytest.raises(ValueError, match='decreases at variant line 3 on chromosome 1'):
        summary.check_cm_non_decreasing(bim([('1', 0.0), ('1', 2.0), ('1', 1.5)]))


def test_finalise_bim_rejects_a_map_that_decreases():
    with pytest.raises(ValueError, match='decreases'):
        summary.finalise_bim(bim([('X', 0.0), ('X', 5.0), ('X', 4.0)]), ['X'])


def test_map_length_is_the_sum_of_chromosome_spans():
    table = bim([('1', 0.0), ('1', 100.0), ('2', 0.5), ('2', 80.5), ('3', 10.0), ('3', 10.0)])
    assert summary.map_length_cm(table) == pytest.approx(100.0 + 80.0 + 0.0)


def test_finalise_bim_measures_the_corrected_map():
    rows = [(str(c), cm) for c in summary.AUTOSOMES for cm in (-0.3, 10.0, 50.0)]
    shifted, length = summary.finalise_bim(bim(rows), summary.AUTOSOMES)
    assert (shifted['CM'] >= 0).all()
    assert length == pytest.approx(22 * 50.3)  # unchanged by the shift; clamping would give 22 * 50.0


def test_chromosome_names_must_match_what_ibis_expects():
    summary.check_chromosome_blocks(bim([('X', 0.0), ('X', 1.0)]), ['X'])
    for wrong in ('chrX', '23'):
        with pytest.raises(ValueError, match='chromosome blocks'):
            summary.check_chromosome_blocks(bim([(wrong, 0.0), (wrong, 1.0)]), ['X'])
    with pytest.raises(ValueError, match='chromosome blocks'):  # chromosome 2 missing
        summary.check_chromosome_blocks(bim([('1', 0.0), ('3', 0.0)]), ['1', '2', '3'])
    with pytest.raises(ValueError, match='chromosome blocks'):  # chromosome 1 split in two blocks
        summary.check_chromosome_blocks(bim([('1', 0.0), ('2', 0.0), ('1', 5.0)]), ['1', '2'])


def test_a_chromosome_without_a_map_is_rejected():
    summary.check_map_applied(bim([('1', 0.0), ('1', 5.0), ('2', 1.0), ('2', 6.0)]))
    with pytest.raises(ValueError, match=r'\[.X.\] span 0 cM'):
        summary.check_map_applied(bim([('1', 0.0), ('1', 5.0), ('X', 0.0), ('X', 0.0)]))


def test_bim_round_trip(tmp_path: Path):
    path = tmp_path / 'a.bim'
    path.write_text('1\tv1\t0.002054\t2000\tB\tA\n1\tv2\t-0.5\t3000\tNA\tG\n')
    table = summary.read_bim(path)
    assert table['A1'].tolist() == ['B', 'NA']  # an allele literally called NA is not read as missing
    out = tmp_path / 'b.bim'
    summary.write_bim(table.assign(CM=[0.0, 1.5]), out)
    assert out.read_text() == '1\tv1\t0\t2000\tB\tA\n1\tv2\t1.5\t3000\tNA\tG\n'


# -- the two job-side steps end to end ---------------------------------------------------------------


def test_prepare_maps_and_summarise_steps(tmp_path: Path):
    autosomal_rows = [
        f'{c}\tv{c}_{i}\t{cm}\t{1000 * (i + 1)}\tA\tG'
        for c in summary.AUTOSOMES
        for i, cm in enumerate((-0.1, 0.4, 20.0))
    ]
    (tmp_path / 'auto.bim').write_text('\n'.join(autosomal_rows) + '\n')
    (tmp_path / 'x.bim').write_text('X\tvx_0\t-0.2\t1000\tA\tG\nX\tvx_1\t30.0\t2000\tA\tG\n')
    summary.main(
        [
            'prepare-maps',
            '--auto-bim', str(tmp_path / 'auto.bim'),
            '--x-bim', str(tmp_path / 'x.bim'),
            '--out-auto-bim', str(tmp_path / 'auto_fixed.bim'),
            '--out-x-bim', str(tmp_path / 'x_fixed.bim'),
            '--out-stats', str(tmp_path / 'stats.json'),
        ],
    )  # fmt: skip
    stats = json.loads((tmp_path / 'stats.json').read_text())
    assert stats['autosomal_map_length_cm'] == pytest.approx(22 * 20.1)
    assert stats['x_map_length_cm'] == pytest.approx(30.2)
    assert stats['n_autosomal_variants_before_map_start'] == 22
    assert stats['n_x_variants_before_map_start'] == 1
    assert stats['autosomal_max_cm_shift'] == pytest.approx(0.1)
    assert stats['x_cm_shift'] == pytest.approx(0.2)
    assert (tmp_path / 'x_fixed.bim').read_text() == 'X\tvx_0\t0\t1000\tA\tG\nX\tvx_1\t30.2\t2000\tA\tG\n'

    (tmp_path / 'auto.fam').write_text('0 SAMPLE_A 0 0 1 -9\n0 SAMPLE_B 0 0 2 -9\n0 SAMPLE_C 0 0 2 -9\n')
    (tmp_path / 'x.fam').write_text('0 SAMPLE_A 0 0 1 -9\n0 SAMPLE_B 0 0 2 -9\n')
    (tmp_path / 'ibis_auto.log').write_text(f'Total Genetic Length in use:{22 * 20.1:f}\n')
    seg_line = '0:{a}\t0:{b}\t1\t1000\t2000\tIBD1\t0.0\t{cm}\t{cm}\t500\t0\t0.0\n'
    ibd2_line = '0:SAMPLE_B\t0:SAMPLE_A\t1\t2000\t3000\tIBD2\t60.0\t70.0\t10.0\t500\t0\t0.0\n'  # abuts the 60 cM IBD1
    (tmp_path / 'ibis_auto.seg').write_text(seg_line.format(a='SAMPLE_B', b='SAMPLE_A', cm=60.0) + ibd2_line)
    (tmp_path / 'ibis_x.seg').write_text(seg_line.format(a='SAMPLE_A', b='SAMPLE_B', cm=12.0).replace('\t1\t', '\tX\t'))
    (tmp_path / 'king.seg').write_text(
        'FID1\tID1\tFID2\tID2\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType\n0\tSAMPLE_A\t0\tSAMPLE_C\t0.3\t0.0\t0.15\t3rd\n',
    )
    summary.main(
        [
            'summarise',
            '--ibis-auto-seg', str(tmp_path / 'ibis_auto.seg'),
            '--ibis-x-seg', str(tmp_path / 'ibis_x.seg'),
            '--ibis-auto-log', str(tmp_path / 'ibis_auto.log'),
            '--king-seg', str(tmp_path / 'king.seg'),
            '--fam', str(tmp_path / 'auto.fam'),
            '--fam-x', str(tmp_path / 'x.fam'),
            '--map-stats', str(tmp_path / 'stats.json'),
            '--min-seg-cm', '5',
            '--join-gap-cm', '0.5',
            '--cutoffs-cm', '5,10',
            '--run-info', json.dumps({'inputs': {'king_seg': 'gs://bucket/king.seg'}}),
            '--out-table', str(tmp_path / 'pairs.tsv.gz'),
            '--out-sidecar', str(tmp_path / 'pairs.json'),
        ],
    )  # fmt: skip
    table = pd.read_csv(tmp_path / 'pairs.tsv.gz', sep='\t')
    assert table[['ID1', 'ID2']].to_numpy().tolist() == [['SAMPLE_A', 'SAMPLE_B'], ['SAMPLE_A', 'SAMPLE_C']]
    ab = table.iloc[0]
    assert ab['IBIS_N_SEG'] == 1  # the IBD1 and IBD2 segments are joined into one 70 cM run
    assert ab['IBIS_TOTAL_CM'] == pytest.approx(70.0)
    expected_kinship = (60 / 4 + 10 / 2) / (22 * 20.1)
    assert ab['IBIS_KINSHIP_GE10CM'] == pytest.approx(expected_kinship, abs=1e-6)  # the TSV keeps six decimals
    assert ab['X_TOTAL_CM'] == pytest.approx(12.0)
    assert table.iloc[1]['KING_INFTYPE'] == '3rd'
    sidecar = json.loads((tmp_path / 'pairs.json').read_text())
    assert sidecar['autosomal_map_length_cm'] == pytest.approx(22 * 20.1)
    assert sidecar['min_seg_cm'] == 5.0
    assert sidecar['join_gap_cm'] == 0.5
    assert sidecar['run']['inputs']['king_seg'] == 'gs://bucket/king.seg'
    assert sidecar['n_samples'] == 3
    assert sidecar['n_samples_x'] == 2


def test_summarise_rejects_king_ids_missing_from_the_fam(tmp_path: Path):
    (tmp_path / 'stats.json').write_text(json.dumps({'autosomal_map_length_cm': 440.0, 'x_map_length_cm': 30.0}))
    (tmp_path / 'auto.fam').write_text('0 SAMPLE_A 0 0 1 -9\n0 SAMPLE_B 0 0 2 -9\n')
    (tmp_path / 'ibis_auto.log').write_text('Total Genetic Length in use:440.000000\n')
    (tmp_path / 'empty.seg').write_text('')
    (tmp_path / 'king.seg').write_text(
        'FID1\tID1\tFID2\tID2\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType\n0\tSAMPLE_A\t0\tSAMPLE_Z\t0.3\t0.0\t0.15\t3rd\n',
    )
    args = Namespace(
        ibis_auto_seg=str(tmp_path / 'empty.seg'),
        ibis_x_seg=str(tmp_path / 'empty.seg'),
        ibis_auto_log=str(tmp_path / 'ibis_auto.log'),
        king_seg=str(tmp_path / 'king.seg'),
        fam=str(tmp_path / 'auto.fam'),
        fam_x=str(tmp_path / 'auto.fam'),
        map_stats=str(tmp_path / 'stats.json'),
        min_seg_cm=5.0,
        join_gap_cm=0.5,
        cutoffs_cm='5,10',
        run_info='{}',
        out_table=str(tmp_path / 'pairs.tsv.gz'),
        out_sidecar=str(tmp_path / 'pairs.json'),
    )
    with pytest.raises(ValueError, match=r'1 KING sample IDs are not in the \.fam'):
        summary.summarise(args)


def test_summarise_rejects_ibis_ids_missing_from_the_fam(tmp_path: Path):
    """IBIS IDs are `<FID>:<IID>`; after stripping they must all be samples IBIS was given."""
    (tmp_path / 'stats.json').write_text(json.dumps({'autosomal_map_length_cm': 440.0, 'x_map_length_cm': 30.0}))
    (tmp_path / 'auto.fam').write_text('0 SAMPLE_A 0 0 1 -9\n0 SAMPLE_B 0 0 2 -9\n')
    (tmp_path / 'ibis_auto.log').write_text('Total Genetic Length in use:440.000000\n')
    (tmp_path / 'empty.seg').write_text('')
    (tmp_path / 'ibis.seg').write_text(IBIS_LINE.replace('SAMPLE_B', 'SAMPLE_Z'))
    (tmp_path / 'king.seg').write_text('FID1\tID1\tFID2\tID2\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType\n')
    args = Namespace(
        ibis_auto_seg=str(tmp_path / 'ibis.seg'),
        ibis_x_seg=str(tmp_path / 'empty.seg'),
        ibis_auto_log=str(tmp_path / 'ibis_auto.log'),
        king_seg=str(tmp_path / 'king.seg'),
        fam=str(tmp_path / 'auto.fam'),
        fam_x=str(tmp_path / 'auto.fam'),
        map_stats=str(tmp_path / 'stats.json'),
        min_seg_cm=5.0,
        join_gap_cm=0.5,
        cutoffs_cm='5,10',
        run_info='{}',
        out_table=str(tmp_path / 'pairs.tsv.gz'),
        out_sidecar=str(tmp_path / 'pairs.json'),
    )
    with pytest.raises(ValueError, match=r'1 IBIS autosomal sample IDs are not in the \.fam'):
        summary.summarise(args)
