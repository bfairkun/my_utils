"""
Tests for genometracks_utils module.
"""

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from my_utils.genometracks_utils import (
    _intersect_block_intervals,
    _lookup_value,
    get_bed12_block_coords,
    get_block_signal_df,
)


class TestGetBed12BlockCoords:
    """Tests for get_bed12_block_coords."""

    def test_single_block(self):
        line = "chr1\t1000\t2000\tgene1\t0\t+\t1000\t2000\t0\t1\t1000,\t0,"
        blocks = get_bed12_block_coords(line)
        assert blocks == [(1000, 2000)]

    def test_two_blocks(self):
        line = "chr1\t1000\t2000\tgene2\t0\t+\t1000\t2000\t0\t2\t200,500,\t0,500,"
        blocks = get_bed12_block_coords(line)
        assert blocks == [(1000, 1200), (1500, 2000)]

    def test_three_blocks(self):
        line = "chr2\t5000\t8000\tgene3\t0\t-\t5000\t8000\t0\t3\t100,200,150,\t0,1000,2850,"
        blocks = get_bed12_block_coords(line)
        assert blocks == [(5000, 5100), (6000, 6200), (7850, 8000)]

    def test_accepts_list_of_fields(self):
        fields = [
            "chr1",
            "1000",
            "2000",
            "gene1",
            "0",
            "+",
            "1000",
            "2000",
            "0",
            "1",
            "1000,",
            "0,",
        ]
        blocks = get_bed12_block_coords(fields)
        assert blocks == [(1000, 2000)]

    def test_no_trailing_comma(self):
        line = "chr1\t1000\t2000\tgene2\t0\t+\t1000\t2000\t0\t2\t200,500\t0,500"
        blocks = get_bed12_block_coords(line)
        assert blocks == [(1000, 1200), (1500, 2000)]


class TestLookupValue:
    """Tests for _lookup_value."""

    def test_hit(self):
        intervals = [(100, 200, 5.0), (200, 300, 3.0)]
        starts = [100, 200]
        assert _lookup_value(intervals, starts, 150) == 5.0
        assert _lookup_value(intervals, starts, 200) == 3.0

    def test_miss(self):
        intervals = [(100, 200, 5.0)]
        starts = [100]
        assert np.isnan(_lookup_value(intervals, starts, 50))
        assert np.isnan(_lookup_value(intervals, starts, 200))

    def test_boundary(self):
        intervals = [(100, 200, 5.0)]
        starts = [100]
        # start is inclusive
        assert _lookup_value(intervals, starts, 100) == 5.0
        # end is exclusive
        assert np.isnan(_lookup_value(intervals, starts, 200))


class TestIntersectBlockIntervals:
    """Tests for _intersect_block_intervals."""

    def test_full_coverage(self):
        intervals = [(100, 200, 1.0), (200, 300, 2.0)]
        starts = [100, 200]
        result = list(_intersect_block_intervals(intervals, starts, 100, 300))
        assert result == [(100, 200, 1.0), (200, 300, 2.0)]

    def test_partial_overlap(self):
        intervals = [(100, 200, 1.0)]
        starts = [100]
        result = list(_intersect_block_intervals(intervals, starts, 150, 250))
        assert len(result) == 2
        assert result[0] == (150, 200, 1.0)
        # gap
        assert result[1][0] == 200
        assert result[1][1] == 250
        assert np.isnan(result[1][2])

    def test_no_overlap(self):
        intervals = [(100, 200, 1.0)]
        starts = [100]
        result = list(_intersect_block_intervals(intervals, starts, 300, 400))
        assert len(result) == 1
        assert result[0][0] == 300
        assert result[0][1] == 400
        assert np.isnan(result[0][2])

    def test_gap_between_intervals(self):
        intervals = [(100, 150, 1.0), (180, 200, 2.0)]
        starts = [100, 180]
        result = list(_intersect_block_intervals(intervals, starts, 100, 200))
        assert result[0] == (100, 150, 1.0)
        assert result[1][0] == 150
        assert result[1][1] == 180
        assert np.isnan(result[1][2])
        assert result[2] == (180, 200, 2.0)


@pytest.fixture
def _signal_fixture(tmp_path):
    """Create minimal track JSON and BED12 files for get_block_signal_df tests.

    BED12 entry has two blocks: (1000,1200) and (1500,2000) on chr1.
    Signal track covers 900-2100 with value 1.0.
    """
    track_data = {
        "genome": "hg38",
        "sig": {"chr1": [[900, 2100, 1.0]]},
    }
    tracks_path = tmp_path / "tracks.json.gz"
    with gzip.open(tracks_path, "wt") as f:
        json.dump(track_data, f)

    bed_path = tmp_path / "test.bed"
    bed_path.write_text(
        "chr1\t1000\t2000\tgene1\t0\t+\t1000\t2000\t0\t2\t200,500\t0,500\n"
    )
    return str(tracks_path), str(bed_path)


class TestGetBlockSignalDfPadding:
    """Tests for padding parameter of get_block_signal_df."""

    def test_padding_zero_in_block_all_true(self, _signal_fixture):
        tracks_path, bed_path = _signal_fixture
        df = get_block_signal_df(tracks_path, bed_path, padding=0)
        assert df["in_block"].all()

    def test_padding_zero_sequential_block_ids(self, _signal_fixture):
        tracks_path, bed_path = _signal_fixture
        df = get_block_signal_df(tracks_path, bed_path, padding=0)
        block_ids = df["block_id"].unique().tolist()
        assert block_ids == ["block_001", "block_002"]

    def test_padding_positive_flanks_not_in_block(self, _signal_fixture):
        tracks_path, bed_path = _signal_fixture
        df = get_block_signal_df(tracks_path, bed_path, padding=50)
        # Flanking rows should have in_block=False
        flank_rows = df[~df["in_block"]]
        assert len(flank_rows) > 0
        # All flank rows should be outside the original blocks
        for _, row in flank_rows.iterrows():
            inside = (row["start"] >= 1000 and row["end"] <= 1200) or (
                row["start"] >= 1500 and row["end"] <= 2000
            )
            assert not inside

    def test_padding_merges_adjacent_blocks(self, _signal_fixture):
        tracks_path, bed_path = _signal_fixture
        # Gap between blocks is 300bp (1200-1500). padding=150 makes them touch.
        df = get_block_signal_df(tracks_path, bed_path, padding=150)
        block_ids = df["block_id"].unique().tolist()
        assert block_ids == ["block_001"]

    def test_negative_padding_raises(self, _signal_fixture):
        tracks_path, bed_path = _signal_fixture
        with pytest.raises(ValueError, match="padding must be an int >= 0"):
            get_block_signal_df(tracks_path, bed_path, padding=-1)

    def test_float_padding_raises(self, _signal_fixture):
        tracks_path, bed_path = _signal_fixture
        with pytest.raises(ValueError, match="padding must be an int >= 0"):
            get_block_signal_df(tracks_path, bed_path, padding=1.5)

    def test_output_columns(self, _signal_fixture):
        tracks_path, bed_path = _signal_fixture
        df = get_block_signal_df(tracks_path, bed_path, padding=0)
        expected = [
            "bed_name",
            "chrom",
            "block_id",
            "start",
            "end",
            "width",
            "in_block",
            "sig",
        ]
        assert list(df.columns) == expected


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
