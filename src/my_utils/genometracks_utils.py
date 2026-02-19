"""
Utilities for working with UCSC Genome Browser track data.

Functions for loading gzip-compressed JSON exports from the UCSC Genome Browser
(bigWig signal tracks) and computing aggregate signals over BED12 blocked regions.

Typical workflow:
    1. parse_ucsc_tracks() — load JSON
    2. get_block_signal_df() — intersect tracks with BED12 blocks → per-interval DataFrame
    3. Aggregate the DataFrame however you like (groupby + weighted mean, etc.)
"""

import gzip
import json
from bisect import bisect_left, bisect_right

import numpy as np
import pandas as pd

from my_utils.bed_utils import merge_blocks


def parse_ucsc_tracks(path):
    """
    Load a gzip-compressed JSON file exported from the UCSC Genome Browser.

    Args:
        path: Path to the .json.gz file

    Returns:
        Dict with two keys:
        - 'metadata': dict with genome, chrom, region start/end, trackType, etc.
        - 'tracks': dict mapping track_name -> {chrom: [(start, end, value), ...]}
    """
    with gzip.open(path, "rt") as f:
        data = json.load(f)

    metadata = {}
    tracks = {}

    for key, value in data.items():
        # Track data is a dict mapping chrom -> list of [start, end, value]
        if (isinstance(value, dict)
                and all(isinstance(v, list) for v in value.values())):
            track = {}
            for chrom, intervals in value.items():
                track[chrom] = [(iv[0], iv[1], iv[2]) for iv in intervals]
            tracks[key] = track
        else:
            metadata[key] = value

    return {"metadata": metadata, "tracks": tracks}


def get_bed12_block_coords(bed12_line):
    """
    Parse a BED12 line and return absolute genomic coordinates for each block.

    Args:
        bed12_line: A string (tab-delimited BED12 line) or list of fields

    Returns:
        List of (start, end) tuples in absolute genomic coordinates
    """
    if isinstance(bed12_line, str):
        fields = bed12_line.strip().split("\t")
    else:
        fields = bed12_line

    chrom_start = int(fields[1])
    block_count = int(fields[9])
    block_sizes = [int(x) for x in fields[10].rstrip(",").split(",")]
    block_starts = [int(x) for x in fields[11].rstrip(",").split(",")]

    blocks = []
    for i in range(block_count):
        abs_start = chrom_start + block_starts[i]
        abs_end = abs_start + block_sizes[i]
        blocks.append((abs_start, abs_end))

    return blocks


def _lookup_value(intervals, interval_starts, pos):
    """Look up the signal value at a single position using binary search.

    Returns the value if pos falls within an interval, else NaN.
    """
    idx = bisect_right(interval_starts, pos) - 1
    if idx >= 0:
        iv_start, iv_end, iv_value = intervals[idx]
        if iv_start <= pos < iv_end:
            return iv_value
    return np.nan


def _intersect_block_intervals(intervals, interval_starts, block_start, block_end):
    """Yield (start, end, value) sub-intervals where a single block overlaps signal.

    Also yields gap sub-intervals with value=NaN for block regions with no signal.
    """
    n = len(intervals)
    # Find candidate intervals overlapping [block_start, block_end)
    left = bisect_right(interval_starts, block_start) - 1
    if left < 0:
        left = 0
    right = bisect_left(interval_starts, block_end)

    cursor = block_start

    for i in range(left, min(right, n)):
        iv_start, iv_end, iv_value = intervals[i]
        overlap_start = max(iv_start, block_start)
        overlap_end = min(iv_end, block_end)

        if overlap_start >= overlap_end:
            continue

        # Gap before this overlap
        if cursor < overlap_start:
            yield (cursor, overlap_start, np.nan)

        yield (overlap_start, overlap_end, iv_value)
        cursor = overlap_end

    # Trailing gap
    if cursor < block_end:
        yield (cursor, block_end, np.nan)


def get_block_signal_df(tracks_path, bed12_path, track_names=None, padding=0):
    """
    Build a per-sub-interval DataFrame of signal values across tracks and BED12 blocks.

    For each BED12 entry's blocks, computes the union of interval breakpoints across
    all tracks to produce sub-intervals. Each row is one sub-interval with columns
    for every track's signal value (NaN where the track has no data).

    Args:
        tracks_path: Path to gzip-compressed UCSC JSON file
        bed12_path: Path to BED12 file
        track_names: List of track names to include, or None for all tracks
        padding: Number of bases to extend each block on both sides. Padded
            blocks are merged when they overlap. Sub-intervals within original
            (unpadded) blocks get ``in_block=True``; flanking padding gets
            ``in_block=False``.

    Returns:
        pandas DataFrame with columns:
        - bed_name: parent BED12 feature name
        - chrom: chromosome
        - block_id: ordered merged-block label ("block_001", ...)
        - start, end: sub-interval coordinates
        - width: end - start
        - in_block: True if sub-interval is within an original (unpadded) block
        - one column per track: signal value (or NaN)
    """
    if not isinstance(padding, int) or padding < 0:
        raise ValueError(f"padding must be an int >= 0, got {padding!r}")

    parsed = parse_ucsc_tracks(tracks_path)
    tracks = parsed["tracks"]

    if track_names is not None:
        tracks = {k: v for k, v in tracks.items() if k in track_names}

    track_name_list = list(tracks.keys())

    # Pre-compute interval_starts lookup arrays per track per chrom
    track_interval_starts = {}
    for tname, tdata in tracks.items():
        track_interval_starts[tname] = {}
        for chrom, intervals in tdata.items():
            track_interval_starts[tname][chrom] = [iv[0] for iv in intervals]

    rows = []

    with open(bed12_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 12:
                continue

            chrom = fields[0]
            bed_name = fields[3]
            original_blocks = get_bed12_block_coords(fields)

            # Pad and merge blocks
            padded = [(max(0, s - padding), e + padding) for s, e in original_blocks]
            merged = merge_blocks(padded)

            for block_idx, (merged_start, merged_end) in enumerate(merged, start=1):
                block_id = f"block_{block_idx:03d}"

                # Collect all breakpoints across all tracks for this merged block
                breakpoints = {merged_start, merged_end}
                # Also add original block boundaries so sub-intervals align
                for orig_s, orig_e in original_blocks:
                    if orig_s >= merged_start and orig_s <= merged_end:
                        breakpoints.add(orig_s)
                    if orig_e >= merged_start and orig_e <= merged_end:
                        breakpoints.add(orig_e)

                for tname in track_name_list:
                    intervals = tracks[tname].get(chrom, [])
                    iv_starts = track_interval_starts[tname].get(chrom, [])
                    for _start, _end, _val in _intersect_block_intervals(
                            intervals, iv_starts, merged_start, merged_end):
                        breakpoints.add(_start)
                        breakpoints.add(_end)

                breakpoints = sorted(breakpoints)

                # Build rows for each sub-interval
                for j in range(len(breakpoints) - 1):
                    sub_start = breakpoints[j]
                    sub_end = breakpoints[j + 1]
                    if sub_start >= sub_end:
                        continue

                    # Check if sub-interval is fully within any original block
                    in_block = any(
                        sub_start >= orig_s and sub_end <= orig_e
                        for orig_s, orig_e in original_blocks
                    )

                    row = {
                        "bed_name": bed_name,
                        "chrom": chrom,
                        "block_id": block_id,
                        "start": sub_start,
                        "end": sub_end,
                        "width": sub_end - sub_start,
                        "in_block": in_block,
                    }
                    for tname in track_name_list:
                        intervals = tracks[tname].get(chrom, [])
                        iv_starts = track_interval_starts[tname].get(chrom, [])
                        row[tname] = _lookup_value(intervals, iv_starts, sub_start)

                    rows.append(row)

    col_order = ["bed_name", "chrom", "block_id", "start", "end", "width",
                 "in_block"] + track_name_list
    return pd.DataFrame(rows, columns=col_order)


def _parse_region(region_str):
    """Parse 'chr2:1,850,700-1,886,852' → ('chr2', 1850700, 1886852)."""
    chrom, coords = region_str.split(":")
    start_str, end_str = coords.split("-")
    return chrom, int(start_str.replace(",", "")), int(end_str.replace(",", ""))


def plot_block_signal(df, output_path, track_names=None, region=None,
                      bed_names=None, ylim=None, colors=None):
    """
    Plot signal across blocks from a get_block_signal_df() DataFrame.

    Creates one row of subplots per bed_name, one column per block_id.
    Within each subplot, signal is drawn as a step plot. Regions where
    ``in_block=True`` are shaded to distinguish exonic from padding regions.

    Args:
        df: DataFrame from get_block_signal_df()
        output_path: Path to save the figure (e.g. "signal.pdf")
        track_names: List of signal columns to overlay, or None for all
            signal columns (auto-detected as columns after the 7 metadata cols)
        region: Genomic region string like "chr2:1,850,700-1,886,852".
            Blocks entirely outside the region are dropped.
        bed_names: List of bed_name values to plot, or None for all
        ylim: (min, max) tuple for shared y-axis, or None to auto-compute
        colors: List of colors (one per track), or None for matplotlib tab10/tab20
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    _meta_cols = {"bed_name", "chrom", "block_id", "start", "end", "width", "in_block"}

    if track_names is None:
        track_names = [c for c in df.columns if c not in _meta_cols]

    has_block_id = "block_id" in df.columns
    has_in_block = "in_block" in df.columns

    # Auto-detect blocks from gaps when block_id is missing
    if not has_block_id:
        df = df.copy()
        df["block_id"] = ""
        for _, grp in df.groupby("bed_name"):
            grp = grp.sort_values("start")
            idx = grp.index
            ends = grp["end"].values[:-1]
            starts = grp["start"].values[1:]
            block_num = 1
            nums = [block_num]
            for i in range(len(starts)):
                if starts[i] > ends[i]:
                    block_num += 1
                nums.append(block_num)
            df.loc[idx, "block_id"] = [f"block_{n:03d}" for n in nums]

    if bed_names is not None:
        df = df[df["bed_name"].isin(bed_names)]

    # Region filtering: drop blocks entirely outside the region
    if region is not None:
        r_chrom, r_start, r_end = _parse_region(region)
        keep_mask = pd.Series(False, index=df.index)
        for (_bn, _bid), grp in df.groupby(["bed_name", "block_id"]):
            block_chrom = grp["chrom"].iloc[0]
            block_min = grp["start"].min()
            block_max = grp["end"].max()
            if block_chrom == r_chrom and block_min < r_end and block_max > r_start:
                keep_mask.loc[grp.index] = True
        df = df[keep_mask]

    # Auto y-axis limits
    if ylim is None:
        all_vals = df[track_names].values.ravel()
        all_vals = all_vals[~np.isnan(all_vals)]
        if len(all_vals) > 0:
            ylim = (float(np.min(all_vals)), float(np.max(all_vals)))
        else:
            ylim = (0, 1)

    # Default colors from tab10/tab20
    if colors is None:
        if len(track_names) <= 10:
            cmap = plt.cm.tab10
        else:
            cmap = plt.cm.tab20
        colors = [cmap(i % cmap.N) for i in range(len(track_names))]

    bed_name_list = df["bed_name"].unique()

    with PdfPages(output_path) as pdf:
        for bed_name in bed_name_list:
            sub = df[df["bed_name"] == bed_name]
            block_ids = sub["block_id"].unique()
            n_blocks = len(block_ids)

            fig, axes = plt.subplots(
                1, n_blocks, figsize=(4 * n_blocks, 3), squeeze=False
            )
            fig.suptitle(bed_name, fontsize=10)

            for col_idx, bid in enumerate(block_ids):
                ax = axes[0, col_idx]
                bdf = sub[sub["block_id"] == bid].sort_values("start")

                # Overlay each track
                for t_idx, tname in enumerate(track_names):
                    xs, ys = [], []
                    for _, row in bdf.iterrows():
                        xs.extend([row["start"], row["end"]])
                        ys.extend([row[tname], row[tname]])
                    ax.plot(xs, ys, linewidth=0.8, color=colors[t_idx],
                            label=tname if col_idx == 0 else None)

                # Shade in_block regions
                if has_in_block:
                    for _, row in bdf[bdf["in_block"]].iterrows():
                        ax.axvspan(
                            row["start"], row["end"],
                            alpha=0.15, color="orange", linewidth=0,
                        )

                ax.set_ylim(ylim)
                ax.set_title(bid, fontsize=8)
                ax.set_xlabel("position", fontsize=7)
                ax.tick_params(labelsize=6)
                if col_idx == 0:
                    ax.set_ylabel("signal", fontsize=7)

            fig.legend(loc="upper right", fontsize=6)
            fig.tight_layout(rect=[0, 0, 1, 0.93])
            pdf.savefig(fig)
            plt.close(fig)


def plot_block_signal_cli():
    """Command-line tool: compute block signal from UCSC JSON + BED12, write TSV and PDF."""
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Compute per-block signal from a UCSC Genome Browser JSON export "
                    "and a BED12 file, then write a TSV table and/or a PDF of plots."
    )
    parser.add_argument("tracks_json", help="Path to gzip-compressed UCSC JSON file")
    parser.add_argument("bed12", help="Path to BED12 file")
    parser.add_argument("-o", "--output-prefix", required=True,
                        help="Output prefix (writes PREFIX.tsv and PREFIX.pdf)")
    parser.add_argument("-t", "--track-names", nargs="+", default=None,
                        help="Track names to include (default: all)")
    parser.add_argument("-p", "--padding", type=int, default=0,
                        help="Bases to pad each block (default: 0)")
    parser.add_argument("-r", "--region", default=None,
                        help="Genomic region to filter, e.g. 'chr2:1,850,700-1,886,852'")
    parser.add_argument("-b", "--bed-names", nargs="+", default=None,
                        help="BED feature names to include (default: all)")
    parser.add_argument("--ylim", nargs=2, type=float, default=None, metavar=("MIN", "MAX"),
                        help="Y-axis limits for plot (default: auto)")
    parser.add_argument("--no-tsv", action="store_true",
                        help="Skip writing the TSV file")
    parser.add_argument("--no-pdf", action="store_true",
                        help="Skip writing the PDF file")
    args = parser.parse_args()

    if args.no_tsv and args.no_pdf:
        parser.error("Cannot use both --no-tsv and --no-pdf")

    df = get_block_signal_df(
        args.tracks_json, args.bed12,
        track_names=args.track_names, padding=args.padding,
    )

    tsv_path = f"{args.output_prefix}.tsv"
    pdf_path = f"{args.output_prefix}.pdf"

    if not args.no_tsv:
        df.to_csv(tsv_path, sep="\t", index=False)
        print(f"Wrote {tsv_path}", file=sys.stderr)

    if not args.no_pdf:
        ylim_tuple = tuple(args.ylim) if args.ylim else None
        plot_block_signal(
            df, pdf_path,
            track_names=args.track_names, region=args.region,
            bed_names=args.bed_names, ylim=ylim_tuple,
        )
        print(f"Wrote {pdf_path}", file=sys.stderr)
