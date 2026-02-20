"""
Utilities for working with BED format files.

Simple conversion functions based on BedToBed12.py and Bed6_to_Bed12.py scripts.
"""

import sys
from collections import defaultdict


def bed_to_bed12(fields, defaults=None):
    """
    Convert BED3/4/6/9 to BED12 format with a single block.

    Args:
        fields: List of BED fields from a line
        defaults: Dict with optional default values for 'name', 'score', 'strand', 'color'

    Returns:
        List of BED12 fields, or None if invalid
    """
    if defaults is None:
        defaults = {}

    # Need at least chrom, start, end
    if len(fields) < 3:
        return None

    chrom = fields[0]
    try:
        start = int(fields[1])
        end = int(fields[2])
    except ValueError:
        return None

    # Optional fields from BED4/6/9
    name = fields[3] if len(fields) >= 4 else defaults.get("name", ".")
    score = fields[4] if len(fields) >= 5 else defaults.get("score", ".")
    strand = fields[5] if len(fields) >= 6 else defaults.get("strand", ".")

    # Thick coordinates (BED9)
    if len(fields) >= 8:
        thick_start = int(fields[6])
        thick_end = int(fields[7])
    else:
        thick_start = start
        thick_end = end

    # Color (BED9)
    if len(fields) >= 9:
        item_rgb = fields[8]
    else:
        item_rgb = defaults.get("color", ".")

    # Single block spanning entire interval
    block_count = 1
    block_sizes = str(end - start)
    block_starts = "0"

    return [
        chrom,
        str(start),
        str(end),
        name,
        score,
        strand,
        str(thick_start),
        str(thick_end),
        item_rgb,
        str(block_count),
        block_sizes,
        block_starts,
    ]


def merge_blocks(blocks):
    """
    Merge overlapping or adjacent blocks.

    Args:
        blocks: List of (start, end) tuples

    Returns:
        List of merged (start, end) tuples
    """
    if not blocks:
        return []

    # Sort by start position
    blocks = sorted(blocks)

    merged = [blocks[0]]
    for current_start, current_end in blocks[1:]:
        last_start, last_end = merged[-1]

        # Merge if overlapping or adjacent
        if current_start <= last_end:
            merged[-1] = (last_start, max(last_end, current_end))
        else:
            merged.append((current_start, current_end))

    return merged


def blocks_to_bed12(chrom, name, strand, blocks, score=".", color="0,0,0"):
    """
    Convert a set of blocks to BED12 format.

    Args:
        chrom: Chromosome name
        name: Feature name
        strand: Strand ('+', '-', '.')
        blocks: List of (start, end) tuples
        score: Score value
        color: RGB color string

    Returns:
        BED12 format string
    """
    if not blocks:
        return ""

    blocks = merge_blocks(blocks)

    # Overall coordinates
    start = min(b[0] for b in blocks)
    end = max(b[1] for b in blocks)

    # Thick coordinates
    thick_start = start
    thick_end = end

    block_count = len(blocks)

    # Calculate block sizes and relative starts
    block_sizes = [str(b[1] - b[0]) for b in blocks]
    block_starts = [str(b[0] - start) for b in blocks]

    block_sizes_str = ",".join(block_sizes) + ","
    block_starts_str = ",".join(block_starts) + ","

    return "\t".join(
        [
            chrom,
            str(start),
            str(end),
            name,
            score,
            strand,
            str(thick_start),
            str(thick_end),
            color,
            str(block_count),
            block_sizes_str,
            block_starts_str,
        ]
    )


def group_bed6_by_name(lines, max_distance=100000):
    """
    Group BED6 lines by name, chromosome, and strand.

    Args:
        lines: Iterator of BED6 lines
        max_distance: Maximum distance between features to group together

    Yields:
        BED12 format strings
    """
    # Group features by (chrom, strand, name)
    groups = defaultdict(list)

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        fields = line.split("\t")
        if len(fields) < 6:
            continue

        chrom = fields[0]
        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError:
            continue

        name = fields[3]
        score = fields[4]
        strand = fields[5]

        key = (chrom, strand, name)
        groups[key].append((start, end, score))

    # Convert each group to BED12
    for (chrom, strand, name), features in groups.items():
        # Sort by start position
        features.sort(key=lambda x: x[0])

        # Split into subgroups based on distance
        subgroups = []
        current_blocks = [(features[0][0], features[0][1])]
        current_score = features[0][2]

        for start, end, score in features[1:]:
            last_end = current_blocks[-1][1]

            if start - last_end <= max_distance:
                current_blocks.append((start, end))
            else:
                # Output current subgroup
                subgroups.append((current_blocks, current_score))
                current_blocks = [(start, end)]
                current_score = score

        # Add last subgroup
        subgroups.append((current_blocks, current_score))

        # Generate BED12 for each subgroup
        for i, (blocks, score) in enumerate(subgroups):
            if len(subgroups) > 1:
                subgroup_name = f"{name}_{i + 1}"
            else:
                subgroup_name = name

            yield blocks_to_bed12(chrom, subgroup_name, strand, blocks, score=score)


# CLI entry points


def bed6_to_bed12_cli():
    """Command-line tool: Convert BED6 to BED12 by grouping features with the same name."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert BED6 to BED12 by grouping features with the same name."
    )
    parser.add_argument(
        "--max-distance",
        type=int,
        default=100000,
        help="Maximum distance (bp) between features to group (default: 100000)",
    )
    args = parser.parse_args()

    for bed12 in group_bed6_by_name(sys.stdin, max_distance=args.max_distance):
        print(bed12)


def bedtobed12_cli():
    """Command-line tool: Convert BED3/4/6/9 to BED12 with single blocks."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert BED3/4/6/9 to BED12 format with single blocks."
    )
    parser.add_argument("--color", default=None, help="Default itemRgb color")
    parser.add_argument("--name", default=None, help="Default name")
    parser.add_argument("--score", default=None, help="Default score")
    parser.add_argument("--strand", default=None, help="Default strand")
    args = parser.parse_args()

    defaults = {
        "color": args.color,
        "name": args.name,
        "score": args.score,
        "strand": args.strand,
    }

    for line in sys.stdin:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        fields = line.split("\t")
        bed12_fields = bed_to_bed12(fields, defaults)

        if bed12_fields:
            print("\t".join(bed12_fields))
