"""
Junction Usage Index (JUI) utilities.

Quantifies how exclusively a splice junction is used relative to competing
junctions and intron retention (IR) at the same splice sites.

JUI is on [0, 1]; higher = more exclusive junction usage.

Two input modes:
  Mode A: BAM file (parses CIGAR for junctions + pileup for IR)
  Mode B: junction table (SJ.out.tab or regtools BED12) + optional bigwig for IR

Coordinate convention:
  donor_pos    – 0-based, last exonic base on the 5' side of the intron
  acceptor_pos – 0-based, first exonic base on the 3' side of the intron
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass


# ── Coordinate / format helpers ───────────────────────────────────────────────

def _parse_region(region: str) -> tuple[str, int, int]:
    """Parse 'chr7:98100000-98250000' → (chrom, start, end) 0-based half-open."""
    chrom, coords = region.split(":", 1)
    start_s, end_s = coords.split("-", 1)
    return chrom, int(start_s), int(end_s)


def _in_region(chrom: str, pos: int, region: tuple[str, int, int] | None) -> bool:
    if region is None:
        return True
    r_chrom, r_start, r_end = region
    return chrom == r_chrom and r_start <= pos < r_end


# ── BAM parsing (Mode A) ──────────────────────────────────────────────────────

BAM_CREF_SKIP = 3  # CIGAR 'N'
BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8

_CONSUMES_REF = {BAM_CMATCH, BAM_CDEL, BAM_CREF_SKIP, BAM_CEQUAL, BAM_CDIFF}


def _strand_from_read(read, library_strand: str) -> str:
    """Return '+', '-', or '.' for a read given library strandedness."""
    if library_strand == "unstranded":
        # Fall back to XS tag if present (STAR / HISAT2 sets this for spliced reads)
        xs = read.get_tag("XS") if read.has_tag("XS") else None
        return xs if xs in ("+", "-") else "."
    if library_strand == "forward":
        if read.is_read1:
            return "-" if read.is_reverse else "+"
        return "+" if read.is_reverse else "-"
    # reverse
    if read.is_read1:
        return "+" if read.is_reverse else "-"
    return "-" if read.is_reverse else "+"


def parse_junctions_from_bam(
    bam_path: str,
    region: str | None = None,
    library_strand: str = "unstranded",
    min_intron_size: int = 70,
    max_intron_size: int = 500_000,
    min_anchor: int = 8,
    include_secondary: bool = False,
) -> dict[tuple[str, int, int, str], int]:
    """
    Walk CIGAR tuples and extract splice junctions from a BAM.

    Filters match regtools junction extract defaults:
      min_intron_size: minimum N-op length to call a junction (default 70, avoids indels)
      max_intron_size: maximum N-op length (default 500,000)
      min_anchor: minimum exonic bases on each side of the junction (default 8).
                  Applied per-junction: a junction passes if *any* supporting read
                  has >= min_anchor bases on both sides.

    Returns dict[(chrom, donor_pos, acceptor_pos, strand)] → read count.
    donor_pos    = 0-based last exonic base on 5' side
    acceptor_pos = 0-based first exonic base on 3' side
    """
    import pysam

    counts: dict[tuple[str, int, int, str], int] = defaultdict(int)
    # Track per-junction max anchor lengths across reads for post-hoc filtering
    max_left_anchor: dict[tuple[str, int, int, str], int] = defaultdict(int)
    max_right_anchor: dict[tuple[str, int, int, str], int] = defaultdict(int)

    fetch_kwargs: dict = {}
    parsed_region = None
    if region is not None:
        parsed_region = _parse_region(region)
        fetch_kwargs = {
            "contig": parsed_region[0],
            "start": parsed_region[1],
            "stop": parsed_region[2],
        }

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(**fetch_kwargs):
            if read.is_unmapped or read.is_supplementary:
                continue
            if read.is_secondary and not include_secondary:
                continue
            if read.cigartuples is None:
                continue

            strand = _strand_from_read(read, library_strand)
            chrom = read.reference_name
            ref_pos = read.reference_start  # 0-based

            # Collect exon block lengths and candidate junctions in one pass.
            # exon_lengths[i] is the left anchor for junctions[i] and the right
            # anchor for junctions[i-1].
            exon_lengths: list[int] = []
            junc_candidates: list[tuple[str, int, int, str] | None] = []
            current_exon = 0

            for op, length in read.cigartuples:
                if op == BAM_CREF_SKIP:
                    exon_lengths.append(current_exon)
                    if min_intron_size <= length <= max_intron_size:
                        donor_pos = ref_pos - 1
                        acceptor_pos = ref_pos + length
                        key = (chrom, donor_pos, acceptor_pos, strand)
                        if parsed_region is None or _in_region(chrom, donor_pos, parsed_region):
                            junc_candidates.append(key)
                        else:
                            junc_candidates.append(None)
                    else:
                        junc_candidates.append(None)  # filtered by size
                    current_exon = 0
                    ref_pos += length
                elif op in _CONSUMES_REF:
                    current_exon += length
                    ref_pos += length
            exon_lengths.append(current_exon)  # final exon block

            # Record counts and update max anchors per junction
            for i, key in enumerate(junc_candidates):
                if key is None:
                    continue
                left_anchor = exon_lengths[i]
                right_anchor = exon_lengths[i + 1]
                counts[key] += 1
                if left_anchor > max_left_anchor[key]:
                    max_left_anchor[key] = left_anchor
                if right_anchor > max_right_anchor[key]:
                    max_right_anchor[key] = right_anchor

    # Apply per-junction anchor filter
    return {
        key: count
        for key, count in counts.items()
        if max_left_anchor[key] >= min_anchor and max_right_anchor[key] >= min_anchor
    }


def get_ir_counts_from_bam(
    bam_path: str,
    splice_sites: list[tuple[str, int, str]],  # (chrom, pos, strand)
    ir_overhang: int = 3,
    region: str | None = None,
) -> dict[tuple[str, int, str], float]:
    """
    For each splice site position, count reads that are truly aligned (not
    spliced across) at pos ± ir_overhang using pysam pileup.

    Returns dict[(chrom, pos, strand)] → IR read count.
    """
    import pysam

    ir_counts: dict[tuple[str, int, str], float] = {}

    fetch_region: tuple[str, int, int] | None = None
    if region is not None:
        fetch_region = _parse_region(region)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, pos, strand in splice_sites:
            # Donor: query pos + ir_overhang (inside intron)
            # Acceptor: query pos - ir_overhang (inside intron)
            # The caller provides pos as either donor_pos or acceptor_pos;
            # they handle the offset logic — we just query `pos` here.
            if fetch_region and not _in_region(chrom, pos, fetch_region):
                ir_counts[(chrom, pos, strand)] = 0.0
                continue

            count = 0
            try:
                for pileup_col in bam.pileup(
                    chrom, pos, pos + 1, truncate=True, stepper="nofilter"
                ):
                    if pileup_col.reference_pos == pos:
                        for pileup_read in pileup_col.pileups:
                            if not pileup_read.is_refskip and not pileup_read.is_del:
                                count += 1
            except ValueError:
                pass
            ir_counts[(chrom, pos, strand)] = float(count)

    return ir_counts


def build_ir_dicts_from_bam(
    bam_path: str,
    junction_counts: dict[tuple[str, int, int, str], int],
    ir_overhang: int = 3,
    region: str | None = None,
) -> tuple[
    dict[tuple[str, int, str], float],
    dict[tuple[str, int, str], float],
]:
    """
    For all donor and acceptor positions from junction_counts, fetch IR
    coverage from BAM at the intronic overhang positions.

    Returns (ir_at_donors, ir_at_acceptors):
      ir_at_donors   dict[(chrom, donor_pos, strand)] → IR count at donor+ir_overhang
      ir_at_acceptors dict[(chrom, acceptor_pos, strand)] → IR count at acceptor-ir_overhang
    """
    donor_sites: set[tuple[str, int, str]] = set()
    acceptor_sites: set[tuple[str, int, str]] = set()
    for chrom, donor, acceptor, strand in junction_counts:
        donor_sites.add((chrom, donor + ir_overhang, strand))
        acceptor_sites.add((chrom, acceptor - ir_overhang - 1, strand))

    # map query positions back to the original splice site positions
    donor_query_to_site = {
        (chrom, donor + ir_overhang, strand): (chrom, donor, strand)
        for chrom, donor, _acc, strand in junction_counts
    }
    acceptor_query_to_site = {
        (chrom, acceptor - ir_overhang - 1, strand): (chrom, acceptor, strand)
        for chrom, _don, acceptor, strand in junction_counts
    }

    raw_donor = get_ir_counts_from_bam(
        bam_path, list(donor_sites), ir_overhang=0, region=region
    )
    raw_acceptor = get_ir_counts_from_bam(
        bam_path, list(acceptor_sites), ir_overhang=0, region=region
    )

    ir_at_donors: dict[tuple[str, int, str], float] = {}
    for qkey, site_key in donor_query_to_site.items():
        ir_at_donors[site_key] = raw_donor.get(qkey, 0.0)

    ir_at_acceptors: dict[tuple[str, int, str], float] = {}
    for qkey, site_key in acceptor_query_to_site.items():
        ir_at_acceptors[site_key] = raw_acceptor.get(qkey, 0.0)

    return ir_at_donors, ir_at_acceptors


# ── Junction table parsing (Mode B) ──────────────────────────────────────────

_STRAND_MAP_SJTAB = {"0": ".", "1": "+", "2": "-"}


def detect_junction_format(path: str) -> str:
    """Auto-detect 'sj_out_tab' or 'bed12' by inspecting the first data line."""
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= 12:
                # BED12 has strand in col 5 as +/-/.
                if fields[5] in ("+", "-", "."):
                    return "bed12"
            # SJ.out.tab col 3 is 0, 1, or 2
            if len(fields) >= 4 and fields[3] in ("0", "1", "2"):
                return "sj_out_tab"
            # Fall back
            return "bed12"
    return "bed12"


def parse_sj_out_tab(
    path: str,
    region: tuple[str, int, int] | None = None,
    use_multi: bool = False,
) -> dict[tuple[str, int, int, str], int]:
    """
    Parse STAR SJ.out.tab file.

    Coord conversion (STAR uses 1-based inclusive intron coordinates):
      donor_pos    = intron_start - 2  (0-based last exon base, 5' side)
      acceptor_pos = intron_end        (1-based intron end = 0-based first exon base, 3' side)
    """
    counts: dict[tuple[str, int, int, str], int] = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom = fields[0]
            try:
                intron_start = int(fields[1])
                intron_end = int(fields[2])
                unique_reads = int(fields[6])
                multi_reads = int(fields[7])
            except ValueError:
                continue

            strand = _STRAND_MAP_SJTAB.get(fields[3], ".")
            donor_pos = intron_start - 2
            acceptor_pos = intron_end  # 1-based last intron base = 0-based first exon base

            if region and not _in_region(chrom, donor_pos, region):
                continue

            count = unique_reads + (multi_reads if use_multi else 0)
            key = (chrom, donor_pos, acceptor_pos, strand)
            counts[key] = counts.get(key, 0) + count

    return counts


def parse_bed12_junctions(
    path: str,
    region: tuple[str, int, int] | None = None,
) -> dict[tuple[str, int, int, str], int]:
    """
    Parse regtools BED12 junction file.

    Coord conversion (2-block BED12 with exon anchors):
      donor_pos    = start + blockSizes[0] - 1
      acceptor_pos = start + blockStarts[1]
    Read count from score column.
    """
    counts: dict[tuple[str, int, int, str], int] = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                continue
            chrom = fields[0]
            try:
                start = int(fields[1])
                score = int(fields[4])
                strand = fields[5]
                block_sizes = [int(x) for x in fields[10].rstrip(",").split(",")]
                block_starts = [int(x) for x in fields[11].rstrip(",").split(",")]
            except (ValueError, IndexError):
                continue

            if len(block_sizes) < 2 or len(block_starts) < 2:
                continue

            donor_pos = start + block_sizes[0] - 1
            acceptor_pos = start + block_starts[1]

            if region and not _in_region(chrom, donor_pos, region):
                continue

            key = (chrom, donor_pos, acceptor_pos, strand)
            counts[key] = counts.get(key, 0) + score

    return counts


# ── BigWig IR (Mode B optional) ───────────────────────────────────────────────

def get_ir_from_bigwig(
    bw_path: str,
    junction_counts: dict[tuple[str, int, int, str], int],
    ir_offset: int = 3,
) -> tuple[
    dict[tuple[str, int, str], float],
    dict[tuple[str, int, str], float],
]:
    """
    Query bigwig at donor+ir_offset and acceptor-ir_offset-1 for each junction.

    Returns (ir_at_donors, ir_at_acceptors).
    """
    try:
        import pyBigWig
    except ImportError as e:
        msg = "pyBigWig is required for bigwig IR. Install with: pip install pyBigWig"
        raise ImportError(msg) from e

    ir_at_donors: dict[tuple[str, int, str], float] = {}
    ir_at_acceptors: dict[tuple[str, int, str], float] = {}

    bw = pyBigWig.open(bw_path)
    try:
        for chrom, donor, acceptor, strand in junction_counts:
            d_query = donor + ir_offset
            a_query = acceptor - ir_offset - 1
            if a_query < 0:
                a_query = 0

            d_key = (chrom, donor, strand)
            a_key = (chrom, acceptor, strand)

            if d_key not in ir_at_donors:
                try:
                    val = bw.stats(chrom, d_query, d_query + 1, type="mean")[0]
                    ir_at_donors[d_key] = float(val) if val is not None else 0.0
                except RuntimeError:
                    ir_at_donors[d_key] = 0.0

            if a_key not in ir_at_acceptors:
                try:
                    val = bw.stats(chrom, a_query, a_query + 1, type="mean")[0]
                    ir_at_acceptors[a_key] = float(val) if val is not None else 0.0
                except RuntimeError:
                    ir_at_acceptors[a_key] = 0.0
    finally:
        bw.close()

    return ir_at_donors, ir_at_acceptors


# ── Core JUI computation ──────────────────────────────────────────────────────

def compute_jui(
    junction_counts: dict[tuple[str, int, int, str], int],
    ir_at_donors: dict[tuple[str, int, str], float],
    ir_at_acceptors: dict[tuple[str, int, str], float],
    min_count: int = 1,
    decompose: bool = False,
    ir_available: bool = True,
    prior: float = 0.5,
) -> list[dict]:
    """
    Compute JUI for every junction in junction_counts.

    Returns a list of dicts (one per junction), suitable for TSV output.
    """
    # Pre-compute donor and acceptor competition sums
    donor_total: dict[tuple[str, int, str], int] = defaultdict(int)
    acceptor_total: dict[tuple[str, int, str], int] = defaultdict(int)
    for (chrom, donor, acceptor, strand), k in junction_counts.items():
        donor_total[(chrom, donor, strand)] += k
        acceptor_total[(chrom, acceptor, strand)] += k

    rows = []
    for (chrom, donor, acceptor, strand), k in sorted(junction_counts.items()):
        d_key = (chrom, donor, strand)
        a_key = (chrom, acceptor, strand)

        C_D = donor_total[d_key] - k
        C_A = acceptor_total[a_key] - k
        R_D = ir_at_donors.get(d_key, 0.0)
        R_A = ir_at_acceptors.get(a_key, 0.0)

        n_D = k + C_D + R_D
        n_A = k + C_A + R_A

        row: dict = {
            "chrom": chrom,
            "donor_pos": donor,
            "acceptor_pos": acceptor,
            "strand": strand,
            "k": k,
            "C_D": C_D,
            "C_A": C_A,
            "R_D": R_D,
            "R_A": R_A,
            "n_D": n_D,
            "n_A": n_A,
            "JUI_D": "NA",
            "JUI_A": "NA",
            "JUI": "NA",
            "JUI_shrunk": "NA",
            "bottleneck_side": "NA",
            "ir_available": ir_available,
        }

        if n_D >= min_count and n_A >= min_count:
            JUI_D = k / n_D
            JUI_A = k / n_A
            JUI = min(JUI_D, JUI_A)
            JUI_shrunk = min((prior + k) / (2 * prior + n_D), (prior + k) / (2 * prior + n_A))
            bottleneck = "D" if JUI_D <= JUI_A else "A"
            row.update(
                JUI_D=round(JUI_D, 6),
                JUI_A=round(JUI_A, 6),
                JUI=round(JUI, 6),
                JUI_shrunk=round(JUI_shrunk, 6),
                bottleneck_side=bottleneck,
            )

            if decompose:
                splice_D_denom = k + C_D
                splice_A_denom = k + C_A
                JUI_D_splice = k / splice_D_denom if splice_D_denom > 0 else "NA"
                JUI_A_splice = k / splice_A_denom if splice_A_denom > 0 else "NA"
                JUI_D_retention = (k + C_D) / n_D if n_D > 0 else "NA"
                JUI_A_retention = (k + C_A) / n_A if n_A > 0 else "NA"

                if JUI_D_retention != "NA" and JUI_A_retention != "NA":
                    retention_index = 1 - min(JUI_D_retention, JUI_A_retention)
                else:
                    retention_index = "NA"

                if JUI_D_splice != "NA" and JUI_A_splice != "NA":
                    competition_index = 1 - min(JUI_D_splice, JUI_A_splice)
                else:
                    competition_index = "NA"

                row.update(
                    JUI_D_splice=_fmt(JUI_D_splice),
                    JUI_A_splice=_fmt(JUI_A_splice),
                    JUI_D_retention=_fmt(JUI_D_retention),
                    JUI_A_retention=_fmt(JUI_A_retention),
                    RetentionIndex=_fmt(retention_index),
                    CompetitionIndex=_fmt(competition_index),
                )

        rows.append(row)
    return rows


def _fmt(v) -> str:
    if v == "NA":
        return "NA"
    return str(round(float(v), 6))


# ── Output formatters ──────────────────────────────────────────────────────────

TSV_COLUMNS = [
    "chrom", "donor_pos", "acceptor_pos", "strand",
    "k", "C_D", "C_A", "R_D", "R_A", "n_D", "n_A",
    "JUI_D", "JUI_A", "JUI", "JUI_shrunk", "bottleneck_side",
    "input_mode", "ir_available",
]

DECOMPOSE_COLUMNS = [
    "JUI_D_splice", "JUI_A_splice",
    "JUI_D_retention", "JUI_A_retention",
    "RetentionIndex", "CompetitionIndex",
]


def write_tsv(rows: list[dict], fh, input_mode: str, decompose: bool = False) -> None:
    cols = list(TSV_COLUMNS)
    if decompose:
        cols += DECOMPOSE_COLUMNS
    print("\t".join(cols), file=fh)
    for row in rows:
        row["input_mode"] = input_mode
        print("\t".join(str(row.get(c, "NA")) for c in cols), file=fh)


def rows_to_bed12(rows: list[dict], bed_score: str = "raw") -> list[str]:
    """
    Convert JUI rows to BED12 lines for IGV junction arc visualization.
    Each junction becomes a 2-block BED12 (10bp anchors at each end).
    Score = round(JUI * 1000), clamped [0, 1000]; 0 if JUI is NA.
    Includes 'track graphType=junctions' header so IGV renders arcs.
    """
    lines = ["track graphType=junctions"]
    anchor = 10
    for row in rows:
        chrom = row["chrom"]
        donor = int(row["donor_pos"])
        acceptor = int(row["acceptor_pos"])
        strand = row["strand"]

        jui_val = row.get("JUI_shrunk" if bed_score == "shrunk" else "JUI", "NA")
        if jui_val == "NA":
            score = 0
        else:
            score = min(50, max(0, round(float(jui_val) * 50)))

        name = f"{chrom}:{donor}-{acceptor}"
        start = donor
        end = acceptor + 1
        thick_start = start
        thick_end = end
        rgb = "0,0,0"
        block_count = 2
        block_sizes = f"{anchor},{anchor}"
        block_starts = f"0,{acceptor - donor - anchor + 1}"

        lines.append(
            "\t".join(
                str(x)
                for x in [
                    chrom, start, end, name, score, strand,
                    thick_start, thick_end, rgb,
                    block_count, block_sizes, block_starts,
                ]
            )
        )
    return lines


# ── CLI ────────────────────────────────────────────────────────────────────────

def compute_jui_cli() -> None:
    parser = argparse.ArgumentParser(
        prog="compute-jui",
        description=(
            "Compute Junction Usage Index (JUI) from RNA-seq data.\n\n"
            "Two input modes (one required):\n"
            "  Mode A -- BAM file:          --bam <file> [--region] [--strand] [--ir-overhang] [--multi-mappers]\n"
            "  Mode B -- junction table:    --junc-table <file> [--bigwig <file>] [--junc-format] [--ir-offset] [--multi-mappers]\n\n"
            "JUI_shrunk uses a Beta(prior, prior) Bayesian prior, shrinking toward 0.5 "
            "(maximum uncertainty). Larger --prior = stronger shrinkage toward 0.5 at low counts. "
            "At large n, JUI_shrunk ≈ JUI regardless of prior."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--bam",
        metavar="BAM",
        help="[Mode A] Sorted, indexed BAM file. Junctions extracted from CIGAR; IR from pileup.",
    )
    input_group.add_argument(
        "--junc-table",
        metavar="FILE",
        help="[Mode B] Junction table: STAR SJ.out.tab or regtools BED12 (.junc). Format auto-detected.",
    )

    modeB = parser.add_argument_group("Mode B options (junction table input)")
    modeB.add_argument(
        "--bigwig",
        metavar="BW",
        help="Coverage bigwig for intron-retention counts (optional; if omitted IR=0 for all junctions).",
    )
    modeB.add_argument(
        "--junc-format",
        choices=["sj_out_tab", "bed12"],
        default=None,
        help="Junction table format. Default: auto-detect from file contents.",
    )
    modeB.add_argument(
        "--ir-offset",
        type=int,
        default=3,
        metavar="BP",
        help="Intronic bp offset from splice site for bigwig IR query (default: 3).",
    )

    shared = parser.add_argument_group("Shared options")
    shared.add_argument(
        "--multi-mappers",
        action="store_true",
        help=(
            "Include multi-mapping reads. "
            "Mode A: include secondary BAM alignments (non-primary hits). "
            "Mode B (SJ.out.tab): sum unique + multi-mapping read columns."
        ),
    )
    shared.add_argument(
        "--region",
        metavar="CHR:START-END",
        default=None,
        help="Restrict to a genomic region, e.g. chr16:185000-233000.",
    )
    shared.add_argument(
        "--strand",
        choices=["forward", "reverse", "unstranded"],
        default="unstranded",
        help=(
            "Library strandedness. Mode A: used to assign strand to reads. "
            "Mode B: informational only (strand already in file). Default: unstranded."
        ),
    )
    shared.add_argument(
        "--ir-overhang",
        type=int,
        default=3,
        metavar="BP",
        help="[Mode A] Intronic bp offset from splice site for pileup IR query (default: 3).",
    )
    shared.add_argument(
        "--min-intron",
        type=int,
        default=70,
        metavar="BP",
        help="[Mode A] Minimum intron size to call a junction from CIGAR (default: 70, matching regtools).",
    )
    shared.add_argument(
        "--max-intron",
        type=int,
        default=500_000,
        metavar="BP",
        help="[Mode A] Maximum intron size (default: 500,000, matching regtools).",
    )
    shared.add_argument(
        "--min-anchor",
        type=int,
        default=8,
        metavar="BP",
        help=(
            "[Mode A] Minimum exonic anchor on each side of junction (default: 8, matching regtools). "
            "Applied per-junction: kept if any supporting read meets the threshold on both sides."
        ),
    )
    shared.add_argument(
        "--min-count",
        type=int,
        default=1,
        metavar="N",
        help="Minimum reads (n_D and n_A) required to compute JUI; otherwise NA (default: 1).",
    )
    shared.add_argument(
        "--prior",
        type=float,
        default=0.5,
        metavar="A",
        help=(
            "Beta(prior, prior) prior pseudocount for JUI_shrunk. "
            "Shrinks toward 0.5; larger = stronger shrinkage at low counts (default: 0.5 = Jeffreys prior)."
        ),
    )
    shared.add_argument("--decompose", action="store_true", help="Add decomposed JUI columns (RetentionIndex, CompetitionIndex, etc.).")

    output = parser.add_argument_group("Output options")
    output.add_argument("--output-tsv", default="-", metavar="FILE", help="TSV output path (default: stdout).")
    output.add_argument("--output-bed", default=None, metavar="FILE", help="BED12 junction arc file for IGV (score = round(JUI × 50)).")
    output.add_argument(
        "--bed-score",
        choices=["raw", "shrunk"],
        default="raw",
        help="Which JUI value to encode as BED arc score: raw JUI or JUI_shrunk (default: raw).",
    )

    args = parser.parse_args()

    region_tuple: tuple[str, int, int] | None = None
    if args.region:
        region_tuple = _parse_region(args.region)

    # ── Mode A: BAM ────────────────────────────────────────────────────────────
    if args.bam:
        junction_counts = parse_junctions_from_bam(
            args.bam,
            region=args.region,
            library_strand=args.strand,
            min_intron_size=args.min_intron,
            max_intron_size=args.max_intron,
            min_anchor=args.min_anchor,
            include_secondary=args.multi_mappers,
        )
        ir_at_donors, ir_at_acceptors = build_ir_dicts_from_bam(
            args.bam,
            junction_counts,
            ir_overhang=args.ir_overhang,
            region=args.region,
        )
        input_mode = "bam"
        ir_available = True

    # ── Mode B: junction table ─────────────────────────────────────────────────
    else:
        fmt = args.junc_format or detect_junction_format(args.junc_table)
        if fmt == "sj_out_tab":
            junction_counts = parse_sj_out_tab(
                args.junc_table,
                region=region_tuple,
                use_multi=args.multi_mappers,
            )
        else:
            junction_counts = parse_bed12_junctions(args.junc_table, region=region_tuple)

        if args.bigwig:
            ir_at_donors, ir_at_acceptors = get_ir_from_bigwig(
                args.bigwig, junction_counts, ir_offset=args.ir_offset
            )
            ir_available = True
        else:
            ir_at_donors = {}
            ir_at_acceptors = {}
            ir_available = False

        input_mode = fmt

    rows = compute_jui(
        junction_counts,
        ir_at_donors,
        ir_at_acceptors,
        min_count=args.min_count,
        decompose=args.decompose,
        ir_available=ir_available,
        prior=args.prior,
    )

    # TSV output
    if args.output_tsv == "-":
        write_tsv(rows, sys.stdout, input_mode=input_mode, decompose=args.decompose)
    else:
        with open(args.output_tsv, "w") as fh:
            write_tsv(rows, fh, input_mode=input_mode, decompose=args.decompose)

    # BED12 output
    if args.output_bed:
        bed_lines = rows_to_bed12(rows, bed_score=args.bed_score)
        with open(args.output_bed, "w") as fh:
            for line in bed_lines:
                print(line, file=fh)
