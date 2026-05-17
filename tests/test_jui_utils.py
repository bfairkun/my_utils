"""Tests for jui_utils.py"""

from __future__ import annotations

import os
import pytest

from my_utils.jui_utils import (
    compute_jui,
    detect_junction_format,
    parse_bed12_junctions,
    parse_junctions_from_bam,
    parse_sj_out_tab,
    rows_to_bed12,
    _parse_region,
)

# ── Test BAMs (skip if not present) ───────────────────────────────────────────

DRUG_BAM = (
    "/project/yangili1/bjf79/20260310_diversesm_dr/code/rna_seq/Alignments/"
    "BPN15477_24h_noStim_25000nM_rep1/Aligned.sortedByCoord.out.bam"
)
DMSO_BAM = (
    "/project/yangili1/bjf79/20260310_diversesm_dr/code/rna_seq/Alignments/"
    "DMSO_24h_noStim_0nM_rep1/Aligned.sortedByCoord.out.bam"
)
DRUG_JUNC = (
    "/project/yangili1/bjf79/20260310_diversesm_dr/code/rna_seq/SplicingAnalysis/"
    "juncfiles/BPN15477_24h_noStim_25000nM_rep1.junc"
)
LUC7L_REGION = "chr7:98100000-98250000"

BAMS_AVAILABLE = os.path.exists(DRUG_BAM) and os.path.exists(DRUG_BAM + ".bai")


# ── Unit tests: coordinate parsing ────────────────────────────────────────────

def test_parse_region():
    chrom, start, end = _parse_region("chr7:98100000-98250000")
    assert chrom == "chr7"
    assert start == 98100000
    assert end == 98250000


def test_sj_out_tab_coord_conversion(tmp_path):
    """
    STAR SJ.out.tab uses 1-based inclusive intron coords.
    intron_start=100, intron_end=200 → donor=98, acceptor=199
    """
    tab = tmp_path / "test.sj"
    # cols: chrom intron_start intron_end strand motif annotated unique multi max_overhang
    tab.write_text("chr1\t100\t200\t1\t2\t1\t15\t3\t29\n")
    counts = parse_sj_out_tab(str(tab))
    assert len(counts) == 1
    (chrom, donor, acceptor, strand), k = next(iter(counts.items()))
    assert chrom == "chr1"
    assert donor == 98    # intron_start(100) - 2
    assert acceptor == 200  # intron_end(200) in 1-based = 0-based first exon base
    assert strand == "+"
    assert k == 15  # unique reads


def test_sj_out_tab_multi_mappers(tmp_path):
    tab = tmp_path / "test.sj"
    tab.write_text("chr1\t100\t200\t2\t2\t1\t10\t5\t29\n")
    counts_no_multi = parse_sj_out_tab(str(tab), use_multi=False)
    counts_with_multi = parse_sj_out_tab(str(tab), use_multi=True)
    k_no = next(iter(counts_no_multi.values()))
    k_with = next(iter(counts_with_multi.values()))
    assert k_no == 10
    assert k_with == 15


def test_sj_out_tab_strand_codes(tmp_path):
    tab = tmp_path / "test.sj"
    # strand codes: 0=undef, 1=+, 2=-
    tab.write_text(
        "chr1\t100\t200\t0\t0\t0\t5\t0\t20\n"
        "chr1\t300\t400\t1\t0\t0\t5\t0\t20\n"
        "chr1\t500\t600\t2\t0\t0\t5\t0\t20\n"
    )
    counts = parse_sj_out_tab(str(tab))
    strands = {strand for (_, _, _, strand) in counts}
    assert "." in strands
    assert "+" in strands
    assert "-" in strands


def test_bed12_coord_conversion(tmp_path):
    """
    regtools BED12: start=83393, blockSizes=152,144, blockStarts=0,533
    donor_pos  = 83393 + 152 - 1 = 83544
    acceptor_pos = 83393 + 533   = 83926
    """
    bed = tmp_path / "test.junc"
    bed.write_text(
        "chr1\t83393\t83927\tname\t10\t+\t83393\t83927\t0,0,0\t2\t152,144\t0,533\n"
    )
    counts = parse_bed12_junctions(str(bed))
    assert len(counts) == 1
    (chrom, donor, acceptor, strand), k = next(iter(counts.items()))
    assert chrom == "chr1"
    assert donor == 83544
    assert acceptor == 83926
    assert strand == "+"
    assert k == 10


def test_detect_format_bed12(tmp_path):
    bed = tmp_path / "junc.bed"
    bed.write_text("chr1\t100\t200\tname\t5\t+\t100\t200\t0\t2\t10,10\t0,90\n")
    assert detect_junction_format(str(bed)) == "bed12"


def test_detect_format_sj_out_tab(tmp_path):
    tab = tmp_path / "SJ.out.tab"
    tab.write_text("chr1\t100\t200\t1\t2\t1\t15\t3\t29\n")
    assert detect_junction_format(str(tab)) == "sj_out_tab"


# ── Unit tests: compute_jui ────────────────────────────────────────────────────

def test_single_junction_no_competition_no_ir():
    """Single junction, no competitors, no IR → JUI should be 1.0."""
    junctions = {("chr1", 100, 200, "+"): 20}
    ir_d = {}
    ir_a = {}
    rows = compute_jui(junctions, ir_d, ir_a, min_count=5)
    assert len(rows) == 1
    row = rows[0]
    assert row["JUI"] == 1.0
    assert row["JUI_D"] == 1.0
    assert row["JUI_A"] == 1.0
    assert row["bottleneck_side"] == "D"  # tie goes to D


def test_min_count_filters():
    """Junction below min_count → JUI = NA."""
    junctions = {("chr1", 100, 200, "+"): 3}
    rows = compute_jui(junctions, {}, {}, min_count=5)
    assert rows[0]["JUI"] == "NA"


def test_competing_junctions():
    """Two junctions sharing a donor: each should have JUI_D = 0.5."""
    junctions = {
        ("chr1", 100, 200, "+"): 10,
        ("chr1", 100, 300, "+"): 10,
    }
    rows = compute_jui(junctions, {}, {}, min_count=5)
    by_acc = {r["acceptor_pos"]: r for r in rows}
    assert float(by_acc[200]["JUI_D"]) == pytest.approx(0.5)
    assert float(by_acc[300]["JUI_D"]) == pytest.approx(0.5)


def test_ir_reduces_jui():
    """IR at donor should reduce JUI_D."""
    junctions = {("chr1", 100, 200, "+"): 10}
    ir_d = {("chr1", 100, "+"): 10.0}  # equal IR as junction reads
    ir_a = {}
    rows = compute_jui(junctions, ir_d, ir_a, min_count=5)
    row = rows[0]
    assert float(row["JUI_D"]) == pytest.approx(0.5)
    assert float(row["JUI_A"]) == pytest.approx(1.0)
    assert row["bottleneck_side"] == "D"


def test_bayesian_shrinkage():
    """JUI_shrunk should equal (0.5 + k) / (1 + n) when n_D == n_A."""
    k = 10
    junctions = {("chr1", 100, 200, "+"): k}
    rows = compute_jui(junctions, {}, {}, min_count=5)
    row = rows[0]
    n = k  # no competition, no IR
    expected = (0.5 + k) / (1 + n)
    assert float(row["JUI_shrunk"]) == pytest.approx(expected)


def test_decompose_columns():
    """decompose=True should add extra columns."""
    junctions = {("chr1", 100, 200, "+"): 10}
    rows = compute_jui(junctions, {}, {}, min_count=5, decompose=True)
    row = rows[0]
    assert "RetentionIndex" in row
    assert "CompetitionIndex" in row
    assert "JUI_D_splice" in row


def test_decompose_no_competition_no_ir():
    """Single junction with no competition/IR: RetentionIndex=0, CompetitionIndex=0."""
    junctions = {("chr1", 100, 200, "+"): 10}
    rows = compute_jui(junctions, {}, {}, min_count=5, decompose=True)
    row = rows[0]
    assert row["RetentionIndex"] == "0.0"
    assert row["CompetitionIndex"] == "0.0"


# ── BED12 output format ────────────────────────────────────────────────────────

def test_bed12_output_columns():
    junctions = {("chr7", 1000, 2000, "+"): 20}
    rows = compute_jui(junctions, {}, {}, min_count=5)
    bed_lines = rows_to_bed12(rows)
    # First line is the track header
    assert bed_lines[0] == "track graphType=junctions"
    fields = bed_lines[1].split("\t")
    assert len(fields) == 12, f"Expected 12 BED12 columns, got {len(fields)}"


def test_bed12_score_range():
    """BED score must be in [0, 1000]."""
    junctions = {
        ("chr7", 1000, 2000, "+"): 50,
        ("chr7", 3000, 4000, "+"): 2,  # below min_count → JUI NA → score 0
    }
    rows = compute_jui(junctions, {}, {}, min_count=5)
    bed_lines = rows_to_bed12(rows)
    for line in bed_lines[1:]:  # skip header
        score = int(line.split("\t")[4])
        assert 0 <= score <= 50


def test_bed12_block_structure():
    """Block sizes should be 10,10 and starts should span the junction."""
    junctions = {("chr7", 1000, 2000, "+"): 20}
    rows = compute_jui(junctions, {}, {}, min_count=5)
    bed_lines = rows_to_bed12(rows)
    fields = bed_lines[1].split("\t")  # skip header
    assert fields[9] == "2"       # blockCount
    assert fields[10] == "10,10"  # blockSizes
    # blockStart[1] = acceptor - donor - 9 = 2000 - 1000 - 9 = 991
    starts = fields[11].split(",")
    assert starts[0] == "0"
    assert int(starts[1]) == 991


# ── BAM integration tests ──────────────────────────────────────────────────────

@pytest.mark.skipif(not BAMS_AVAILABLE, reason="Test BAMs not available")
def test_bam_junction_parsing_returns_results():
    counts = parse_junctions_from_bam(DRUG_BAM, region=LUC7L_REGION)
    assert len(counts) > 0


@pytest.mark.skipif(not BAMS_AVAILABLE, reason="Test BAMs not available")
def test_bam_jui_values_in_range():
    counts = parse_junctions_from_bam(DRUG_BAM, region=LUC7L_REGION)
    rows = compute_jui(counts, {}, {}, min_count=5)
    numeric_rows = [r for r in rows if r["JUI"] != "NA"]
    assert len(numeric_rows) > 0
    for row in numeric_rows:
        jui = float(row["JUI"])
        assert 0.0 <= jui <= 1.0, f"JUI out of range: {jui}"


@pytest.mark.skipif(not BAMS_AVAILABLE, reason="Test BAMs not available")
def test_bed12_junc_file_parsing():
    counts = parse_bed12_junctions(DRUG_JUNC, region=_parse_region(LUC7L_REGION))
    assert len(counts) > 0
