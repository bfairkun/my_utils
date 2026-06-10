"""
Coordinate-correctness tests for indel-aware spliceO_predictions / fetch_sequence.

These do NOT require SpliceAI or a GPU: spliceO_predictions is driven by a
deterministic mock `predictor_fn` whose `donor_prob`/`acceptor_prob` encode the
sequence index, so we can assert exactly which haplotype index (and base) each
site of interest reads. A tiny synthetic FASTA gives full control of the bases.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from my_utils.spliceai_utils import (  # noqa: E402
    CoordinateMapper,
    Variant,
    fetch_sequence,
    spliceO_predictions,
)

REF = "ACGT" * 25  # 100 bp synthetic reference; base at 0-based i is "ACGT"[i % 4]


@pytest.fixture(scope="module")
def fasta(tmp_path_factory):
    """Write the synthetic reference as chrT and index it."""
    d = tmp_path_factory.mktemp("ref")
    fa = d / "chrT.fa"
    fa.write_text(">chrT\n" + REF + "\n")
    pysam.faidx(str(fa))
    return str(fa)


def mock_predictor(seq: str) -> pd.DataFrame:
    """donor_prob == acceptor_prob == sequence index; seq column = the base.

    Lets a test assert both the index a site mapped to (via the prob value) and
    the base actually present there (via the seq value).
    """
    n = len(seq)
    return pd.DataFrame(
        {
            "donor_prob": np.arange(n, dtype=float),
            "acceptor_prob": np.arange(n, dtype=float),
            "seq": list(seq),
        },
        index=range(n),
    )


def revcomp(s: str) -> str:
    return str(Seq(s).reverse_complement())


# --------------------------------------------------------------------------- #
# CoordinateMapper round-trips
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "variants",
    [
        [],
        [Variant("chrT", 50, "C", "A")],  # SNV
        [Variant("chrT", 50, "C", "CC")],  # insertion (+1)
        [Variant("chrT", 50, "CG", "C")],  # deletion (-1)
    ],
)
def test_coordinatemapper_roundtrip_outside_variants(variants):
    """ref_to_hap then hap_to_ref is identity for positions away from variants."""
    m = CoordinateMapper(region_start=0, variants=variants)
    for gpos in (10, 30, 80, 95):  # all < or > the variant locus, never inside it
        assert m.hap_to_ref(m.ref_to_hap(gpos)) == gpos


def test_coordinatemapper_insertion_shifts_downstream():
    """A +1 insertion at pos 50 shifts plus-hap indices of pos > 50 by +1."""
    m = CoordinateMapper(region_start=0, variants=[Variant("chrT", 50, "C", "CC")])
    assert m.ref_to_hap(40) == 39  # before insertion: unchanged
    assert m.ref_to_hap(70) == 70  # after insertion: 69 + 1


# --------------------------------------------------------------------------- #
# fetch_sequence: SNV / insertion / deletion on + and - strand
# --------------------------------------------------------------------------- #
def test_fetch_sequence_snv_plus(fasta):
    seq = fetch_sequence(fasta, "chrT", 0, 100, "+", [Variant("chrT", 50, "C", "A")])
    assert seq == REF[:49] + "A" + REF[50:]
    assert len(seq) == 100


def test_fetch_sequence_insertion_plus(fasta):
    seq = fetch_sequence(fasta, "chrT", 0, 100, "+", [Variant("chrT", 50, "C", "CC")])
    assert seq == REF[:49] + "CC" + REF[50:]
    assert len(seq) == 101


def test_fetch_sequence_deletion_plus(fasta):
    seq = fetch_sequence(fasta, "chrT", 0, 100, "+", [Variant("chrT", 50, "CG", "C")])
    assert seq == REF[:49] + "C" + REF[51:]
    assert len(seq) == 99


def test_fetch_sequence_insertion_minus(fasta):
    seq = fetch_sequence(fasta, "chrT", 0, 100, "-", [Variant("chrT", 50, "C", "CC")])
    assert seq == revcomp(REF[:49] + "CC" + REF[50:])
    assert len(seq) == 101


def test_fetch_sequence_validate_ref_mismatch(fasta):
    """Wrong REF base must raise."""
    with pytest.raises(ValueError):
        fetch_sequence(fasta, "chrT", 0, 100, "+", [Variant("chrT", 50, "A", "T")])


# --------------------------------------------------------------------------- #
# spliceO_predictions: no-variant baseline on both strands
# --------------------------------------------------------------------------- #
def _wt_value(df, site, col="donor_prob"):
    row = df[(df["MaskStart"] == "WT") & (df["Site"] == site)]
    assert len(row) == 1
    return row.iloc[0][col]


def _walk_kwargs(fasta, strand, sites, variants=None):
    return dict(
        fasta_path=fasta,
        region_interval="chrT:1-100",
        walk_interval=(40, 60),
        walk_step=5,
        substitution_length=6,
        sites_of_interest=sites,
        predictor_fn=mock_predictor,
        strand=strand,
        variants=variants,
    )


def test_spliceO_no_variant_plus_indices(fasta):
    df = spliceO_predictions(**_walk_kwargs(fasta, "+", {"a": 30, "b": 70}))
    # plus, no indel: genomic g -> seq index g-1, mock donor_prob == index
    assert _wt_value(df, "a") == 29
    assert _wt_value(df, "b") == 69
    # 'seq' base must equal the reference base at that genomic position
    arow = df[(df["MaskStart"] == "WT") & (df["Site"] == "a")].iloc[0]
    assert arow["seq"] == REF[29]
    assert arow["GenomicPos"] == 30


def test_spliceO_no_variant_minus_indices(fasta):
    df = spliceO_predictions(**_walk_kwargs(fasta, "-", {"a": 30, "b": 70}))
    # minus, no indel: seq index = 100 - g
    assert _wt_value(df, "a") == 70
    assert _wt_value(df, "b") == 30
    arow = df[(df["MaskStart"] == "WT") & (df["Site"] == "a")].iloc[0]
    assert arow["seq"] == revcomp(REF)[70]
    assert arow["GenomicPos"] == 30


# --------------------------------------------------------------------------- #
# spliceO_predictions: the decisive dupG-style case (minus strand insertion)
# --------------------------------------------------------------------------- #
def test_spliceO_insertion_minus_donor_shift(fasta):
    """Mirrors the CLN3 dupG: insertion within the exon shifts the donor (low
    genomic, minus strand) seq index by +1 while the acceptor (high genomic) is
    unchanged — matching the prior notebook's donor9_idx_pat = donor9_idx + 1."""
    sites = {"acceptor": 70, "donor": 30}  # minus strand: acceptor high, donor low
    wt = spliceO_predictions(**_walk_kwargs(fasta, "-", sites))
    ins = spliceO_predictions(
        **_walk_kwargs(fasta, "-", sites, variants=[Variant("chrT", 50, "C", "CC")])
    )
    # acceptor (genomic 70 > insertion locus 50) — seq index unchanged
    assert _wt_value(wt, "acceptor") == _wt_value(ins, "acceptor") == 30
    # donor (genomic 30 < 50) — seq index shifts +1 under the +1 insertion
    assert _wt_value(wt, "donor") == 70
    assert _wt_value(ins, "donor") == 71


def test_spliceO_snv_changes_site_base(fasta):
    """An SNV at a monitored site changes the base read there, not the index."""
    sites = {"s": 50}
    wt = spliceO_predictions(**_walk_kwargs(fasta, "+", sites))
    snv = spliceO_predictions(
        **_walk_kwargs(fasta, "+", sites, variants=[Variant("chrT", 50, "C", "A")])
    )
    assert _wt_value(wt, "s", "donor_prob") == _wt_value(snv, "s", "donor_prob") == 49
    assert wt[(wt["MaskStart"] == "WT") & (wt["Site"] == "s")].iloc[0]["seq"] == "C"
    assert snv[(snv["MaskStart"] == "WT") & (snv["Site"] == "s")].iloc[0]["seq"] == "A"


# --------------------------------------------------------------------------- #
# spliceO_predictions: mask reporting, ASO sequence, columns
# --------------------------------------------------------------------------- #
def test_spliceO_mask_genomic_and_aso(fasta):
    df = spliceO_predictions(**_walk_kwargs(fasta, "+", {"a": 50}))
    masked = df[df["MaskStart"] != "WT"]
    assert len(masked) > 0
    # Required columns present (consumed by summarize_spliceO_walk / score_fn)
    for col in ["MaskStart", "MaskEnd", "donor_prob", "acceptor_prob",
                "Site", "Mask_Context", "ASO_Sequence"]:
        assert col in df.columns
    # MaskEnd - MaskStart + 1 == substitution_length (genomic, plus strand)
    r = masked.iloc[0]
    assert int(r["MaskEnd"]) - int(r["MaskStart"]) + 1 == 6
    # ASO_Sequence == reverse complement of the reference window it masked
    start0 = int(r["MaskStart"]) - 1
    assert r["ASO_Sequence"] == revcomp(REF[start0:start0 + 6])


def test_spliceO_immutable_ranges_preserve_bases(fasta):
    """Bases inside immutable_ranges are not masked to N."""
    df = spliceO_predictions(
        **_walk_kwargs(fasta, "+", {"a": 50}),
    )  # baseline
    df_prot = spliceO_predictions(
        fasta_path=fasta,
        region_interval="chrT:1-100",
        walk_interval=(40, 60),
        walk_step=1,
        substitution_length=6,
        sites_of_interest={"a": 50},
        predictor_fn=mock_predictor,
        strand="+",
        immutable_ranges=[(50, 51)],  # protect genomic 50-51
    )
    # Some mask context should retain a non-N base from the protected window
    ctxs = df_prot[df_prot["MaskStart"] != "WT"]["Mask_Context"].tolist()
    assert any(set(c) - {"N"} for c in ctxs)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
