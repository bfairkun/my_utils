"""
SpliceAI utilities for variant analysis and splice site prediction.

This module provides tools for:
- Variant representation and validation
- Sequence fetching with variant application
- Coordinate mapping between reference and haplotype space
- Primer walk predictions with SpliceAI
"""

from typing import NamedTuple

import numpy as np
import pandas as pd
import pysam
from Bio.Seq import Seq
from importlib.resources import files
from keras.models import load_model

# Local re-implementation of spliceai.utils.one_hot_encode.
# Upstream uses np.fromstring(seq, np.int8) which was removed in numpy>=2 —
# importing/calling it raises ValueError on modern envs. This version uses
# np.frombuffer instead, and also maps any non-ACGTN character to the
# all-zero row (upstream silently maps them to arbitrary rows via ord(c) % 5).
_ONE_HOT_MAP = np.asarray(
    [[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
    dtype=np.int8,
)


def one_hot_encode(seq: str) -> np.ndarray:
    """One-hot encode a DNA sequence. A=row 0, C=1, G=2, T=3; N and any
    other character → all-zero row. Output shape: (len(seq), 4), dtype int8.
    """
    buf = bytearray(seq.upper().encode("latin-1"))
    # Map ACGT to 1..4; everything else (including N and any non-ACGT char) to 0.
    table = bytearray(256)
    table[ord("A")] = 1
    table[ord("C")] = 2
    table[ord("G")] = 3
    table[ord("T")] = 4
    idx = np.frombuffer(buf.translate(table), dtype=np.uint8)
    return _ONE_HOT_MAP[idx]


class Variant(NamedTuple):
    """Simple variant representation compatible with pysam.VariantRecord"""

    chrom: str
    pos: int  # 1-based position
    ref: str
    alt: str
    id: str | None = None

    @property
    def start(self):
        """0-based start position (pysam convention)"""
        return self.pos - 1

    @property
    def end(self):
        """0-based end position (pysam convention)"""
        return self.pos - 1 + len(self.ref)

    @property
    def is_snv(self):
        return len(self.ref) == 1 and len(self.alt) == 1

    @property
    def is_insertion(self):
        return len(self.ref) < len(self.alt)

    @property
    def is_deletion(self):
        return len(self.ref) > len(self.alt)

    @property
    def indel_length(self):
        """Returns positive for insertion, negative for deletion, 0 for SNV"""
        return len(self.alt) - len(self.ref)

    def validate_ref(self, fasta_path: str) -> bool:
        """Check that REF matches reference genome"""
        with pysam.FastaFile(fasta_path) as fasta:
            actual_ref = fasta.fetch(self.chrom, self.start, self.end)
            if actual_ref.upper() != self.ref.upper():
                raise ValueError(
                    f"REF mismatch at {self.chrom}:{self.pos}: "
                    f"expected {self.ref}, found {actual_ref}"
                )
        return True

    @classmethod
    def from_pysam(cls, record: pysam.VariantRecord, alt_index: int = 0):
        """Convert pysam.VariantRecord to Variant (for compatibility)"""
        return cls(
            chrom=record.chrom,
            pos=record.pos,
            ref=record.ref,
            alt=record.alts[alt_index],
            id=record.id,
        )

    def to_tuple(self):
        """For backward compatibility with (pos, ref, alt) tuples"""
        return (self.pos, self.ref, self.alt)


class CoordinateMapper:
    """Maps coordinates between reference and haplotype space"""

    def __init__(self, region_start: int, variants: list[Variant]):
        """
        Args:
            region_start: 0-based start of region
            variants: List of variants (must be sorted by position)
        """
        self.region_start = region_start
        self.variants = sorted(variants, key=lambda v: v.pos)
        self._build_mapping()

    def _build_mapping(self):
        """Build coordinate mapping arrays"""
        # Calculate cumulative offset at each variant position
        self.variant_offsets = {}
        cumulative_offset = 0

        for v in self.variants:
            self.variant_offsets[v.pos] = cumulative_offset
            cumulative_offset += v.indel_length

    def ref_to_hap(self, ref_pos: int) -> int:
        """Convert reference position (1-based) to haplotype sequence index (0-based)"""
        # Find cumulative offset up to this position
        offset = 0
        for v in self.variants:
            if v.pos <= ref_pos:
                offset += v.indel_length
            else:
                break

        # Convert: 1-based ref → 0-based ref → apply offset → hap index
        return (ref_pos - 1 - self.region_start) + offset

    def hap_to_ref(self, hap_idx: int) -> int | None:
        """Convert haplotype sequence index (0-based) to reference position (1-based).

        Walks the haplotype left to right, tracking the haplotype cursor and the
        corresponding 0-based reference cursor across each variant. Bases inside an
        insertion's alt allele map to the variant's reference position.
        """
        # 0-based reference position aligned with `hap_cursor` in haplotype space.
        ref_cursor = self.region_start
        hap_cursor = 0

        for v in self.variants:
            # Reference segment before this variant maps 1:1 to haplotype space.
            pre_len = v.start - ref_cursor
            if hap_idx < hap_cursor + pre_len:
                return ref_cursor + (hap_idx - hap_cursor) + 1  # -> 1-based

            # The variant's alt allele occupies `len(v.alt)` haplotype bases.
            alt_hap_start = hap_cursor + pre_len
            alt_hap_end = alt_hap_start + len(v.alt)
            if hap_idx < alt_hap_end:
                if v.is_insertion:
                    return v.pos  # inserted base -> insertion point
                # SNV / first base(s) of a complex alt -> aligned ref base(s)
                return v.pos + (hap_idx - alt_hap_start)

            # Advance past the variant.
            hap_cursor = alt_hap_end
            ref_cursor = v.end  # 0-based end of the variant's ref span

        # Trailing reference segment after all variants.
        return ref_cursor + (hap_idx - hap_cursor) + 1


def fetch_sequence(
    fasta_path: str,
    chrom: str,
    start: int,
    end: int,
    strand: str = "+",
    variants: list[Variant] | None = None,
) -> str:
    """
    Fetch sequence from FASTA file with optional variants applied.

    Supports SNVs and indels. Variants use a left-anchored, VCF-style representation
    (ref and alt share their leading anchor base(s)): an insertion has len(alt) >
    len(ref), a deletion has len(alt) < len(ref), an SNV has len(ref) == len(alt) == 1.
    Variants must lie fully within [start, end). The same representation is assumed by
    CoordinateMapper in spliceO_predictions when remapping genomic <-> haplotype indices.

    Args:
        fasta_path: Path to FASTA file
        chrom: Chromosome name
        start: 0-based start position (inclusive)
        end: 0-based end position (exclusive)
        strand: '+' or '-'
        variants: Optional list of Variant objects (SNVs and/or indels)

    Returns:
        sequence: str, the (possibly mutated) sequence
    """
    with pysam.FastaFile(fasta_path) as fasta:
        sequence = fasta.fetch(chrom, start, end)

    if variants:
        # Filter variants fully contained in this region (plus-strand coords)
        variants_in_region = []
        for v in variants:
            if v.chrom != chrom:
                continue
            if v.start < start or v.end > end:
                continue
            # REF must match the reference genome
            v.validate_ref(fasta_path)
            variants_in_region.append(v)

        # Apply variants (SNVs and indels) on the plus strand. Apply right-to-left
        # (descending start) so earlier offsets stay valid as lengths change.
        if variants_in_region:
            seq_list = list(sequence)
            for v in sorted(variants_in_region, key=lambda v: v.start, reverse=True):
                i = v.start - start
                seq_list[i : i + len(v.ref)] = list(v.alt)
            sequence = "".join(seq_list)

    # Apply strand
    if strand == "-":
        sequence = str(Seq(sequence).reverse_complement())

    return sequence


def walk_sequence_with_Ns(
    sequence: str,
    walk_range: range,
    substitution_length: int = 20,
    immutable_ranges: list[tuple[int, int]] | None = None,
) -> list[tuple[tuple[int, int], str]]:
    """
    Generate sequences with N substitutions at specified positions.

    Args:
        sequence: Input sequence to mutate
        walk_range: Iterable of positions where substitutions should start
        substitution_length: Length of N substitution (default 20)
        immutable_ranges: Optional list of (start, end) tuples defining regions that
                         cannot be substituted. Each tuple is a half-open interval
                         [start, end). Bases within any of these ranges will be
                         preserved while surrounding bases are masked.

    Returns:
        List of ((start_idx, end_idx), mutated_sequence) tuples
    """
    # Pre-build set of protected positions for O(1) lookup
    protected_positions: set[int] = set()
    if immutable_ranges is not None:
        for immutable_start, immutable_end in immutable_ranges:
            protected_positions.update(range(immutable_start, immutable_end))

    substituted_sequences = []

    for substitution_start_index in walk_range:
        substitution_end_index = substitution_start_index + substitution_length

        # Check if substitution is out of bounds
        if substitution_start_index < 0 or substitution_end_index > len(sequence):
            continue

        # Build masked sequence: replace non-protected positions with N
        seq_list = list(sequence)
        for i in range(substitution_start_index, substitution_end_index):
            if i not in protected_positions:
                seq_list[i] = "N"
        mutated_sequence = "".join(seq_list)

        substituted_sequences.append(
            ((substitution_start_index, substitution_end_index), mutated_sequence)
        )

    return substituted_sequences


def load_spliceai_models():
    """Load the SpliceAI models"""
    paths = (f"models/spliceai{x}.h5" for x in range(1, 6))
    return [load_model(str(files("spliceai") / x)) for x in paths]


def predict_splice_sites(
    sequence: str,
    models: list,
    start_pos: int | None = None,
    end_pos: int | None = None,
    pad_size: int = 5000,
    shift_to_intron: bool = True,
) -> pd.DataFrame:
    """
    Predict splice sites for a sequence using SpliceAI.

    Args:
        sequence: DNA sequence
        models: List of loaded SpliceAI models
        start_pos: Start position in sequence (0-based, inclusive)
        end_pos: End position in sequence (0-based, exclusive)
        pad_size: Padding size for context (default 5000)
        shift_to_intron: Shift predictions by 1bp to align with intron boundaries

    Returns:
        DataFrame with columns: donor_prob, acceptor_prob, seq
        Index is 0-based position in input sequence
    """
    # Pad the sequence with 'N' on both sides for context
    padded_sequence = "N" * pad_size + sequence + "N" * pad_size

    # One-hot encode the padded sequence
    x = one_hot_encode(padded_sequence)[None, :]

    # Predict using all models and average the results
    y = np.mean([model.predict(x, verbose=0) for model in models], axis=0)

    # y has shape (1, L, 3) where L is the length of predictions
    # SpliceAI predictions are shorter than input by 2*pad_size
    # So if input is N bases, predictions are for N - 2*pad_size bases
    # But we padded with pad_size on each side, so predictions cover the original sequence

    # Default to entire unpadded sequence if range not specified
    if start_pos is None:
        start_pos = 0
    if end_pos is None:
        end_pos = len(sequence)

    # Extract predictions for the requested range
    # The predictions y[0] correspond to the original sequence (unpadded)
    # because we added pad_size padding, and SpliceAI drops pad_size from each end
    acceptor_prob = y[0, start_pos:end_pos, 1]
    donor_prob = y[0, start_pos:end_pos, 2]

    if shift_to_intron:
        # Shift donor right by 1 (to first base of intron)
        # Shift acceptor left by 1 (to last base of intron)
        donor_prob = np.roll(donor_prob, 1)
        donor_prob[0] = 0  # Clear wraparound
        acceptor_prob = np.roll(acceptor_prob, -1)
        acceptor_prob[-1] = 0  # Clear wraparound

    # Build DataFrame with 0-based indices relative to the unpadded sequence
    idx = np.arange(start_pos, end_pos)
    seq_slice = sequence[start_pos:end_pos]

    df = pd.DataFrame(
        {
            "donor_prob": donor_prob,
            "acceptor_prob": acceptor_prob,
            "seq": list(seq_slice),
        },
        index=idx,
    )

    return df


def interval_notation_to_genomic_pos(interval: str) -> tuple[str, int, int]:
    """
    Convert interval notation to genomic positions.

    Args:
        interval: String like 'chr2:1,840,760-1,887,609' (1-based, fully closed)

    Returns:
        Tuple of (chrom, start, end) where start is 0-based inclusive, end is 0-based exclusive
    """
    chrom, positions = interval.split(":")
    start, end = map(int, positions.replace(",", "").split("-"))
    return chrom, start - 1, end


def spliceO_predictions(
    fasta_path: str,
    region_interval: str,
    walk_interval: tuple[int, int],
    walk_step: int,
    substitution_length: int,
    sites_of_interest: dict[str, int],
    models: list | None = None,
    strand: str = "+",
    variants: list[Variant] | None = None,
    immutable_ranges: list[tuple[int, int]] | None = None,
    predictor_fn=None,
) -> pd.DataFrame:
    """
    Perform SpliceO primer walk predictions with optional variants (SNVs and indels).
    Variants are applied to the haplotype sequence and all genomic <-> sequence-index
    conversions are routed through CoordinateMapper, so sites_of_interest, walk_interval,
    immutable_ranges, and the reported MaskStart/MaskEnd stay in genomic coordinates even
    when an indel changes the sequence length.

    Args:
        fasta_path: Path to reference FASTA
        region_interval: String like 'chr2:1,840,760-1,887,609'
        walk_interval: Tuple of (start, end) genomic positions (1-based) for walk
        walk_step: Step size for walk
        substitution_length: Length of N-mer substitution
        sites_of_interest: Dict of {name: genomic_position} (1-based)
        models: Loaded SpliceAI models (optional when predictor_fn is provided)
        strand: '+' or '-'
        variants: Optional list of Variant objects (SNVs and/or indels) to apply
        immutable_ranges: Optional list of (start, end) genomic position tuples (1-based,
                         fully closed) that should not be masked with N's. Bases within
                         any of these ranges will be preserved while surrounding bases
                         are masked.
        predictor_fn: Optional callable(seq_str) -> DataFrame(donor_prob, acceptor_prob, seq).
                      When provided, used instead of SpliceAI models.

    Returns:
        DataFrame with predictions for sites of interest across all masked positions
    """
    if predictor_fn is None and models is None:
        raise ValueError("Either models or predictor_fn must be provided")
    # Parse region and fetch the (possibly variant-applied) haplotype sequence
    chrom, region_start, region_end = interval_notation_to_genomic_pos(region_interval)
    seq = fetch_sequence(
        fasta_path, chrom, region_start, region_end, strand=strand, variants=variants
    )
    seq_len = len(seq)

    # Coordinate mapping between 1-based genomic positions and haplotype sequence
    # indices, accounting for any indels in `variants`. CoordinateMapper works in
    # plus-strand haplotype space; for the minus strand the SpliceAI sequence is the
    # reverse complement, so seq index s <-> plus-hap index h via h = seq_len - 1 - s.
    # With no indels this reduces to the original 1:1 index<->genomic mapping.
    mapper = CoordinateMapper(region_start, list(variants) if variants else [])

    def gpos_to_seqidx(gpos: int) -> int:
        h = mapper.ref_to_hap(gpos)
        return h if strand == "+" else (seq_len - 1 - h)

    def seqidx_to_gpos(sidx: int) -> int:
        h = sidx if strand == "+" else (seq_len - 1 - sidx)
        return mapper.hap_to_ref(h)

    if predictor_fn is not None:
        predictor = predictor_fn
    else:
        predictor = lambda s: predict_splice_sites(s, models)  # noqa: E731

    # Map sites of interest (1-based genomic) -> haplotype sequence indices
    site_seqidx = {}
    for name, gpos in sites_of_interest.items():
        sidx = gpos_to_seqidx(gpos)
        if not 0 <= sidx < seq_len:
            raise ValueError(
                f"Site '{name}' at {chrom}:{gpos} maps outside the region "
                f"(or falls within a deletion); seq index {sidx}."
            )
        site_seqidx[name] = sidx

    # Wildtype (unmasked) predictions on the (possibly variant-applied) sequence
    predictions_wt = predictor(seq).reset_index(drop=True)

    # Convert genomic walk_interval -> sequence index range (orientation-agnostic)
    walk_start_idx, walk_end_idx = sorted(
        (gpos_to_seqidx(walk_interval[0]), gpos_to_seqidx(walk_interval[1]))
    )

    # Convert genomic immutable_ranges -> half-open sequence index ranges
    immutable_ranges_seq = None
    if immutable_ranges is not None:
        immutable_ranges_seq = []
        for g0, g1 in immutable_ranges:
            lo, hi = sorted((gpos_to_seqidx(g0), gpos_to_seqidx(g1)))
            immutable_ranges_seq.append((lo, hi + 1))

    # Perform primer walks (mask with N's)
    walk_substitutions = walk_sequence_with_Ns(
        seq,
        range(walk_start_idx, walk_end_idx, walk_step),
        substitution_length,
        immutable_ranges=immutable_ranges_seq,
    )

    def site_rows(preds_df, mask_start, mask_end, mask_ctx, aso_seq):
        rows = []
        for name, sidx in site_seqidx.items():
            r = preds_df.iloc[sidx]
            rows.append(
                {
                    "MaskStart": mask_start,
                    "MaskEnd": mask_end,
                    "GenomicPos": seqidx_to_gpos(sidx),
                    "donor_prob": r["donor_prob"],
                    "acceptor_prob": r["acceptor_prob"],
                    "seq": r["seq"],
                    "Site": name,
                    "Mask_Context": mask_ctx,
                    "ASO_Sequence": aso_seq,
                }
            )
        return rows

    records = []
    for (start_idx, end_idx), seq_mut in walk_substitutions:
        # Genomic bounds of the masked window (min/max over covered positions, so the
        # values are correct even when the mask spans an indel).
        covered = [seqidx_to_gpos(s) for s in range(start_idx, end_idx)]
        mask_start_genomic = min(covered)
        mask_end_genomic = max(covered)

        # Context sequence (±5 bp around mask, with N's)
        context_start = max(0, start_idx - 5)
        context_end = min(seq_len, end_idx + 5)
        mask_context = seq_mut[context_start:context_end]

        # ASO sequence = reverse complement of the masked region
        masked_seq = seq[start_idx:end_idx]
        aso_seq = str(Seq(masked_seq).reverse_complement())

        preds = predictor(seq_mut).reset_index(drop=True)
        records.extend(
            site_rows(preds, mask_start_genomic, mask_end_genomic, mask_context, aso_seq)
        )

    # Wildtype control rows (MaskStart/MaskEnd = "WT")
    records.extend(site_rows(predictions_wt, "WT", "WT", "", ""))

    return pd.DataFrame.from_records(records)


def summarize_spliceO_walk(
    df: pd.DataFrame,
    summary_func: callable,
    score_col: str = "score",
) -> pd.DataFrame:
    """
    Reshape spliceO_predictions output to one row per ASO with an aggregated score.

    Args:
        df: Output DataFrame from spliceO_predictions (or a saved/loaded version).
        summary_func: Callable that takes a sub-DataFrame (all Site rows for a
            single ASO, or the WT rows) and returns a scalar score.
            The sub-DataFrame has the same columns as df.
            Examples:
              lambda d: d.loc[d["Site"]=="cryptic_donor", "donor_prob"].iloc[0]
              lambda d: (d.loc[d["Site"]=="E18_donor", "donor_prob"].iloc[0]
                         + d.loc[d["Site"]=="E18_acceptor", "acceptor_prob"].iloc[0])
              lambda d: (d.loc[d["Site"]=="cryptic_donor", "donor_prob"].iloc[0]
                         / d.loc[d["Site"]=="E18_donor", "donor_prob"].iloc[0])
        score_col: Name for the score column in the output. Default "score".

    Returns:
        DataFrame with one row per ASO (unique MaskStart/MaskEnd), columns:
          - MaskStart (int)
          - MaskEnd (int)
          - ASO_Sequence
          - Mask_Context
          - {score_col}     : summary_func applied to this ASO's rows
          - {score_col}_wt  : summary_func applied to the WT rows (constant)
    """
    wt_df = df[df["MaskStart"] == "WT"]
    aso_df = df[df["MaskStart"] != "WT"]

    wt_score = summary_func(wt_df)

    group_keys = ["MaskStart", "MaskEnd", "ASO_Sequence", "Mask_Context"]
    rows = []
    for (mask_start, mask_end, aso_seq, mask_ctx), group in aso_df.groupby(
        group_keys, sort=False
    ):
        rows.append(
            {
                "MaskStart": int(mask_start),
                "MaskEnd": int(mask_end),
                "ASO_Sequence": aso_seq,
                "Mask_Context": mask_ctx,
                score_col: summary_func(group),
            }
        )

    result = pd.DataFrame(rows)
    result[f"{score_col}_wt"] = wt_score
    return result
