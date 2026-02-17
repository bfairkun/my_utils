"""
SpliceAI utilities for variant analysis and splice site prediction.

This module provides tools for:
- Variant representation and validation
- Sequence fetching with variant application
- Coordinate mapping between reference and haplotype space
- Primer walk predictions with SpliceAI
"""

from typing import NamedTuple, Optional, List, Tuple, Dict
import pysam
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode


class Variant(NamedTuple):
    """Simple variant representation compatible with pysam.VariantRecord"""
    chrom: str
    pos: int  # 1-based position
    ref: str
    alt: str
    id: Optional[str] = None
    
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
            id=record.id
        )
    
    def to_tuple(self):
        """For backward compatibility with (pos, ref, alt) tuples"""
        return (self.pos, self.ref, self.alt)


class CoordinateMapper:
    """Maps coordinates between reference and haplotype space"""
    
    def __init__(self, region_start: int, variants: List[Variant]):
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
    
    def hap_to_ref(self, hap_idx: int) -> Optional[int]:
        """Convert haplotype sequence index (0-based) to reference position (1-based)
        Returns None if position falls within a deletion"""
        # This is more complex - need to account for insertions/deletions
        ref_offset = self.region_start
        hap_offset = 0
        
        for v in self.variants:
            variant_hap_start = hap_offset + (v.start - self.region_start)
            variant_hap_end = variant_hap_start + len(v.alt)
            
            if hap_idx < variant_hap_start:
                # Before this variant
                return ref_offset + hap_idx - hap_offset + 1  # Convert to 1-based
            elif hap_idx < variant_hap_end:
                # Within variant alt allele
                if v.is_insertion:
                    # Position is in inserted sequence, map to insertion point
                    return v.pos
                else:
                    # SNV or within alt of complex variant
                    offset_in_variant = hap_idx - variant_hap_start
                    return v.pos + offset_in_variant
            
            # Update offsets
            hap_offset += v.indel_length
            ref_offset = v.end
        
        # After all variants
        return ref_offset + (hap_idx - hap_offset) + 1


def fetch_sequence(fasta_path: str, chrom: str, start: int, end: int, 
                   strand: str = '+', variants: Optional[List[Variant]] = None) -> str:
    """
    Fetch sequence from FASTA file with optional SNV variants applied.
    Note: Only supports SNVs (single nucleotide variants). Indels will raise an error.
    
    Args:
        fasta_path: Path to FASTA file
        chrom: Chromosome name
        start: 0-based start position (inclusive)
        end: 0-based end position (exclusive)
        strand: '+' or '-'
        variants: Optional list of Variant objects (SNVs only)
    
    Returns:
        sequence: str, the (possibly mutated) sequence
    """
    with pysam.FastaFile(fasta_path) as fasta:
        sequence = fasta.fetch(chrom, start, end)
    
    if variants:
        # Filter variants to this region and validate
        variants_in_region = []
        for v in variants:
            if v.chrom != chrom:
                continue
            if v.start < start or v.end > end:
                continue
            
            # Only allow SNVs
            if not v.is_snv:
                raise ValueError(f"Only SNVs are supported. Got: {v}")
            
            # Validate REF
            v.validate_ref(fasta_path)
            variants_in_region.append(v)
        
        # Apply SNVs (no need to worry about offsets since length doesn't change)
        if variants_in_region:
            seq_list = list(sequence)
            for v in variants_in_region:
                idx = v.start - start
                seq_list[idx] = v.alt
            sequence = ''.join(seq_list)
    
    # Apply strand
    if strand == "-":
        sequence = str(Seq(sequence).reverse_complement())
    
    return sequence


def walk_sequence_with_Ns(sequence: str, walk_range: range, 
                          substitution_length: int = 20,
                          immutable_range: Optional[Tuple[int, int]] = None) -> List[Tuple[Tuple[int, int], str]]:
    """
    Generate sequences with N substitutions at specified positions.
    
    Args:
        sequence: Input sequence to mutate
        walk_range: Iterable of positions where substitutions should start
        substitution_length: Length of N substitution (default 20)
        immutable_range: Optional tuple (start, end) defining a region that cannot be substituted.
                        Immutable bases will be preserved while surrounding bases are masked.
    
    Returns:
        List of ((start_idx, end_idx), mutated_sequence) tuples
    """
    substituted_sequences = []
    
    for substitution_start_index in walk_range:
        substitution_end_index = substitution_start_index + substitution_length
        
        # Check if substitution is out of bounds
        if substitution_start_index < 0 or substitution_end_index > len(sequence):
            continue
        
        # Build the mutated sequence
        if immutable_range is not None:
            immutable_start, immutable_end = immutable_range
            
            # Check if there's overlap with immutable region
            overlap_start = max(substitution_start_index, immutable_start)
            overlap_end = min(substitution_end_index, immutable_end)
            
            if overlap_start < overlap_end:  # There is overlap
                # Mask before immutable region
                before_mask = 'N' * (overlap_start - substitution_start_index)
                # Keep immutable region
                preserved = sequence[overlap_start:overlap_end]
                # Mask after immutable region
                after_mask = 'N' * (substitution_end_index - overlap_end)
                
                mutated_sequence = (
                    sequence[:substitution_start_index] +
                    before_mask +
                    preserved +
                    after_mask +
                    sequence[substitution_end_index:]
                )
            else:  # No overlap
                mutated_sequence = (
                    sequence[:substitution_start_index] +
                    'N' * substitution_length +
                    sequence[substitution_end_index:]
                )
        else:
            mutated_sequence = (
                sequence[:substitution_start_index] +
                'N' * substitution_length +
                sequence[substitution_end_index:]
            )
        
        substituted_sequences.append(((substitution_start_index, substitution_end_index), mutated_sequence))
    
    return substituted_sequences


def load_spliceai_models():
    """Load the SpliceAI models"""
    paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
    return [load_model(resource_filename('spliceai', x)) for x in paths]


def predict_splice_sites(sequence: str, models: list, 
                         start_pos: Optional[int] = None,
                         end_pos: Optional[int] = None,
                         pad_size: int = 5000,
                         shift_to_intron: bool = True) -> pd.DataFrame:
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
    padded_sequence = 'N' * pad_size + sequence + 'N' * pad_size
    
    # One-hot encode the padded sequence
    x = one_hot_encode(padded_sequence)[None, :]
    
    # Predict using all models and average the results
    y = np.mean([model.predict(x) for model in models], axis=0)
    
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
    
    df = pd.DataFrame({
        'donor_prob': donor_prob,
        'acceptor_prob': acceptor_prob,
        'seq': list(seq_slice)
    }, index=idx)
    
    return df


def interval_notation_to_genomic_pos(interval: str) -> Tuple[str, int, int]:
    """
    Convert interval notation to genomic positions.
    
    Args:
        interval: String like 'chr2:1,840,760-1,887,609' (1-based, fully closed)
    
    Returns:
        Tuple of (chrom, start, end) where start is 0-based inclusive, end is 0-based exclusive
    """
    chrom, positions = interval.split(':')
    start, end = map(int, positions.replace(',', '').split('-'))
    return chrom, start - 1, end


def spliceO_predictions(
    fasta_path: str,
    region_interval: str,
    walk_interval: Tuple[int, int],
    walk_step: int,
    substitution_length: int,
    sites_of_interest: Dict[str, int],
    models: list,
    strand: str = "+",
    variants: Optional[List[Variant]] = None,
    immutable_range: Optional[Tuple[int, int]] = None
) -> pd.DataFrame:
    """
    Perform SpliceO primer walk predictions with optional SNV variants.
    Note: Only SNVs are supported.
    
    Args:
        fasta_path: Path to reference FASTA
        region_interval: String like 'chr2:1,840,760-1,887,609'
        walk_interval: Tuple of (start, end) genomic positions (1-based) for walk
        walk_step: Step size for walk
        substitution_length: Length of N-mer substitution
        sites_of_interest: Dict of {name: genomic_position} (1-based)
        models: Loaded SpliceAI models
        strand: '+' or '-'
        variants: Optional list of SNV Variant objects to apply
        immutable_range: Optional tuple of (start, end) genomic positions (1-based, fully closed)
                        that should not be masked with N's. Bases around this range will be masked.
    
    Returns:
        DataFrame with predictions for sites of interest across all masked positions
    """
    # Parse region and fetch sequence (with optional SNVs applied)
    chrom, region_start, region_end = interval_notation_to_genomic_pos(region_interval)
    seq = fetch_sequence(fasta_path, chrom, region_start, region_end, 
                        strand=strand, variants=variants)
    
    # Build genomic position array (maps each sequence index to genomic coordinate)
    seq_len = len(seq)
    if strand == "-":
        genomic_pos_array = list(range(region_end, region_end - seq_len, -1))
    else:
        genomic_pos_array = list(range(region_start + 1, region_start + 1 + seq_len))
    
    # Get predictions for the (possibly mutated) sequence
    predictions_wt = predict_splice_sites(seq, models)
    predictions_wt = predictions_wt.reset_index(drop=True)
    predictions_wt['GenomicPos'] = genomic_pos_array
    predictions_wt = predictions_wt.set_index('GenomicPos').sort_index()
    
    # Convert genomic walk_interval to sequence indices
    if strand == "+":
        walk_start_idx = walk_interval[0] - 1 - region_start
        walk_end_idx = walk_interval[1] - 1 - region_start
    else:
        walk_start_idx = region_end - walk_interval[1]
        walk_end_idx = region_end - walk_interval[0]
    
    # Convert genomic immutable_range to sequence indices if provided
    immutable_range_seq = None
    if immutable_range is not None:
        if strand == "+":
            immutable_start_idx = immutable_range[0] - 1 - region_start
            immutable_end_idx = immutable_range[1] - region_start  # Exclusive end
        else:
            # For minus strand, flip the coordinates
            immutable_start_idx = region_end - immutable_range[1]
            immutable_end_idx = region_end - immutable_range[0] + 1  # Exclusive end
        immutable_range_seq = (immutable_start_idx, immutable_end_idx)
    
    # Perform primer walks (mask with N's)
    walk_substitutions = walk_sequence_with_Ns(
        seq, range(walk_start_idx, walk_end_idx, walk_step), substitution_length,
        immutable_range=immutable_range_seq
    )
    
    # Make predictions for each masked sequence
    predictions = {}
    mask_sequences = {}
    aso_sequences = {}
    
    for (start_idx, end_idx), seq_mut in walk_substitutions:
        # Convert sequence indices to genomic coordinates for the mask
        if strand == "+":
            mask_start_genomic = region_start + start_idx + 1
            mask_end_genomic = region_start + end_idx
        else:
            mask_start_genomic = region_end - end_idx + 1
            mask_end_genomic = region_end - start_idx
        
        # Extract context sequence (±5 bp around mask with N's)
        context_start = max(0, start_idx - 5)
        context_end = min(len(seq), end_idx + 5)
        mask_context = seq_mut[context_start:context_end]
        
        # Get ASO sequence (reverse complement of masked region)
        masked_seq = seq[start_idx:end_idx]
        aso_seq = str(Seq(masked_seq).reverse_complement())
        
        # Make predictions for this masked sequence
        preds = predict_splice_sites(seq_mut, models)
        preds = preds.reset_index(drop=True)
        preds['GenomicPos'] = genomic_pos_array
        preds = preds.set_index('GenomicPos').sort_index()
        
        predictions[(mask_start_genomic, mask_end_genomic)] = preds
        mask_sequences[(mask_start_genomic, mask_end_genomic)] = mask_context
        aso_sequences[(mask_start_genomic, mask_end_genomic)] = aso_seq
    
    # Concatenate all predictions
    predictions = pd.concat(predictions.values(), keys=predictions.keys())
    predictions.index.names = ['MaskStart', 'MaskEnd', 'GenomicPos']
    
    # Select only sites of interest
    selected = predictions.loc[pd.IndexSlice[:, :, list(sites_of_interest.values())], :]
    selected_reset = selected.reset_index()
    
    # Add mask context and ASO sequences
    selected_reset['Mask_Context'] = selected_reset.apply(
        lambda row: mask_sequences.get((row['MaskStart'], row['MaskEnd']), ''), axis=1
    )
    selected_reset['ASO_Sequence'] = selected_reset.apply(
        lambda row: aso_sequences.get((row['MaskStart'], row['MaskEnd']), ''), axis=1
    )
    
    # Add wildtype control (unmasked sequence)
    wt_selected = predictions_wt.loc[list(sites_of_interest.values())].reset_index()
    wt_selected['Site'] = [k for k, v in sites_of_interest.items() 
                           if v in wt_selected['GenomicPos'].values]
    wt_selected['MaskStart'] = 'WT'
    wt_selected['MaskEnd'] = 'WT'
    wt_selected['Mask_Context'] = ''
    wt_selected['ASO_Sequence'] = ''
    
    # Combine masked and wildtype results
    final_df = pd.concat([
        selected_reset.assign(Site=selected_reset['GenomicPos'].map(
            {v: k for k, v in sites_of_interest.items()}
        )),
        wt_selected
    ], ignore_index=True)
    
    return final_df