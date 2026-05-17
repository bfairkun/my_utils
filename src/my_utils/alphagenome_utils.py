"""
AlphaGenome utilities for splice site prediction.

AlphaGenome (Google DeepMind, v0.5.1) predicts per-position splice site
probability/usage across many tissues via an HTTP API.

Supported sequence lengths: 16,384 / 131,072 / 524,288 / 1,048,576 bp
(sequence is padded/centered to the nearest supported length)

Install (in py_general — already available):
    pip install alphagenome

Usage:
    import os
    from my_utils import alphagenome_utils
    client = alphagenome_utils.create_alphagenome_client(
        os.environ["ALPHAGENOME_API_KEY"]
    )
    predictor = alphagenome_utils.make_alphagenome_predictor(client)
    # Pass to spliceO_predictions(predictor_fn=predictor)
"""

from __future__ import annotations

import warnings
from typing import Callable

import numpy as np
import pandas as pd

from alphagenome.data import ontology
from alphagenome.data.ontology import OntologyTerm
from alphagenome.models.dna_client import (
    Organism,
    OutputType,
    create as _ag_create,
)

# AlphaGenome supported sequence lengths (must match exactly)
_SUPPORTED_LENGTHS = [16_384, 131_072, 524_288, 1_048_576]


def create_alphagenome_client(api_key: str):
    """
    Initialize and return an AlphaGenome API client.

    Args:
        api_key: AlphaGenome API key.

    Returns:
        AlphaGenome client object.
    """
    return _ag_create(api_key=api_key)


def list_alphagenome_tissues(client, output_type: str = "splice_site_usage") -> pd.DataFrame:
    """
    Return a DataFrame of available tissue tracks and their ontology CURIEs.

    Args:
        client: AlphaGenome client from create_alphagenome_client().
        output_type: "splice_site_usage" or "splice_sites"

    Returns:
        DataFrame with columns: name, ontology_curie (if available)
    """
    otype = OutputType.SPLICE_SITE_USAGE if output_type == "splice_site_usage" else OutputType.SPLICE_SITES

    # Use a short dummy sequence to probe available tracks
    dummy_seq = "A" * 16_384
    try:
        output = client.predict_sequence(
            dummy_seq,
            organism=Organism.HOMO_SAPIENS,
            requested_outputs=[otype],
            ontology_terms=[],
        )
        data = output.splice_site_usage if output_type == "splice_site_usage" else output.splice_sites
        if hasattr(data, "metadata") and data.metadata is not None:
            rows = []
            for i, meta in enumerate(data.metadata):
                rows.append({
                    "index": i,
                    "name": getattr(meta, "name", None) or getattr(meta, "tissue", None) or str(i),
                    "ontology_curie": getattr(meta, "ontology_term", None),
                })
            return pd.DataFrame(rows)
        else:
            n = data.values.shape[-1] if hasattr(data, "values") else "unknown"
            return pd.DataFrame({"index": range(n), "name": [f"track_{i}" for i in range(n)]})
    except Exception as e:
        raise RuntimeError(f"Could not probe AlphaGenome tissue list: {e}") from e


def _pad_sequence_centered(sequence: str, target_length: int) -> tuple[str, int]:
    """
    Pad sequence to target_length by centering, adding N's on each side.

    Returns:
        padded_sequence, left_pad_size
    """
    seq_len = len(sequence)
    if seq_len > target_length:
        raise ValueError(
            f"Sequence length {seq_len} exceeds target_length {target_length}"
        )
    pad_total = target_length - seq_len
    left_pad = pad_total // 2
    right_pad = pad_total - left_pad
    return "N" * left_pad + sequence + "N" * right_pad, left_pad


def _choose_target_length(seq_len: int) -> int:
    """Choose smallest AlphaGenome-supported length that fits sequence."""
    for length in _SUPPORTED_LENGTHS:
        if seq_len <= length:
            return length
    raise ValueError(
        f"Sequence length {seq_len} exceeds maximum AlphaGenome input "
        f"({_SUPPORTED_LENGTHS[-1]} bp)"
    )


def predict_splice_sites_alphagenome(
    sequence: str,
    client,
    tissue: str | None = None,
    output_type: str = "usage",
    target_length: int | None = None,
    shift_to_intron: bool = True,
    sense_strand: str | None = "+",
    aggregation: str = "median",
) -> pd.DataFrame:
    """
    Predict splice sites using AlphaGenome.

    Args:
        sequence: DNA sequence (ACGTN). Will be padded to AlphaGenome-supported length.
            The caller supplies the sequence in the intended orientation; AlphaGenome
            tracks for the matching strand (`sense_strand`) are what apply to it.
        client: AlphaGenome client from create_alphagenome_client().
        tissue: Ontology CURIE (e.g. "UBERON:0000955" for brain) or None →
            average all tissues on the selected strand. When None a UserWarning
            is emitted, since multi-tissue averaging hides the aggregation choice.
        output_type: "usage" (SPLICE_SITE_USAGE) or "probability" (SPLICE_SITES).
        target_length: AlphaGenome input length to use (auto-selected if None).
        shift_to_intron: Shift donor +1, acceptor -1 to align to intron boundaries.
        sense_strand: "+" (default), "-", or None. Restricts tracks to the matching
            DNA strand so the opposite strand's zeros do not dilute the signal.
            None keeps all tracks (legacy behaviour).
        aggregation: How to collapse the tissue/track dimension. One of:
            - "median" (default): robust to tissue outliers.
            - "mean": arithmetic mean of raw 0–1 values.
            - "logit_mean": sigmoid(mean(logit(x))) — additive in log-odds space.

    Returns:
        DataFrame with columns donor_prob, acceptor_prob, seq.
        Index is 0-based position in input sequence (non-padded positions only).
    """
    if sense_strand not in ("+", "-", None):
        raise ValueError(f"sense_strand must be '+', '-', or None; got {sense_strand!r}")
    if aggregation not in ("mean", "median", "logit_mean"):
        raise ValueError(
            f"aggregation must be one of 'mean', 'median', 'logit_mean'; got {aggregation!r}"
        )
    if tissue is None:
        warnings.warn(
            "tissue=None averages all AlphaGenome tissue tracks on the selected "
            f"strand (sense_strand={sense_strand!r}) using aggregation={aggregation!r}. "
            "For stranded genes consider specifying a tissue CURIE "
            "(e.g. 'UBERON:0000955' for brain).",
            UserWarning,
            stacklevel=2,
        )

    seq_len = len(sequence)
    if target_length is None:
        target_length = _choose_target_length(seq_len)

    padded_seq, left_pad = _pad_sequence_centered(sequence, target_length)

    otype = OutputType.SPLICE_SITE_USAGE if output_type == "usage" else OutputType.SPLICE_SITES

    # Optionally restrict to a specific tissue
    ontology_terms = None
    if tissue is not None:
        try:
            # Use from_curie() — OntologyTerm(curie_str) has wrong signature
            ontology_terms = [ontology.from_curie(tissue)]
        except Exception:
            try:
                ontology_terms = [OntologyTerm(tissue)]
            except Exception:
                ontology_terms = None

    result = client.predict_sequence(
        padded_seq,
        organism=Organism.HOMO_SAPIENS,
        requested_outputs=[otype],
        ontology_terms=ontology_terms,
    )

    data = result.splice_site_usage if output_type == "usage" else result.splice_sites

    # data.values shape: (L, n_tracks). SPLICE_SITE_USAGE returns one value per
    # position per track — donor vs acceptor is distinguished by *position*, not
    # by track index, so the same tracks are used for both.
    vals = np.array(data.values)
    orig_slice = vals[left_pad: left_pad + seq_len]  # (seq_len, n_tracks)

    # --- Strand filter --------------------------------------------------------
    # Prefer reading strand from metadata; fall back to the documented track
    # layout (tracks 0..366 = +, 367..733 = −) if metadata is unavailable.
    if sense_strand is not None and orig_slice.ndim == 2:
        strand_mask = None
        metadata = getattr(data, "metadata", None)
        if metadata is not None:
            try:
                if isinstance(metadata, pd.DataFrame) and "strand" in metadata.columns:
                    strand_mask = (metadata["strand"].values == sense_strand)
            except Exception:
                strand_mask = None
        if strand_mask is None:
            n_tracks = orig_slice.shape[1]
            half = n_tracks // 2
            if sense_strand == "+":
                strand_mask = np.zeros(n_tracks, dtype=bool)
                strand_mask[:half] = True
            else:
                strand_mask = np.zeros(n_tracks, dtype=bool)
                strand_mask[half:] = True
        if strand_mask.any():
            orig_slice = orig_slice[:, strand_mask]

    # Post-hoc tissue selection if ontology filtering failed
    if tissue is not None and ontology_terms is None and orig_slice.ndim == 2:
        try:
            metadata = getattr(data, "metadata", None)
            if metadata is not None:
                names = [
                    getattr(m, "name", None) or getattr(m, "ontology_term", None) or ""
                    for m in metadata
                ]
                matches = [i for i, n in enumerate(names) if tissue.lower() in str(n).lower()]
                if matches:
                    orig_slice = orig_slice[:, [matches[0]]]
        except Exception:
            pass

    # Donor and acceptor read the same tracks; position determines which.
    donor_tracks = orig_slice
    acceptor_tracks = orig_slice

    # --- Aggregate across tissue/track dimension ------------------------------
    def _aggregate(x: np.ndarray) -> np.ndarray:
        if x.ndim <= 1:
            return x
        if aggregation == "mean":
            return x.mean(axis=-1)
        if aggregation == "median":
            return np.median(x, axis=-1)
        # logit_mean
        eps = 1e-6
        xc = np.clip(x, eps, 1 - eps)
        l = np.log(xc / (1 - xc))
        m = l.mean(axis=-1)
        return 1.0 / (1.0 + np.exp(-m))

    donor_prob = _aggregate(donor_tracks)
    acceptor_prob = _aggregate(acceptor_tracks)

    donor_prob = np.clip(donor_prob, 0, 1)
    acceptor_prob = np.clip(acceptor_prob, 0, 1)

    if shift_to_intron:
        donor_prob = np.roll(donor_prob, 1)
        donor_prob[0] = 0
        acceptor_prob = np.roll(acceptor_prob, -1)
        acceptor_prob[-1] = 0

    idx = np.arange(seq_len)
    df = pd.DataFrame(
        {
            "donor_prob": donor_prob,
            "acceptor_prob": acceptor_prob,
            "seq": list(sequence),
        },
        index=idx,
    )
    return df


def make_alphagenome_predictor(
    client,
    tissue: str | None = None,
    output_type: str = "usage",
    target_length: int | None = None,
    shift_to_intron: bool = True,
    sense_strand: str | None = "+",
    aggregation: str = "median",
) -> Callable[[str], pd.DataFrame]:
    """
    Return a predictor function compatible with spliceO_predictions(predictor_fn=...).

    See `predict_splice_sites_alphagenome` for full parameter semantics, including
    `sense_strand` (default "+") and `aggregation` (default "median").

    Returns:
        Callable: sequence_str -> DataFrame(donor_prob, acceptor_prob, seq)
    """
    def predictor(sequence: str) -> pd.DataFrame:
        return predict_splice_sites_alphagenome(
            sequence,
            client,
            tissue=tissue,
            output_type=output_type,
            target_length=target_length,
            shift_to_intron=shift_to_intron,
            sense_strand=sense_strand,
            aggregation=aggregation,
        )
    return predictor


# Common tissue ontology CURIEs for reference
TISSUE_CURIES = {
    "brain": "UBERON:0000955",
    "liver": "UBERON:0002107",
    "heart": "UBERON:0000948",
    "skeletal_muscle": "UBERON:0001134",
    "testis": "UBERON:0000473",
    "kidney": "UBERON:0002113",
    "lung": "UBERON:0002048",
    "blood": "UBERON:0000178",
    "ovary": "UBERON:0000992",
    "thyroid": "UBERON:0002046",
}
