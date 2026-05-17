"""
Pangolin utilities for splice site prediction.

Pangolin (Zeng & Li, Genome Biology 2022) predicts tissue-specific splice site usage.
This module wraps Pangolin to produce DataFrames compatible with spliceai_utils,
enabling drop-in use via spliceO_predictions(predictor_fn=...).

Usage:
    from my_utils import pangolin_utils
    models = pangolin_utils.load_pangolin_models()
    predictor = pangolin_utils.make_pangolin_predictor(models, tissue=None)
    # Pass predictor to spliceO_predictions(predictor_fn=predictor)

Pangolin API notes (from source):
    - 12 models total: 3 ensembles × 4 tissue groups
    - Loading: for i in [0,2,4,6]: for j in range(1,4):
      model = Pangolin(L, W, AR); load "models/final.{j}.{i}.3.v2"
    - Input: one-hot encoded (batch=1, 4, seq_len) — channels first
    - Output: (1, 12, seq_len_out), channels = [softmax2, sigmoid, softmax2, sigmoid, ...]
      Tissue group scores are at indices [1, 4, 7, 10] (2nd channel of each softmax pair)
    - CL = 10000 context length removed: output_len = input_len - 10000
      With pad_size=5000 (adds 10000), output_len = seq_len
"""

from __future__ import annotations

import os
import warnings
from typing import Callable

import numpy as np
import pandas as pd

# Tissue group names (indices 0-3 corresponding to i=0,2,4,6 in the loading loop)
# From Pangolin paper (Zeng & Li 2022, Genome Biology): heart, liver, brain, testis
TISSUE_NAMES = ["heart", "liver", "brain", "testis"]

# One-hot encoding lookup table matching Pangolin source
# A=1, C=2, G=3, T=4, N=0 → IN_MAP indexed
_IN_MAP = np.asarray(
    [
        [0, 0, 0, 0],  # N → all zeros
        [1, 0, 0, 0],  # A
        [0, 1, 0, 0],  # C
        [0, 0, 1, 0],  # G
        [0, 0, 0, 1],  # T
    ],
    dtype=np.float32,
)


def _one_hot_encode_pangolin(sequence: str) -> np.ndarray:
    """
    One-hot encode sequence matching Pangolin's encoding.

    Maps: A→1, C→2, G→3, T→4, N→0, then uses IN_MAP lookup.
    Returns array of shape (4, L) — channels first, for PyTorch Conv1d.
    """
    seq = sequence.upper()
    # Build integer array
    int_arr = np.zeros(len(seq), dtype=np.int8)
    for i, nt in enumerate(seq):
        if nt == "A":
            int_arr[i] = 1
        elif nt == "C":
            int_arr[i] = 2
        elif nt == "G":
            int_arr[i] = 3
        elif nt == "T":
            int_arr[i] = 4
        # else N → stays 0
    # (L, 4) then transpose to (4, L)
    return _IN_MAP[int_arr].T


def load_pangolin_models(model_dir: str | None = None) -> list:
    """
    Load all 12 Pangolin models (3 ensembles × 4 tissue groups).

    Args:
        model_dir: Directory containing Pangolin model weight files.
                   If None, uses the default location from the pangolin package.

    Returns:
        List of 12 loaded Pangolin model objects.
        The list is ordered as: [tissue0_ensemble1, tissue0_ensemble2, tissue0_ensemble3,
                                  tissue1_ensemble1, ..., tissue3_ensemble3]
        Tissue groups are: heart (idx 0), liver (idx 1), brain (idx 2), testis (idx 3)
        Access tissue group t (0-3): models[3*t : 3*t+3]

    Raises:
        ImportError: If pangolin is not installed.
        FileNotFoundError: If model weight files are not found.
    """
    try:
        import torch
        from pangolin.model import Pangolin, L, W, AR
    except ImportError as e:
        raise ImportError(
            "Pangolin is not installed. Install with:\n"
            "  conda run -n pangolin pip install git+https://github.com/tkzeng/Pangolin.git\n"
            "Then use the pangolin conda env to run code that calls this function."
        ) from e

    if model_dir is None:
        import pangolin as _pangolin_pkg
        model_dir = os.path.join(os.path.dirname(_pangolin_pkg.__file__), "models")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    models = []
    # Loading order matches Pangolin source: outer loop i in [0,2,4,6], inner j in range(1,4)
    for i in [0, 2, 4, 6]:
        for j in range(1, 4):
            model_path = os.path.join(model_dir, f"final.{j}.{i}.3.v2")
            if not os.path.exists(model_path):
                raise FileNotFoundError(
                    f"Model file not found: {model_path}\n"
                    f"Files in {model_dir}: {os.listdir(model_dir)[:10]}"
                )
            model = Pangolin(L, W, AR)
            if torch.cuda.is_available():
                model.cuda()
                weights = torch.load(model_path)
            else:
                weights = torch.load(model_path, map_location=torch.device("cpu"))
            model.load_state_dict(weights)
            model.eval()
            models.append(model)

    return models


def predict_splice_sites_pangolin(
    sequence: str,
    models: list,
    start_pos: int | None = None,
    end_pos: int | None = None,
    tissue: str | None = None,
    output: str = "usage",
    pad_size: int = 5000,
    shift_to_intron: bool = True,
) -> pd.DataFrame:
    """
    Predict splice sites using Pangolin.

    Args:
        sequence: DNA sequence string (ACGTN).
        models: List of 12 loaded Pangolin models from load_pangolin_models().
        start_pos: Start of output range (0-based, inclusive). Default: 0.
        end_pos: End of output range (0-based, exclusive). Default: len(sequence).
        tissue: One of TISSUE_NAMES ("heart", "liver", "brain", "testis"), or None
                to average across all 4 tissue groups.
        output: Currently unused; kept for API compatibility.
        pad_size: Padding added to each side of the sequence for context. Default 5000.
                  With pad_size=5000, Pangolin's CL=10000 is exactly cancelled so
                  the output length equals the input sequence length.
        shift_to_intron: If True, shift donor +1bp and acceptor -1bp to align with
                         intron boundaries (matches SpliceAI behavior). Default True.

    Returns:
        DataFrame with columns: donor_prob, acceptor_prob, seq
        Index is 0-based position in the input sequence (same format as
        predict_splice_sites() in spliceai_utils).

    Notes:
        Pangolin does not natively separate donor from acceptor splice sites.
        The splice site usage score is used for both donor_prob and acceptor_prob.
        When using spliceO_predictions(), only the relevant site type is extracted
        per site_of_interest, so this is functionally correct.
    """
    import torch

    if start_pos is None:
        start_pos = 0
    if end_pos is None:
        end_pos = len(sequence)

    # Pad with N's for context
    padded = "N" * pad_size + sequence + "N" * pad_size

    # One-hot encode: (4, padded_len) then unsqueeze to (1, 4, padded_len)
    x_np = _one_hot_encode_pangolin(padded)  # (4, padded_len)
    x_tensor = torch.from_numpy(x_np).unsqueeze(0)  # (1, 4, padded_len)

    if torch.cuda.is_available():
        x_tensor = x_tensor.cuda()

    # Run models and collect per-tissue-group scores
    # tissue_scores[t] = list of score arrays from 3 ensemble models for tissue group t
    tissue_scores = [[] for _ in range(4)]

    with torch.no_grad():
        for group_idx in range(4):
            group_models = models[3 * group_idx : 3 * group_idx + 3]
            out_channel = [1, 4, 7, 10][group_idx]  # sigma channel for this group
            for model in group_models:
                # output shape: (1, 12, output_len)
                pred = model(x_tensor)
                # Extract the splice site score channel for this tissue group
                scores = pred[0, out_channel, :].cpu().numpy()  # (output_len,)
                tissue_scores[group_idx].append(scores)

    # Average across the 3 ensembles per tissue group
    # tissue_avg[t] shape: (output_len,)
    tissue_avg = [np.mean(scores, axis=0) for scores in tissue_scores]

    # The output length = padded_len - CL where CL=10000
    # With pad_size=5000: output_len = len(sequence) + 10000 - 10000 = len(sequence)
    # Extract positions corresponding to the original (unpadded) sequence
    # The model removes CL//2 = 5000 from each end, so the output directly corresponds
    # to the original sequence positions (0 to len(sequence)-1)
    seq_len = len(sequence)

    # Select tissue or average
    if tissue is not None:
        tissue_lower = tissue.lower()
        matches = [i for i, t in enumerate(TISSUE_NAMES) if tissue_lower in t]
        if not matches:
            raise ValueError(
                f"Tissue '{tissue}' not found. Available: {TISSUE_NAMES}"
            )
        t_idx = matches[0]
        scores_arr = tissue_avg[t_idx]
    else:
        # Average across all 4 tissue groups
        scores_arr = np.mean(tissue_avg, axis=0)

    # Verify length matches sequence
    if len(scores_arr) != seq_len:
        # This can happen if pad_size is not exactly CL//2; trim/pad as needed
        if len(scores_arr) > seq_len:
            scores_arr = scores_arr[:seq_len]
        else:
            scores_arr = np.pad(
                scores_arr, (0, seq_len - len(scores_arr)), constant_values=0.0
            )

    scores_arr = np.clip(scores_arr, 0.0, 1.0)

    # Pangolin outputs a single splice-site score per position (no donor/acceptor split)
    # Use the same score for both; the caller selects only the relevant site type
    donor_prob = scores_arr.copy()
    acceptor_prob = scores_arr.copy()

    if shift_to_intron:
        # Shift donor right by 1 (first base of intron), acceptor left by 1 (last base of intron)
        donor_prob = np.roll(donor_prob, 1)
        donor_prob[0] = 0.0
        acceptor_prob = np.roll(acceptor_prob, -1)
        acceptor_prob[-1] = 0.0

    # Slice to requested range
    donor_slice = donor_prob[start_pos:end_pos]
    acceptor_slice = acceptor_prob[start_pos:end_pos]
    seq_slice = sequence[start_pos:end_pos]

    idx = np.arange(start_pos, end_pos)
    df = pd.DataFrame(
        {
            "donor_prob": donor_slice,
            "acceptor_prob": acceptor_slice,
            "seq": list(seq_slice),
        },
        index=idx,
    )
    return df


def make_pangolin_predictor(
    models: list,
    tissue: str | None = None,
    output: str = "usage",
    pad_size: int = 5000,
    shift_to_intron: bool = True,
) -> Callable[[str], pd.DataFrame]:
    """
    Create a predictor function compatible with spliceO_predictions(predictor_fn=...).

    Args:
        models: List of 12 loaded Pangolin models from load_pangolin_models().
        tissue: Tissue name ("heart", "liver", "brain", "testis") or None to average.
        output: "usage" or "probability" (for API compatibility; currently unused).
        pad_size: Context padding size. Default 5000.
        shift_to_intron: Align predictions to intron boundaries. Default True.

    Returns:
        Callable: sequence_str -> DataFrame(donor_prob, acceptor_prob, seq)

    Example:
        models = load_pangolin_models()
        predictor = make_pangolin_predictor(models, tissue="liver")
        results = spliceO_predictions(..., models=None, predictor_fn=predictor)
    """

    def predictor(sequence: str) -> pd.DataFrame:
        return predict_splice_sites_pangolin(
            sequence,
            models,
            tissue=tissue,
            output=output,
            pad_size=pad_size,
            shift_to_intron=shift_to_intron,
        )

    return predictor
