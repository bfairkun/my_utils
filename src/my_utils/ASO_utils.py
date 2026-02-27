"""
ASO off-target analysis utilities.

Two main functions:
  run_blast()              – BLAST a dict of {name: sequence} and return raw hits
  parse_blast_offtargets() – classify hits and return a per-ASO summary table
"""

import os
import subprocess
import tempfile

import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "stitle", "qseq",
]

_BLAST_DB_DEFAULT = (
    "/project2/yangili1/bjf79/ReferenceGenomes/GRCh38_GencodeRelease44Comprehensive/Blastdb/Hg38Genome_And_Gencodev49Transcripts"
)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_blast(aso_dict, blast_db=_BLAST_DB_DEFAULT, max_target_seqs=100,
              blastn_path="blastn"):
    """Run blastn on a dict of {name: sequence} and return raw hit table.

    Parameters
    ----------
    aso_dict : dict
        Mapping of ``name -> sequence`` for each oligo to BLAST.
        Names become the ``qseqid`` in the returned DataFrame.
    blast_db : str
        Path to the BLAST database (without file extension).
    max_target_seqs : int
        Passed to ``-max_target_seqs``.
    blastn_path : str
        Path to the ``blastn`` executable.  Defaults to ``"blastn"`` (i.e.
        uses whatever is on ``$PATH``).  Override with a full path when
        calling from a script that doesn't inherit the conda env PATH.

    Returns
    -------
    pd.DataFrame
        One row per BLAST hit with columns:
        qseqid, sseqid, pident, length, mismatch, gapopen,
        qstart, qend, sstart, send, evalue, bitscore, stitle, qseq.
        Columns match the blastn -outfmt 7 output with no post-processing beyond
        numeric casting.  sstart > send indicates a minus-strand hit (blastn
        convention).
    """
    fd, fasta_path = tempfile.mkstemp(suffix=".fa")
    try:
        with os.fdopen(fd, "w") as fh:
            for name, seq in aso_dict.items():
                fh.write(f">{name}\n{seq}\n")

        cmd = [
            blastn_path,
            "-db", blast_db,
            "-outfmt", "7 std stitle qseq",
            "-max_target_seqs", str(max_target_seqs),
            "-task", "blastn-short",
            "-query", fasta_path,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    finally:
        os.unlink(fasta_path)

    # Parse tabular output — skip comment lines that start with '#'
    lines = [l for l in result.stdout.splitlines() if l and not l.startswith("#")]
    if not lines:
        return pd.DataFrame(columns=_BLAST_COLS)

    blast_df = pd.DataFrame([l.split("\t") for l in lines], columns=_BLAST_COLS)

    # Cast numeric columns
    int_cols = ["length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send"]
    float_cols = ["pident", "evalue", "bitscore"]
    blast_df[int_cols] = blast_df[int_cols].astype(int)
    blast_df[float_cols] = blast_df[float_cols].astype(float)

    return blast_df


def parse_blast_offtargets(blast_df, walk_df, delta_tm=10,
                           FilterTranscriptomeHitsForMinusStrand=True):
    """Classify BLAST hits and build a per-ASO off-target summary.

    Expects ``blast_df["qseqid"]`` to be formatted as ``"{MaskStart}_{MaskEnd}"``
    (i.e. the caller used those keys when building the dict passed to
    :func:`run_blast`).

    For each unique ASO (identified by MaskStart / MaskEnd):

    * **best_hit_tm** – Tm of the best-scoring genome hit.
    * **on-target gene** – any transcriptome gene whose best-transcript Tm
      equals best_hit_tm (i.e. it carries the same target sequence).
    * **off-targets** – all other hits within *delta_tm* °C of best_hit_tm.

    Parameters
    ----------
    blast_df : pd.DataFrame
        Output of :func:`run_blast`.
    walk_df : pd.DataFrame
        Original walk-results table with columns ``MaskStart``, ``MaskEnd``,
        ``ASO_Sequence`` (and any other columns to carry through).
    delta_tm : float
        Hits with ``|Tm − best_hit_tm| <= delta_tm`` are considered
        potentially off-target.  Default: 10 °C.
    FilterTranscriptomeHitsForMinusStrand : bool
        If True (default), only transcriptome hits on the minus strand
        (``strand == "-"``, i.e. raw BLAST sstart > send) are counted as
        off-targets.  ASOs are antisense oligos, so genuine off-target binding
        to a transcript appears as a minus-strand hit against the transcript
        sequence.  Plus-strand hits (sense matches) are excluded.

    Returns
    -------
    pd.DataFrame
        One row per unique ASO.  Contains every column from *walk_df* (first
        occurrence for each MaskStart/MaskEnd) plus:

        * ``best_hit_tm``
        * ``n_genome_offtargets``           – unique genomic loci
        * ``n_transcriptome_offtargets``    – unique off-target genes
        * ``genome_offtarget_locs``         – comma-sep ``chr:start-stop`` strings
        * ``genome_offtarget_tms``          – comma-sep Tms (descending)
        * ``genome_offtarget_strands``      – comma-sep strands for each genome locus
        * ``transcriptome_offtarget_genes`` – comma-sep gene names
        * ``transcriptome_offtarget_tms``   – comma-sep Tms (descending)
        * ``transcriptome_offtarget_strands`` – comma-sep strands for each off-target gene
          (always ``"-"`` when FilterTranscriptomeHitsForMinusStrand is True)
    """
    hits = blast_df.copy()

    # --- Derive strand and normalise sstart/send to BED order ----------------
    # blastn reports sstart > send for minus-strand hits; we record this as
    # strand="-" then sort so sstart is always the smaller coordinate.
    hits["strand"] = (hits["send"] > hits["sstart"]).map({True: "+", False: "-"})
    hits[["sstart", "send"]] = pd.DataFrame(
        {"sstart": hits[["sstart", "send"]].min(axis=1),
         "send":   hits[["sstart", "send"]].max(axis=1)}
    )

    # --- Tm from query-aligned subsequence (gaps stripped) -------------------
    hits["tm"] = hits["qseq"].apply(lambda s: round(mt.Tm_NN(s.replace("-", "")), 2))

    # --- Classify hit source -------------------------------------------------
    # Transcriptome entries: sseqid is the full pipe-delimited FASTA header
    # (no spaces → BLAST treats the whole string as the ID).
    # Format: ENST...|ENSG...|OTTHUMG...|OTTHUMT...|transcript_name|gene_name|len|biotype|
    # Genome entries: sseqid = chromosome name (e.g. "chr2").
    hits["hit_source"] = hits["sseqid"].apply(
        lambda s: "transcriptome" if s.startswith("ENST") else "genome"
    )

    def _ensg_id(sseqid):
        fields = sseqid.split("|")
        return fields[1] if len(fields) > 1 else None

    def _gene_name(sseqid):
        fields = sseqid.split("|")
        return fields[5] if len(fields) > 5 else None

    tr_mask = hits["hit_source"] == "transcriptome"
    hits.loc[tr_mask, "ensg_id"]   = hits.loc[tr_mask, "sseqid"].apply(_ensg_id)
    hits.loc[tr_mask, "gene_name"] = hits.loc[tr_mask, "sseqid"].apply(_gene_name)

    ge_mask = hits["hit_source"] == "genome"
    hits.loc[ge_mask, "genome_loc"] = hits.loc[ge_mask].apply(
        lambda r: f"{r.sseqid}:{r.sstart}-{r.send}",
        axis=1,
    )

    # --- Recover MaskStart / MaskEnd from qseqid ----------------------------
    split = hits["qseqid"].str.split("_", n=1, expand=True)
    hits["MaskStart"] = split[0].astype(int)
    hits["MaskEnd"]   = split[1].astype(int)

    # --- On-target genome locus (exact coordinate match) --------------------
    hits["is_on_target_locus"] = (
        (hits["hit_source"] == "genome")
        & (hits["sstart"] == hits["MaskStart"])
        & (hits["send"]   == hits["MaskEnd"])
    )

    # --- best_hit_tm per query (best genome hit by bitscore) ----------------
    genome_hits = hits[hits["hit_source"] == "genome"].sort_values("bitscore", ascending=False)
    best_genome = (
        genome_hits.groupby("qseqid", as_index=False)
        .first()[["qseqid", "tm"]]
        .rename(columns={"tm": "best_hit_tm"})
    )
    hits = hits.merge(best_genome, on="qseqid", how="left")

    # --- Filter to within delta_tm ------------------------------------------
    hits_close = hits[abs(hits["tm"] - hits["best_hit_tm"]) <= delta_tm].copy()

    # --- Genome off-targets -------------------------------------------------
    # Best Tm per unique locus (take first after sorting by tm desc to also get strand),
    # excluding the on-target locus.
    genome_off = (
        hits_close[
            (hits_close["hit_source"] == "genome")
            & ~hits_close["is_on_target_locus"]
        ]
        .sort_values(["qseqid", "genome_loc", "tm"], ascending=[True, True, False])
        .groupby(["qseqid", "genome_loc"], as_index=False)
        .first()[["qseqid", "genome_loc", "tm", "strand"]]
        .sort_values(["qseqid", "tm"], ascending=[True, False])
    )
    genome_summary = (
        genome_off
        .groupby("qseqid")
        .agg(
            n_genome_offtargets=("genome_loc", "count"),
            genome_offtarget_locs=("genome_loc", ",".join),
            genome_offtarget_tms=("tm", lambda x: ",".join(str(v) for v in x)),
            genome_offtarget_strands=("strand", ",".join),
        )
        .reset_index()
    )

    # --- Transcriptome off-targets ------------------------------------------
    # Best Tm per gene; exclude genes whose best Tm == best_hit_tm (on-target gene)
    tx_hits = hits_close[hits_close["hit_source"] == "transcriptome"]
    if FilterTranscriptomeHitsForMinusStrand:
        tx_hits = tx_hits[tx_hits["strand"] == "-"]
    gene_best = (
        tx_hits
        .sort_values(["qseqid", "ensg_id", "gene_name", "tm"], ascending=[True, True, True, False])
        .groupby(["qseqid", "ensg_id", "gene_name"], as_index=False)
        .first()[["qseqid", "ensg_id", "gene_name", "tm", "strand"]]
        .rename(columns={"tm": "gene_best_tm"})
        .merge(best_genome, on="qseqid", how="left")
    )
    gene_best["is_on_target_gene"] = gene_best["gene_best_tm"] == gene_best["best_hit_tm"]

    tx_off = (
        gene_best[~gene_best["is_on_target_gene"]]
        .sort_values(["qseqid", "gene_best_tm"], ascending=[True, False])
    )
    tx_summary = (
        tx_off
        .groupby("qseqid")
        .agg(
            n_transcriptome_offtargets=("ensg_id", "count"),
            transcriptome_offtarget_genes=("gene_name", ",".join),
            transcriptome_offtarget_tms=("gene_best_tm", lambda x: ",".join(str(v) for v in x)),
            transcriptome_offtarget_strands=("strand", ",".join),
        )
        .reset_index()
    )

    # --- Build final summary ------------------------------------------------
    unique_asos = walk_df.drop_duplicates(subset=["MaskStart", "MaskEnd"]).copy()
    unique_asos["qseqid"] = (
        unique_asos["MaskStart"].astype(str) + "_" + unique_asos["MaskEnd"].astype(str)
    )

    summary = (
        unique_asos
        .merge(best_genome, on="qseqid", how="left")
        .merge(genome_summary, on="qseqid", how="left")
        .merge(tx_summary, on="qseqid", how="left")
    )

    summary["n_genome_offtargets"]        = summary["n_genome_offtargets"].fillna(0).astype(int)
    summary["n_transcriptome_offtargets"] = summary["n_transcriptome_offtargets"].fillna(0).astype(int)
    for col in ["genome_offtarget_locs", "genome_offtarget_tms", "genome_offtarget_strands",
                "transcriptome_offtarget_genes", "transcriptome_offtarget_tms",
                "transcriptome_offtarget_strands"]:
        summary[col] = summary[col].fillna("")

    return summary.drop(columns=["qseqid"])
