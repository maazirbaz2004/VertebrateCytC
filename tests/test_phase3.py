"""Phase 3 pytest suite: validate MSAs, ungapped fidelity, heme motif, and outputs."""
from __future__ import annotations

import re

import pandas as pd
from Bio import AlignIO, SeqIO

from cytc.utils import REPO_ROOT

A = REPO_ROOT / "results/alignments"
T = REPO_ROOT / "results/tables"
F = REPO_ROOT / "results/figures"
P = REPO_ROOT / "data/processed"


def test_equal_length_rows():
    aln = AlignIO.read(A / "primary_all.fasta", "fasta")
    assert len({len(r.seq) for r in aln}) == 1 and len(aln) >= 12

    vert_aln = AlignIO.read(A / "primary_vert.fasta", "fasta")
    assert len({len(r.seq) for r in vert_aln}) == 1 and len(vert_aln) == 15


def test_ungapped_sequences_match_curated():
    cur_all = {r.id: str(r.seq) for r in SeqIO.parse(P / "cytc_curated_treeid.fasta", "fasta")}
    for fname in ("all_clustalo.fasta", "all_muscle.fasta", "primary_all.fasta"):
        aln = AlignIO.read(A / fname, "fasta")
        for r in aln:
            assert str(r.seq).replace("-", "") == cur_all[r.id], f"{fname}:{r.id}"

    cur_vert = {r.id: str(r.seq) for r in SeqIO.parse(P / "cytc_vertebrates.fasta", "fasta")}
    for fname in ("vert_clustalo.fasta", "vert_muscle.fasta", "primary_vert.fasta"):
        aln = AlignIO.read(A / fname, "fasta")
        for r in aln:
            assert str(r.seq).replace("-", "") == cur_vert[r.id], f"{fname}:{r.id}"


def test_heme_motif_columns_consistent():
    for fname in ("primary_all.fasta", "primary_vert.fasta"):
        aln = AlignIO.read(A / fname, "fasta")
        starts = set()
        for r in aln:
            m = re.search(r"C-*.-*.-*C-*H", str(r.seq))
            assert m, f"{fname}:{r.id}"
            assert "-" not in m.group(0), f"gap inside heme motif in {fname}:{r.id}"
            starts.add(m.start())
        assert len(starts) == 1, f"inconsistent heme motif starts in {fname}: {starts}"


def test_fixed_order_and_expected_outputs():
    id_map = pd.read_csv(P / "id_map.csv").sort_values("order_index")
    expected_all = list(id_map["tree_id"])
    expected_vert = list(id_map.loc[id_map["role"] == "ingroup", "tree_id"])

    assert [r.id for r in AlignIO.read(A / "primary_all.fasta", "fasta")] == expected_all
    assert [r.id for r in AlignIO.read(A / "primary_vert.fasta", "fasta")] == expected_vert
    assert [r.id for r in AlignIO.read(A / "primary_all_trimmed.fasta", "fasta")] == expected_all

    for p in (
        A / "primary_all.aln",
        A / "primary_all.phy",
        A / "trimmed_columns.txt",
        A / "alignment_notes.md",
        T / "alignment_stats.csv",
        T / "gap_blocks.csv",
        T / "aligner_agreement.csv",
        F / "fig_alignment_overview.png",
        F / "fig_alignment_overview.svg",
    ):
        assert p.exists() and p.stat().st_size > 0, p

