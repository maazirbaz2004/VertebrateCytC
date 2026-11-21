"""Phase 4 pytest suite: validate conservation metrics, coordinate map, anchors, and outputs."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from cytc.utils import REPO_ROOT

R = REPO_ROOT / "results/conservation"
T = REPO_ROOT / "results/tables"
F = REPO_ROOT / "results/figures"


def test_per_position_shapes_and_ranges():
    for fname in ("per_position_vert.csv", "per_position_all.csv"):
        df = pd.read_csv(R / fname)
        assert (df.identity_frac.between(0, 1)).all(), fname
        assert (df.shannon_entropy_bits >= -1e-9).all(), fname
        assert (df.shannon_entropy_bits <= math.log2(20) + 1e-6).all(), fname
        assert (df.norm_entropy.between(0, 1 + 1e-9)).all(), fname
        assert (df.jsd_score.between(0, 1 + 1e-9)).all(), fname


def test_invariant_means_zero_entropy():
    df = pd.read_csv(R / "per_position_vert.csv")
    inv = df[df["class"] == "invariant"]
    assert len(inv) > 0 and np.allclose(inv.shannon_entropy_bits, 0)


def test_heme_positions_conserved():
    pm = pd.read_csv(T / "position_map.csv")
    df = pd.read_csv(R / "per_position_vert.csv").merge(pm, on="aln_col")
    for classical in (14, 17, 18):  # Cys14, Cys17, His18
        row = df[df.classical_pos == classical].iloc[0]
        assert row["class"] == "invariant", classical
    # Also verify Met80 (classical 80, UniProt 81)
    met_row = df[df.classical_pos == 80].iloc[0]
    assert met_row["human_residue"] == "M" and met_row["class"] == "invariant"


def test_class_counts_and_outputs_exist():
    for fname in ("per_position_vert.csv", "per_position_all.csv"):
        df = pd.read_csv(R / fname)
        non_gap_rich = df[df["class"] != "gap_rich"]
        expected_sum = sum(
            (df["class"] == c).sum()
            for c in ("invariant", "highly_conserved", "moderately_conserved", "variable")
        )
        assert expected_sum == len(non_gap_rich), fname

    for p in (
        T / "position_map.csv",
        T / "conservation_summary.csv",
        T / "lineage_specific_substitutions.csv",
        R / "regions.csv",
        R / "summary.md",
        F / "fig_conservation_profile.png",
        F / "fig_conservation_profile.svg",
        F / "fig_entropy_bar.png",
        F / "fig_entropy_bar.svg",
        F / "fig_sequence_logo_vert.png",
        F / "fig_sequence_logo_vert.svg",
        F / "fig_sequence_logo_all.png",
        F / "fig_sequence_logo_all.svg",
        F / "fig_conservation_class_bar.png",
        F / "fig_conservation_class_bar.svg",
        F / "fig_conservation_heatmap_per_species.png",
        F / "fig_conservation_heatmap_per_species.svg",
    ):
        assert p.exists() and p.stat().st_size > 0, p

