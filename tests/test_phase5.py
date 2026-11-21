"""Phase 5 pytest suite: validate pairwise matrices, symmetry, bounds, long table, and figures."""
from __future__ import annotations

import numpy as np
import pandas as pd

from cytc.utils import REPO_ROOT

S = REPO_ROOT / "results/similarity"
T = REPO_ROOT / "results/tables"
F = REPO_ROOT / "results/figures"
P = REPO_ROOT / "data/processed"


def test_matrix_symmetric_and_diag():
    m = pd.read_csv(S / "identity_matrix.csv", index_col=0)
    assert (m.index == m.columns).all()
    assert np.allclose(m.values, m.values.T)
    assert np.allclose(np.diag(m.values), 100)
    assert ((m.values >= 0) & (m.values <= 100 + 1e-9)).all()

    for dist_name in ("p_distance_matrix.csv", "poisson_distance_matrix.csv", "n_substitutions_matrix.csv"):
        d = pd.read_csv(S / dist_name, index_col=0)
        assert (d.index == d.columns).all(), dist_name
        assert np.allclose(d.values, d.values.T), dist_name
        assert np.allclose(np.diag(d.values), 0), dist_name
        assert (d.values >= -1e-9).all(), dist_name


def test_similarity_geq_identity():
    i = pd.read_csv(S / "identity_matrix.csv", index_col=0)
    s = pd.read_csv(S / "similarity_matrix.csv", index_col=0)
    assert (s.values + 1e-9 >= i.values).all()


def test_long_table_complete_and_order():
    id_map = pd.read_csv(P / "id_map.csv").sort_values("order_index")
    m = pd.read_csv(S / "identity_matrix.csv", index_col=0)
    assert list(m.index) == list(id_map["tree_id"])

    n = len(m)
    long = pd.read_csv(S / "pairwise_long.csv")
    assert len(long) == n * (n - 1) // 2
    assert long["relationship"].notna().all()
    assert (long["n_compared"] > 0).all()


def test_expected_outputs_exist():
    for p in (
        S / "identity_matrix.csv",
        S / "similarity_matrix.csv",
        S / "n_substitutions_matrix.csv",
        S / "p_distance_matrix.csv",
        S / "poisson_distance_matrix.csv",
        S / "blosum62_score_matrix.csv",
        S / "pairwise_long.csv",
        S / "metric_definitions.md",
        S / "crosscheck.md",
        S / "summary.md",
        T / "similarity_group_summary.csv",
        T / "identity_to_human.csv",
        F / "fig_identity_heatmap.png",
        F / "fig_identity_heatmap.svg",
        F / "fig_similarity_heatmap.png",
        F / "fig_similarity_heatmap.svg",
        F / "fig_identity_clustermap.png",
        F / "fig_identity_clustermap.svg",
        F / "fig_identity_to_human.png",
        F / "fig_identity_to_human.svg",
        F / "fig_identity_vs_time.png",
        F / "fig_identity_vs_time.svg",
        F / "fig_identity_distribution_by_relationship.png",
        F / "fig_identity_distribution_by_relationship.svg",
    ):
        assert p.exists() and p.stat().st_size > 0, p
