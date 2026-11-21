"""Phase 6 pytest suite: validate NJ, FastTree, IQ-TREE trees, rooting, and exports."""
from __future__ import annotations

import pandas as pd
from Bio import Phylo

from cytc.utils import REPO_ROOT

T = REPO_ROOT / "results/trees"
TAB = REPO_ROOT / "results/tables"
FIG = REPO_ROOT / "results/figures"
ID_MAP = pd.read_csv(REPO_ROOT / "data/processed/id_map.csv")
EXPECTED = set(ID_MAP.tree_id)
EXPECTED_VERT = set(ID_MAP.loc[ID_MAP["role"] == "ingroup", "tree_id"])


def leaves(path):
    return {t.name for t in Phylo.read(path, "newick").get_terminals()}


def test_ml_tree_has_all_leaves():
    assert leaves(T / "iqtree/cytc_all.treefile") == EXPECTED
    assert leaves(T / "iqtree/cytc_all.collapsed.nwk") == EXPECTED
    assert leaves(T / "iqtree/cytc_vert.treefile") == EXPECTED_VERT


def test_nj_tree_has_all_leaves():
    assert leaves(T / "nj/nj_poisson.nwk") == EXPECTED
    assert leaves(T / "nj/nj_blosum62.nwk") == EXPECTED


def test_fasttree_has_all_leaves():
    assert leaves(T / "fasttree/fasttree_lg_gamma.nwk") == EXPECTED


def test_outgroup_rooting():
    tree = Phylo.read(T / "iqtree/cytc_all.treefile", "newick")
    kids = [{t.name for t in c.get_terminals()} for c in tree.root.clades]
    assert {"Saccharomyces_cerevisiae"} in kids


def test_phase6_tables_itol_and_figures_exist():
    for p in (
        TAB / "iqtree_summary.csv",
        TAB / "ml_clade_support.csv",
        T / "sensitivity_summary.csv",
        T / "phylogeny_summary.md",
        T / "itol/cytc_ml_rooted.nwk",
        T / "itol/itol_colorstrip_clade.txt",
        T / "itol/itol_labels_common_names.txt",
        T / "itol/README_itol.md",
        FIG / "fig_tree_ml_rectangular.png",
        FIG / "fig_tree_ml_rectangular.svg",
        FIG / "fig_tree_nj_rectangular.png",
        FIG / "fig_tree_nj_rectangular.svg",
        FIG / "fig_tree_comparison_NJ_vs_ML.png",
        FIG / "fig_tree_comparison_NJ_vs_ML.svg",
        FIG / "fig_tree_with_identity_heatmap.png",
        FIG / "fig_tree_with_identity_heatmap.svg",
    ):
        assert p.exists() and p.stat().st_size > 0, p

