"""Phase 7 pytest suite: validate reference tree, functional checklist, RF & clade tables, and RQ answers."""
from __future__ import annotations

import pandas as pd
from Bio import Phylo

from cytc.utils import REPO_ROOT

T = REPO_ROOT / "results/tables"
I = REPO_ROOT / "results/interpretation"
F = REPO_ROOT / "results/figures"


def test_reference_tree_leaves_match_dataset():
    ref = {t.name for t in Phylo.read(REPO_ROOT / "config/reference_species_tree.nwk", "newick").get_terminals()}
    data = set(pd.read_csv(REPO_ROOT / "data/processed/id_map.csv").tree_id)
    assert ref == data


def test_required_outputs_exist():
    for rel in [
        "results/tables/clade_recovery.csv",
        "results/tables/rf_vs_reference.csv",
        "results/tables/functional_checklist_results.csv",
        "results/tables/functional_sites_conservation.csv",
        "results/tables/variable_sites_analysis.csv",
        "results/interpretation/uniprot_features_human.csv",
        "results/interpretation/uniprot_comments.md",
        "results/interpretation/reference_tree_notes.md",
        "results/interpretation/research_questions_answers.md",
        "results/interpretation/limitations.md",
    ]:
        p = REPO_ROOT / rel
        assert p.exists() and p.stat().st_size > 0, rel

    for stem in (
        "fig_functional_sites_on_conservation",
        "fig_tree_vs_reference",
        "fig_clade_support_summary",
        "fig_rf_distance_summary",
    ):
        for ext in (".png", ".svg"):
            p = F / f"{stem}{ext}"
            assert p.exists() and p.stat().st_size > 0, p


def test_rq_file_answers_all_four():
    txt = (I / "research_questions_answers.md").read_text(encoding="utf-8")
    for q in ("RQ1", "RQ2", "RQ3", "RQ4", "Hypothesis"):
        assert q in txt

    lim = (I / "limitations.md").read_text(encoding="utf-8")
    for num in ("1.", "2.", "3.", "4.", "5.", "6.", "7.", "8.", "9."):
        assert num in lim
