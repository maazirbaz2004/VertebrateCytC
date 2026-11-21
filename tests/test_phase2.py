"""Phase 2 pytest suite: validate curated FASTAs, metadata, QC tables, and figures."""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from cytc.utils import REPO_ROOT

P = REPO_ROOT / "data/processed"
T = REPO_ROOT / "results/tables"
F = REPO_ROOT / "results/figures"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def test_counts_and_unique_ids():
    recs = list(SeqIO.parse(P / "cytc_curated_treeid.fasta", "fasta"))
    raw_recs = list(SeqIO.parse(REPO_ROOT / "data/raw/cytc_raw.fasta", "fasta"))
    assert 12 <= len(recs) <= 20
    assert len(recs) == len(raw_recs), "No cross-species deduplication allowed"
    assert len({r.id for r in recs}) == len(recs)


def test_alphabet_and_motif():
    for fname in ("cytc_curated.fasta", "cytc_curated_treeid.fasta", "cytc_vertebrates.fasta"):
        for r in SeqIO.parse(P / fname, "fasta"):
            s = str(r.seq)
            assert set(s) <= set("ACDEFGHIKLMNPQRSTVWY"), f"{fname}:{r.id}"
            assert re.search(r"C..CH", s), f"{fname}:{r.id}"


def test_order_matches_id_map():
    ids = [r.id for r in SeqIO.parse(P / "cytc_curated_treeid.fasta", "fasta")]
    m = pd.read_csv(P / "id_map.csv").sort_values("order_index")
    assert ids == list(m.tree_id)

    curated_ids = [r.id.split("|")[0] for r in SeqIO.parse(P / "cytc_curated.fasta", "fasta")]
    assert curated_ids == ids

    qc = pd.read_csv(T / "seq_qc.csv")
    assert list(qc["species"]) == ids

    aa = pd.read_csv(T / "aa_composition.csv")
    assert list(aa["species"]) == ids


def test_outgroups_present_and_vertebrates_subset():
    m = pd.read_csv(P / "id_map.csv")
    assert {"Saccharomyces_cerevisiae", "Drosophila_melanogaster"} <= set(m.tree_id)

    vert_ids = [r.id for r in SeqIO.parse(P / "cytc_vertebrates.fasta", "fasta")]
    expected_verts = list(m.loc[m["role"] == "ingroup"].sort_values("order_index")["tree_id"])
    assert vert_ids == expected_verts
    assert "Saccharomyces_cerevisiae" not in vert_ids
    assert "Drosophila_melanogaster" not in vert_ids


def test_manifest_checksums_and_qc_outputs():
    manifest = json.loads((P / "dataset_manifest.json").read_text(encoding="utf-8"))
    for fname, info in manifest["files"].items():
        fpath = P / fname
        assert fpath.exists(), fname
        assert _sha256(fpath) == info["sha256"], f"Checksum mismatch for {fname}"

    for tbl in (
        "seq_qc.csv",
        "aa_composition.csv",
        "aa_composition_long.csv",
        "length_summary.csv",
        "identical_sequence_groups.csv",
    ):
        assert (T / tbl).exists(), tbl

    for stem in ("fig_length_by_species", "fig_aa_composition_heatmap", "fig_pI_vs_length"):
        for ext in (".png", ".svg"):
            p = F / f"{stem}{ext}"
            assert p.exists() and p.stat().st_size > 1000, p

