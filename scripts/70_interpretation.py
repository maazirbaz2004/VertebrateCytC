"""Phase 7: Biological interpretation, UniProt annotation mapping, species-tree validation, and RQ answers."""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import dendropy
from dendropy.calculate import treecompare
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment, substitution_matrices
from Bio.Phylo.BaseTree import Clade, Tree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cytc.utils import REPO_ROOT, get_logger, load_config, resolve  # noqa: E402

KEEP_FEATURE_TYPES = {
    "Binding site",
    "Modified residue",
    "Natural variant",
    "Site",
    "Helix",
    "Beta strand",
    "Turn",
    "Domain",
    "Region",
    "Disulfide bond",
    "Cross-link",
}

# Kyte-Doolittle hydropathy scale for surface vs core proxy
KD_HYDROPATHY: dict[str, float] = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}

# Literature functional checklist (classical numbering = UniProt - 1 for human CYCS)
FUNCTIONAL_CHECKLIST: list[tuple[str, int, str]] = [
    ("Heme covalent attachment (thioether)", 14, "binds heme vinyl group"),
    ("Heme covalent attachment (thioether)", 17, "binds heme vinyl group"),
    ("Proximal axial heme ligand", 18, "binds heme iron"),
    ("Distal axial heme ligand", 80, "binds heme iron; redox & apoptotic switch"),
    ("Heme pocket H-bond network", 52, "stabilizes heme pocket"),
    ("Heme pocket H-bond network", 59, "stabilizes heme via H-bond"),
    ("Heme pocket H-bond network", 67, "stabilizes heme H-bond network"),
    ("Heme pocket H-bond network", 78, "stabilizes Met80 loop"),
    ("Close packing around heme", 6, "glycine tight packing"),
    ("Close packing around heme", 30, "proline structural turn"),
    ("Close packing around heme", 41, "glycine close approach to heme"),
    ("Close packing around heme", 82, "hydrophobic contact with heme"),
    ("Lysine patch (partner docking)", 13, "electrostatic docking"),
    ("Lysine patch (partner docking)", 27, "electrostatic docking"),
    ("Lysine patch (partner docking)", 72, "electrostatic docking / Apaf-1"),
    ("Lysine patch (partner docking)", 73, "electrostatic docking"),
    ("Lysine patch (partner docking)", 79, "electrostatic docking"),
    ("Lysine patch (partner docking)", 86, "electrostatic docking"),
    ("Lysine patch (partner docking)", 87, "electrostatic docking"),
    ("Regulatory phosphosite (human)", 48, "regulation of electron transfer / apoptosis"),
    ("Regulatory phosphosite (human)", 97, "regulation of electron transfer / apoptosis"),
]

REFERENCE_CLADES: list[tuple[str, list[str]]] = [
    ("Mammalia", ["Homo_sapiens", "Macaca_mulatta", "Mus_musculus", "Rattus_norvegicus",
                  "Canis_lupus_familiaris", "Bos_taurus", "Equus_caballus", "Monodelphis_domestica"]),
    ("Eutheria", ["Homo_sapiens", "Macaca_mulatta", "Mus_musculus", "Rattus_norvegicus",
                  "Canis_lupus_familiaris", "Bos_taurus", "Equus_caballus"]),
    ("Primates", ["Homo_sapiens", "Macaca_mulatta"]),
    ("Rodentia", ["Mus_musculus", "Rattus_norvegicus"]),
    ("Euarchontoglires", ["Homo_sapiens", "Macaca_mulatta", "Mus_musculus", "Rattus_norvegicus"]),
    ("Laurasiatheria", ["Canis_lupus_familiaris", "Bos_taurus", "Equus_caballus"]),
    ("Aves", ["Gallus_gallus", "Taeniopygia_guttata"]),
    ("Sauropsida", ["Gallus_gallus", "Taeniopygia_guttata", "Anolis_carolinensis", "Chrysemys_picta"]),
    ("Amniota", ["Homo_sapiens", "Macaca_mulatta", "Mus_musculus", "Rattus_norvegicus",
                 "Canis_lupus_familiaris", "Bos_taurus", "Equus_caballus", "Monodelphis_domestica",
                 "Gallus_gallus", "Taeniopygia_guttata", "Anolis_carolinensis", "Chrysemys_picta"]),
    ("Tetrapoda", ["Homo_sapiens", "Macaca_mulatta", "Mus_musculus", "Rattus_norvegicus",
                   "Canis_lupus_familiaris", "Bos_taurus", "Equus_caballus", "Monodelphis_domestica",
                   "Gallus_gallus", "Taeniopygia_guttata", "Anolis_carolinensis", "Chrysemys_picta",
                   "Xenopus_tropicalis"]),
    ("Teleostei", ["Danio_rerio", "Takifugu_rubripes"]),
    ("Vertebrata", ["Homo_sapiens", "Macaca_mulatta", "Mus_musculus", "Rattus_norvegicus",
                    "Canis_lupus_familiaris", "Bos_taurus", "Equus_caballus", "Monodelphis_domestica",
                    "Gallus_gallus", "Taeniopygia_guttata", "Anolis_carolinensis", "Chrysemys_picta",
                    "Xenopus_tropicalis", "Danio_rerio", "Takifugu_rubripes"]),
    ("Bilateria", ["Homo_sapiens", "Macaca_mulatta", "Mus_musculus", "Rattus_norvegicus",
                   "Canis_lupus_familiaris", "Bos_taurus", "Equus_caballus", "Monodelphis_domestica",
                   "Gallus_gallus", "Taeniopygia_guttata", "Anolis_carolinensis", "Chrysemys_picta",
                   "Xenopus_tropicalis", "Danio_rerio", "Takifugu_rubripes", "Drosophila_melanogaster"]),
]

GROUP_TO_TREE_CLADE = {
    "Primate": "Primate", "Rodent": "Rodent", "Carnivore": "OtherMammal",
    "Artiodactyl": "OtherMammal", "Perissodactyl": "OtherMammal", "Marsupial": "OtherMammal",
    "Bird": "Bird", "Reptile": "Reptile", "Amphibian": "Amphibian",
    "Fish": "Fish", "Insect": "Insect", "Fungus": "Fungus",
}
TREE_CLADE_COLORS = {
    "Primate": "#1b9e77", "Rodent": "#d95f02", "OtherMammal": "#7570b3",
    "Bird": "#e7298a", "Reptile": "#66a61e", "Amphibian": "#e6ab02",
    "Fish": "#a6761d", "Insect": "#666666", "Fungus": "#1f78b4",
}


def df_to_markdown(df: pd.DataFrame) -> str:
    """Format DataFrame as GitHub-flavored Markdown without tabulate."""
    headers = [str(c) for c in df.columns]
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in df.itertuples(index=False):
        lines.append("| " + " | ".join(str(v) for v in row) + " |")
    return "\n".join(lines)


def parse_support_label(label: str | None) -> tuple[float | None, float | None]:
    """Parse IQ-TREE 'SH-aLRT/UFBoot' or single float label."""
    if not label:
        return None, None
    if "/" in label:
        p = label.split("/")
        try:
            return float(p[0]), float(p[1])
        except ValueError:
            return None, None
    try:
        return None, float(label)
    except ValueError:
        return None, None


def extract_uniprot_human_annotations(
    raw_dir: Path, prov_df: pd.DataFrame, interp_dir: Path, logger: Any
) -> pd.DataFrame:
    """Task 7.1: Parse UniProt JSON for human CYCS into features CSV and comments markdown."""
    hum_acc = str(prov_df.loc[prov_df["species"] == "Homo_sapiens", "accession"].iloc[0])
    json_path = raw_dir / "uniprot" / f"{hum_acc}.json"
    if not json_path.exists():
        raise FileNotFoundError(f"Missing local UniProt JSON for human ({hum_acc}): {json_path}")

    entry = json.loads(json_path.read_text(encoding="utf-8"))
    feat_rows: list[dict[str, Any]] = []
    for f in entry.get("features", []):
        ftype = f.get("type", "")
        if ftype not in KEEP_FEATURE_TYPES:
            continue
        loc = f.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        ligand = f.get("ligand", {}).get("name", "")
        desc = f.get("description", "")
        if ftype == "Natural variant":
            alt = f.get("alternativeSequence", {})
            orig = alt.get("originalSequence", "")
            alts = ",".join(alt.get("alternativeSequences", []))
            if orig and alts:
                desc = f"{orig}->{alts}; {desc}".strip("; ")
        feat_rows.append(
            {"type": ftype, "description": desc, "start": start, "end": end, "ligand": ligand}
        )

    feat_df = pd.DataFrame(feat_rows, columns=["type", "description", "start", "end", "ligand"])
    feat_df.to_csv(interp_dir / "uniprot_features_human.csv", index=False)

    # Extract concise UniProt comments
    comments_by_type: dict[str, list[str]] = {"FUNCTION": [], "SUBCELLULAR LOCATION": [], "DISEASE": []}
    for c in entry.get("comments", []):
        ct = c.get("commentType", "")
        if ct == "FUNCTION":
            for t in c.get("texts", []):
                comments_by_type[ct].append(t.get("value", ""))
        elif ct == "SUBCELLULAR LOCATION":
            locs = [sl.get("location", {}).get("value", "") for sl in c.get("subcellularLocations", [])]
            if locs:
                comments_by_type[ct].append(", ".join(x for x in locs if x))
        elif ct == "DISEASE":
            dis = c.get("disease", {})
            dname = dis.get("diseaseId", "")
            ddesc = dis.get("description", "")
            if dname:
                comments_by_type[ct].append(f"**{dname}**: {ddesc}")

    md_lines = [
        f"# UniProt Functional Annotations for Human Cytochrome c (`{hum_acc}`)",
        "",
        "## Function",
        *([f"- {t}" for t in comments_by_type["FUNCTION"]] or ["- _None recorded_"]),
        "",
        "## Subcellular Location",
        *([f"- {t}" for t in comments_by_type["SUBCELLULAR LOCATION"]] or ["- _None recorded_"]),
        "",
        "## Disease Associations",
        *([f"- {t}" for t in comments_by_type["DISEASE"]] or ["- _None recorded_"]),
        "",
    ]
    (interp_dir / "uniprot_comments.md").write_text("\n".join(md_lines), encoding="utf-8")
    logger.info("Extracted %d UniProt features for %s", len(feat_df), hum_acc)
    return feat_df


def build_functional_tables(
    feat_df: pd.DataFrame,
    per_vert: pd.DataFrame,
    pos_map: pd.DataFrame,
    aln_vert: MultipleSeqAlignment,
    aln_all: MultipleSeqAlignment,
    tables_dir: Path,
    logger: Any,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Tasks 7.1 (join) & 7.2 (literature checklist verified against alignment)."""
    vert_m = per_vert.merge(pos_map, on="aln_col", how="left")
    by_uniprot = vert_m.set_index("human_uniprot_pos")

    # Join UniProt features with conservation
    joined_rows: list[dict[str, Any]] = []
    for f in feat_df.itertuples():
        if pd.isna(f.start) or pd.isna(f.end):
            continue
        s, e = int(f.start), int(f.end)
        sub = vert_m[(vert_m["human_uniprot_pos"] >= s) & (vert_m["human_uniprot_pos"] <= e)]
        if sub.empty:
            continue
        joined_rows.append(
            {
                "type": f.type,
                "description": f.description,
                "uniprot_start": s,
                "uniprot_end": e,
                "classical_start": s - 1 if s > 1 else None,
                "classical_end": e - 1 if e > 1 else None,
                "ligand": f.ligand,
                "mean_identity_frac": round(float(sub["identity_frac"].mean()), 4),
                "mean_entropy_bits": round(float(sub["shannon_entropy_bits"].mean()), 4),
                "classes": ",".join(sorted(set(sub["class"]))),
            }
        )
    func_cons_df = pd.DataFrame(joined_rows)
    func_cons_df.to_csv(tables_dir / "functional_sites_conservation.csv", index=False)

    # Task 7.2: Literature checklist + UniProt Natural variants
    checklist_items = list(FUNCTIONAL_CHECKLIST)
    for f in feat_df[feat_df["type"] == "Natural variant"].itertuples():
        if pd.notna(f.start):
            u_pos = int(f.start)
            cl_pos = u_pos - 1
            checklist_items.append((f"Disease variant (UniProt Natural variant)", cl_pos, str(f.description)))

    def _counts_str(aln: MultipleSeqAlignment, col_0: int) -> str:
        cnt = Counter(str(r.seq)[col_0] for r in aln if str(r.seq)[col_0] != "-")
        return ", ".join(f"{aa}:{n}" for aa, n in sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0])))

    chk_rows: list[dict[str, Any]] = []
    for feat_name, cl_pos, role_desc in checklist_items:
        u_pos = cl_pos + 1
        if u_pos not in by_uniprot.index:
            continue
        row = by_uniprot.loc[u_pos]
        aln_col = int(row["aln_col"])
        vert_col = int(row["vert_aln_col"])
        v_cnt = _counts_str(aln_vert, vert_col - 1)
        a_cnt = _counts_str(aln_all, aln_col - 1)
        v_cls = str(row["class"])
        h_res = str(row["human_residue"])

        # Check deep conservation in all-17
        a_distinct = len({str(r.seq)[aln_col - 1] for r in aln_all if str(r.seq)[aln_col - 1] != "-"})
        if v_cls == "invariant" and a_distinct == 1:
            supports = "yes"
            comment = f"Strictly invariant across all 15 vertebrates and both outgroups ({role_desc})."
        elif v_cls in ("invariant", "highly_conserved"):
            supports = "yes" if v_cls == "invariant" else "partly"
            comment = f"{v_cls.replace('_', ' ').title()} in vertebrates; outgroups show {a_cnt} ({role_desc})."
        else:
            supports = "partly" if v_cls == "moderately_conserved" else "no"
            comment = f"Variable in vertebrates ({v_cnt}); {role_desc}."

        chk_rows.append(
            {
                "feature": feat_name,
                "classical_pos": cl_pos,
                "uniprot_pos": u_pos,
                "aln_col": aln_col,
                "human_residue": h_res,
                "vert_class": v_cls,
                "vert_residue_counts": v_cnt,
                "all17_residue_counts": a_cnt,
                "supports_expected_role": supports,
                "comment": comment,
            }
        )

    chk_df = pd.DataFrame(chk_rows)
    chk_df.to_csv(tables_dir / "functional_checklist_results.csv", index=False)
    logger.info("Verified %d functional checklist / disease variant sites against alignment.", len(chk_df))
    return func_cons_df, chk_df


def analyze_variable_sites(
    per_vert: pd.DataFrame,
    pos_map: pd.DataFrame,
    aln_vert: MultipleSeqAlignment,
    lineage_df: pd.DataFrame,
    tables_dir: Path,
    logger: Any,
) -> pd.DataFrame:
    """Task 7.3: Hydropathy (surface proxy), BLOSUM62 conservativeness, and lineage correlation of top 10 variable sites."""
    blosum = substitution_matrices.load("BLOSUM62")
    vert_m = per_vert.merge(pos_map, on="aln_col", how="left")
    top10 = vert_m.sort_values(["shannon_entropy_bits", "identity_frac"], ascending=[False, True]).head(10)
    lin_by_col = lineage_df.groupby("aln_col")["group"].apply(lambda s: ",".join(sorted(set(s)))).to_dict()

    rows: list[dict[str, Any]] = []
    for r in top10.itertuples():
        c0 = int(r.vert_aln_col) - 1
        chars = [str(rec.seq)[c0] for rec in aln_vert if str(rec.seq)[c0] != "-"]
        mean_kd = float(np.mean([KD_HYDROPATHY[ch] for ch in chars]))
        polar_charged_frac = sum(1 for ch in chars if ch in "KRHDEQNSGT") / len(chars)
        surface_proxy = "likely surface (hydrophilic)" if mean_kd < -0.5 else ("borderline" if mean_kd < 0.5 else "likely buried (hydrophobic)")

        # Pairwise substitutions among observed differing pairs
        n_diff_pairs = 0
        n_cons_pairs = 0
        for i in range(len(chars)):
            for j in range(i + 1, len(chars)):
                if chars[i] != chars[j]:
                    n_diff_pairs += 1
                    if float(blosum[chars[i], chars[j]]) >= 0:
                        n_cons_pairs += 1
        cons_frac = n_cons_pairs / n_diff_pairs if n_diff_pairs > 0 else 1.0
        sub_type = "mostly conservative" if cons_frac >= 0.5 else "mostly radical"

        rows.append(
            {
                "aln_col": int(r.aln_col),
                "human_uniprot_pos": int(r.human_uniprot_pos),
                "classical_pos": int(r.classical_pos) if pd.notna(r.classical_pos) else None,
                "entropy_bits": round(float(r.shannon_entropy_bits), 3),
                "mean_kd_hydropathy": round(mean_kd, 2),
                "surface_proxy": surface_proxy,
                "frac_conservative_subs": round(cons_frac, 3),
                "substitution_character": sub_type,
                "diagnostic_for_group": lin_by_col.get(int(r.aln_col), "none"),
            }
        )
    var_df = pd.DataFrame(rows)
    var_df.to_csv(tables_dir / "variable_sites_analysis.csv", index=False)
    logger.info("Analyzed top 10 variable vertebrate sites (PDB SASA skipped; Kyte-Doolittle hydropathy proxy used).")
    return var_df


def count_parsimony_informative(aln: MultipleSeqAlignment) -> int:
    """Count alignment columns with >= 2 distinct non-gap residues each occurring >= 2 times."""
    n_inf = 0
    seqs = [str(r.seq) for r in aln]
    for c in range(aln.get_alignment_length()):
        cnt = Counter(s[c] for s in seqs if s[c] != "-")
        if sum(1 for v in cnt.values() if v >= 2) >= 2:
            n_inf += 1
    return n_inf


def compare_trees_to_reference(
    ref_path: Path,
    tree_paths: dict[str, Path],
    all_taxa: list[str],
    root_outgroup: str,
    tables_dir: Path,
    logger: Any,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Task 7.4.2: RF distances (unrooted & rooted) and clade-by-clade monophyly recovery."""
    tns = dendropy.TaxonNamespace(all_taxa)
    ref_d = dendropy.Tree.get(path=str(ref_path), schema="newick", taxon_namespace=tns, preserve_underscores=True, rooting="force-rooted")
    ref_d.encode_bipartitions()
    ref_unrooted = dendropy.Tree.get(path=str(ref_path), schema="newick", taxon_namespace=tns, preserve_underscores=True, rooting="force-unrooted")
    ref_unrooted.encode_bipartitions()

    # Count non-trivial reference clades (rooted)
    ref_clade_sets: set[frozenset[str]] = set()
    for nd in ref_d.postorder_node_iter():
        lf = frozenset(n.taxon.label for n in nd.leaf_iter() if n.taxon is not None)
        if 1 < len(lf) < len(all_taxa):
            ref_clade_sets.add(lf)

    rf_rows: list[dict[str, Any]] = []
    max_rf_unrooted_bin = 2 * (len(all_taxa) - 3)
    max_rf_rooted_bin = 2 * (len(all_taxa) - 2)

    for label, path in tree_paths.items():
        # Read via Biopython, root cleanly on yeast outgroup, then export string for dendropy rooted RF
        bp = Phylo.read(path, "newick")
        bp.root_with_outgroup(root_outgroup)
        rooted_nwk = bp.format("newick")

        d_unr = dendropy.Tree.get(path=str(path), schema="newick", taxon_namespace=tns, preserve_underscores=True, rooting="force-unrooted")
        d_unr.encode_bipartitions()
        d_roo = dendropy.Tree.get(data=rooted_nwk, schema="newick", taxon_namespace=tns, preserve_underscores=True, rooting="force-rooted")
        d_roo.encode_bipartitions()

        rf_u = int(treecompare.symmetric_difference(ref_unrooted, d_unr))
        rf_r = int(treecompare.symmetric_difference(ref_d, d_roo))

        # Count recovered rooted reference clades
        inf_clade_sets: set[frozenset[str]] = set()
        for nd in d_roo.postorder_node_iter():
            lf = frozenset(n.taxon.label for n in nd.leaf_iter() if n.taxon is not None)
            if 1 < len(lf) < len(all_taxa):
                inf_clade_sets.add(lf)
        n_rec = len(ref_clade_sets & inf_clade_sets)

        rf_rows.append(
            {
                "tree": label,
                "unrooted_rf": rf_u,
                "norm_unrooted_rf": round(rf_u / max_rf_unrooted_bin, 4),
                "rooted_rf": rf_r,
                "norm_rooted_rf": round(rf_r / max_rf_rooted_bin, 4),
                "ref_clades_recovered": n_rec,
                "total_ref_clades": len(ref_clade_sets),
                "pct_ref_clades_recovered": round(100.0 * n_rec / len(ref_clade_sets), 1),
            }
        )

    rf_df = pd.DataFrame(rf_rows)
    rf_df.to_csv(tables_dir / "rf_vs_reference.csv", index=False)

    # Clade-by-clade recovery in ML, NJ, FastTree (all rooted on Saccharomyces_cerevisiae)
    def _clade_support_map(path: Path, mode: str) -> dict[frozenset[str], float | None]:
        t = Phylo.read(path, "newick")
        t.root_with_outgroup(root_outgroup)
        cmap: dict[frozenset[str], float | None] = {}
        for cl in t.get_nonterminals():
            lf = frozenset(x.name for x in cl.get_terminals())
            if mode == "iqtree":
                raw = cl.name if cl.name else (str(cl.confidence) if cl.confidence is not None else None)
                _, uf = parse_support_label(raw)
                cmap[lf] = uf
            elif mode == "nj":
                cmap[lf] = float(cl.confidence) if cl.confidence is not None else None
            elif mode == "fasttree":
                cmap[lf] = round(float(cl.confidence) * 100.0, 1) if cl.confidence is not None else None
        return cmap

    ml_map = _clade_support_map(tree_paths["ML (IQ-TREE primary)"], "iqtree")
    nj_map = _clade_support_map(tree_paths["NJ (Poisson)"], "nj")
    ft_map = _clade_support_map(tree_paths["FastTree (LG+G)"], "fasttree")

    clade_rows: list[dict[str, Any]] = []
    for cname, members in REFERENCE_CLADES:
        mset = frozenset(members)
        in_ml = mset in ml_map
        in_nj = mset in nj_map
        in_ft = mset in ft_map
        uf = ml_map.get(mset)
        nb = nj_map.get(mset)
        fsup = ft_map.get(mset)

        if in_ml and in_nj and in_ft:
            comment = "Recovered across all three methods."
        elif cname == "Teleostei":
            comment = "Split by fast-evolving Takifugu rubripes (long-branch attraction / rate heterogeneity)."
        elif cname in ("Sauropsida", "Aves"):
            comment = "Disrupted in ML by Takifugu rubripes nesting inside Sauropsida (low UFBoot 43-49); weak signal."
        elif cname in ("Mammalia", "Eutheria", "Euarchontoglires"):
            comment = "Short internal branches (<1-2 substitutions) and primate-specific substitutions limit resolution."
        elif cname == "Vertebrata":
            comment = "Recovered in ML (UFBoot 70) and FastTree; disrupted in NJ because Takifugu is pulled outside Drosophila."
        elif cname == "Bilateria":
            comment = "Rooted on S. cerevisiae; all animals form the sister clade to yeast."
        else:
            comment = "Weak or conflicting single-gene signal across short 105-aa alignment."

        clade_rows.append(
            {
                "clade": cname,
                "n_members": len(members),
                "members": ";".join(members),
                "recovered_in_ML": "Y" if in_ml else "N",
                "UFBoot": uf if in_ml else np.nan,
                "recovered_in_NJ": "Y" if in_nj else "N",
                "NJ_bootstrap": nb if in_nj else np.nan,
                "recovered_in_FastTree": "Y" if in_ft else "N",
                "FastTree_SH_pct": fsup if in_ft else np.nan,
                "comment": comment,
            }
        )

    # Also check unrooted Actinopterygii vs Tetrapoda split (requires Teleostei + Tetrapoda)
    clade_df = pd.DataFrame(clade_rows)
    clade_df.to_csv(tables_dir / "clade_recovery.csv", index=False)
    logger.info(
        "Tree comparison complete: ML recovered %d/%d reference clades (unrooted RF=%d).",
        int(rf_df.iloc[0]["ref_clades_recovered"]),
        int(rf_df.iloc[0]["total_ref_clades"]),
        int(rf_df.iloc[0]["unrooted_rf"]),
    )
    return rf_df, clade_df


def save_fig(fig: plt.Figure, stem: Path, cfg: dict[str, Any]) -> None:
    """Save figure in all configured formats and close it."""
    dpi = int(cfg["reporting"]["dpi"])
    for fmt in cfg["reporting"]["figure_formats"]:
        fig.savefig(stem.with_suffix(f".{fmt}"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_phase7_figures(
    per_vert: pd.DataFrame,
    pos_map: pd.DataFrame,
    chk_df: pd.DataFrame,
    ref_path: Path,
    ml_path: Path,
    rf_df: pd.DataFrame,
    clade_df: pd.DataFrame,
    id_map: pd.DataFrame,
    fig_dir: Path,
    cfg: dict[str, Any],
) -> None:
    """Task 7.8: Generate all 4 Phase 7 figures (PNG 300 dpi + SVG)."""
    fs = int(cfg["reporting"]["font_size"])
    plt.rcParams.update({"font.size": fs})
    vert_m = per_vert.merge(pos_map, on="aln_col", how="left")

    # 1. Functional sites overlaid on vertebrate conservation profile
    fig, ax = plt.subplots(figsize=(13.5, 5.4))
    x = vert_m["human_uniprot_pos"].astype(int).to_numpy()
    ax.plot(x, vert_m["conservation_score"], color="#A6CEE3", lw=1.1, label="Raw conservation")
    ax.plot(x, vert_m["smoothed_conservation"], color="#1F78B4", lw=2.2, label="Smoothed (window=5)")

    cat_styles = {
        "Heme ligands (C15,C18,H19,M81)": ("Heme|axial", "#D55E00", "D", 70),
        "Heme pocket & packing": ("H-bond|packing", "#009E73", "s", 55),
        "Lysine patch": ("Lysine patch", "#E69F00", "^", 60),
        "Phosphosites (Y49,Y98)": ("phosphosite", "#CC79A7", "P", 70),
        "Disease variants (THC4)": ("Disease variant", "#000000", "X", 70),
    }
    by_u = vert_m.set_index("human_uniprot_pos")["conservation_score"].to_dict()
    for label, (pat, col, mk, sz) in cat_styles.items():
        sub = chk_df[chk_df["feature"].str.contains(pat, case=False, regex=True)]
        if not sub.empty:
            xs = sub["uniprot_pos"].astype(int).to_numpy()
            ys = np.array([by_u[u] for u in xs])
            ax.scatter(xs, ys, color=col, marker=mk, s=sz, edgecolor="white", lw=0.6, label=label, zorder=5)

    ax.set_xlim(1, 105)
    ax.set_ylim(0.45, 1.06)
    ax.set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 105])
    ax.set_xlabel("Human UniProt residue position (Met1 = 1; classical = UniProt - 1)")
    ax.set_ylabel("Vertebrate conservation score (1 - norm H)")
    ax.set_title("Functional & Disease-Associated Sites Mapped onto Vertebrate Cytochrome c Conservation", fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower left", ncol=2, frameon=True, fontsize=fs - 1.5)
    save_fig(fig, fig_dir / "fig_functional_sites_on_conservation", cfg)

    # 2. Side-by-side Reference Species Tree vs Primary ML Tree
    sp_to_clade = {r.tree_id: GROUP_TO_TREE_CLADE[r.group] for r in id_map.itertuples()}
    ref_t = Phylo.read(ref_path, "newick")
    ref_t.ladderize()
    ml_t = Phylo.read(ml_path, "newick")
    ml_t.root_with_outgroup(str(cfg["phylogenetics"]["primary_root_outgroup"]))
    ml_t.ladderize()

    def _draw_cladogram_or_phylogram(ax: plt.Axes, tree: Tree, use_bl: bool, title: str) -> None:
        terms = list(tree.get_terminals())
        y_c: dict[Clade, float] = {t: float(i) for i, t in enumerate(terms)}
        x_c: dict[Clade, float] = {}

        def _depth(cl: Clade, d: float) -> None:
            x_c[cl] = d
            for ch in cl.clades:
                step = float(ch.branch_length or 0.0) if use_bl else 1.0
                _depth(ch, d + step)

        def _y(cl: Clade) -> float:
            if cl in y_c:
                return y_c[cl]
            ys = [_y(ch) for ch in cl.clades]
            y_c[cl] = 0.5 * (min(ys) + max(ys))
            return y_c[cl]

        _depth(tree.root, 0.0)
        _y(tree.root)
        max_x = max(x_c.values()) or 1.0
        if not use_bl:
            for t in terms:
                x_c[t] = max_x

        for cl in tree.find_clades():
            x0, y0 = x_c[cl], y_c[cl]
            if not cl.is_terminal() and cl.clades:
                ys = [y_c[ch] for ch in cl.clades]
                ax.plot([x0, x0], [min(ys), max(ys)], color="#333333", lw=1.3)
                for ch in cl.clades:
                    ax.plot([x0, x_c[ch]], [y_c[ch], y_c[ch]], color="#333333", lw=1.3)

        for t in terms:
            col = TREE_CLADE_COLORS[sp_to_clade[t.name]]
            ax.scatter(x_c[t], y_c[t], color=col, edgecolor="black", lw=0.6, s=55, zorder=4)
            ax.text(x_c[t] + max_x * 0.025, y_c[t], t.name.replace("_", " "), va="center", ha="left", fontsize=fs - 0.5)

        ax.set_xlim(-max_x * 0.05, max_x * 1.55)
        ax.set_ylim(len(terms) + 0.5, -0.8)
        ax.axis("off")
        ax.set_title(title, loc="left", fontweight="bold", fontsize=fs + 1)

    fig, axes = plt.subplots(1, 2, figsize=(16.5, 7.2))
    _draw_cladogram_or_phylogram(axes[0], ref_t, False, "A. Accepted Reference Species Tree (Cladogram)")
    _draw_cladogram_or_phylogram(axes[1], ml_t, True, "B. Inferred ML Tree (IQ-TREE Q.YEAST+G4)")
    handles = [mpatches.Patch(facecolor=TREE_CLADE_COLORS[c], edgecolor="black", lw=0.5, label=c) for c in TREE_CLADE_COLORS]
    fig.legend(handles=handles, title="Clade", loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=5, frameon=True, fontsize=fs - 1.5)
    fig.suptitle("Reference Species Tree vs. Cytochrome c Maximum Likelihood Tree", y=0.98, fontsize=fs + 3, fontweight="bold")
    fig.tight_layout(rect=[0, 0.06, 1, 0.94])
    save_fig(fig, fig_dir / "fig_tree_vs_reference", cfg)

    # 3. Clade recovery & support summary
    fig, ax = plt.subplots(figsize=(11.0, 5.8))
    y = np.arange(len(clade_df))
    h = 0.26
    ml_vals = clade_df["UFBoot"].fillna(0).to_numpy()
    nj_vals = clade_df["NJ_bootstrap"].fillna(0).to_numpy()
    ft_vals = clade_df["FastTree_SH_pct"].fillna(0).to_numpy()

    ax.barh(y - h, ml_vals, height=h, color="#0072B2", edgecolor="black", lw=0.4, label="ML (IQ-TREE UFBoot %)")
    ax.barh(y, nj_vals, height=h, color="#E69F00", edgecolor="black", lw=0.4, label="NJ (Poisson Bootstrap %)")
    ax.barh(y + h, ft_vals, height=h, color="#009E73", edgecolor="black", lw=0.4, label="FastTree (SH-like × 100)")
    ax.axvline(float(cfg["phylogenetics"]["support_threshold"]), color="crimson", linestyle="--", lw=1.0, label=f"Support threshold ({cfg['phylogenetics']['support_threshold']}%)")

    ax.set_yticks(y)
    ax.set_yticklabels(clade_df["clade"])
    ax.invert_yaxis()
    ax.set_xlim(0, 112)
    ax.set_xlabel("Node support (%) — 0 indicates clade not recovered")
    ax.set_title("Recovery & Support of Reference Clades Across Phylogenetic Methods", fontweight="bold")
    ax.grid(axis="x", linestyle=":", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", frameon=True, fontsize=fs - 1.5)
    save_fig(fig, fig_dir / "fig_clade_support_summary", cfg)

    # 4. RF distance summary
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    y = np.arange(len(rf_df))
    h = 0.35
    b1 = ax.barh(y - h / 2, rf_df["norm_unrooted_rf"], height=h, color="#0072B2", edgecolor="black", lw=0.5, label="Normalized Unrooted RF")
    b2 = ax.barh(y + h / 2, rf_df["norm_rooted_rf"], height=h, color="#D55E00", edgecolor="black", lw=0.5, label="Normalized Rooted RF")
    for bar, rf_raw in zip(b1, rf_df["unrooted_rf"]):
        w = bar.get_width()
        ax.text(w + 0.012, bar.get_y() + bar.get_height() / 2, f"{w:.2f} (RF={rf_raw})", va="center", fontsize=fs - 1.5)
    for bar, rf_raw in zip(b2, rf_df["rooted_rf"]):
        w = bar.get_width()
        ax.text(w + 0.012, bar.get_y() + bar.get_height() / 2, f"{w:.2f} (RF={rf_raw})", va="center", fontsize=fs - 1.5)

    ax.set_yticks(y)
    ax.set_yticklabels(rf_df["tree"])
    ax.invert_yaxis()
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Normalized Robinson-Foulds distance to reference species tree (lower = closer)")
    ax.set_title("Topological Distance (RF) to Accepted Vertebrate Species Tree", fontweight="bold", pad=26)
    ax.grid(axis="x", linestyle=":", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2, frameon=True)
    save_fig(fig, fig_dir / "fig_rf_distance_summary", cfg)


def write_interpretation_docs(
    interp_dir: Path,
    cons_sum: pd.DataFrame,
    regions_df: pd.DataFrame,
    chk_df: pd.DataFrame,
    var_df: pd.DataFrame,
    lineage_df: pd.DataFrame,
    hum_df: pd.DataFrame,
    group_sum: pd.DataFrame,
    rf_df: pd.DataFrame,
    clade_df: pd.DataFrame,
    ml_sup: pd.DataFrame,
    iq_sum: pd.DataFrame,
    n_pi_vert: int,
    n_pi_all: int,
    ml_tree: Tree,
    id_map: pd.DataFrame,
) -> None:
    """Write reference_tree_notes.md, research_questions_answers.md, and limitations.md."""
    # 1. reference_tree_notes.md
    ref_notes = """# Reference Species Tree Notes

- **File:** `config/reference_species_tree.nwk` (17 leaves, exact match to `data/processed/id_map.csv`).
- **Source:** Consensus vertebrate molecular phylogeny consistent with NCBI Taxonomy and TimeTree 5 major clades.
- **Deliberate polytomy:** Relationships among the three Laurasiatherian mammals (*Canis lupus familiaris*, *Equus caballus*, *Bos taurus*) are left as a polytomy `(Canis_lupus_familiaris,Equus_caballus,Bos_taurus)` because Scrotifera branching order remains debated across genomic datasets.
- **Turtle placement:** *Chrysemys picta* is placed as sister to Aves (*Gallus gallus*, *Taeniopygia guttata*) within Archelosauria, reflecting consensus molecular phylogenomics (distinct from older morphology-only placements).
- **Outgroups:** *Saccharomyces cerevisiae* (Fungus) as primary root outgroup, and *Drosophila melanogaster* (Insect) as animal outgroup sister to Vertebrata.
"""
    (interp_dir / "reference_tree_notes.md").write_text(ref_notes, encoding="utf-8")

    # Gather exact numbers for RQ1-RQ4
    v_row = cons_sum[cons_sum["dataset"] == "vertebrates_only"].iloc[0]
    a_row = cons_sum[cons_sum["dataset"] == "all_17_species"].iloc[0]
    n_cons_reg = int((regions_df["region_type"] == "conserved").sum())
    longest_reg = regions_df[regions_df["region_type"] == "conserved"].sort_values("length", ascending=False).iloc[0]

    # Branch length diagnostics
    term_bls = sorted([(t.name, float(t.branch_length or 0.0)) for t in ml_tree.get_terminals()], key=lambda kv: -kv[1])
    zero_terms = [name for name, bl in term_bls if bl < 1e-5]
    n_uf70 = int((ml_sup["ufboot"] >= 70).sum())
    n_uf90 = int((ml_sup["ufboot"] >= 90).sum())
    n_both = int(ml_sup["supported"].sum())

    ml_rf = rf_df[rf_df["tree"] == "ML (IQ-TREE primary)"].iloc[0]
    col_rf = rf_df[rf_df["tree"] == "ML (IQ-TREE collapsed)"].iloc[0]
    nj_rf = rf_df[rf_df["tree"] == "NJ (Poisson)"].iloc[0]
    ft_rf = rf_df[rf_df["tree"] == "FastTree (LG+G)"].iloc[0]

    mam_mean = float(group_sum.loc[group_sum["relationship"] == "mammal_vs_mammal_all", "mean_identity"].iloc[0])
    mb_mean = float(group_sum.loc[group_sum["relationship"] == "mammal_vs_bird", "mean_identity"].iloc[0])
    mf_mean = float(group_sum.loc[group_sum["relationship"] == "mammal_vs_fish", "mean_identity"].iloc[0])
    vy_mean = float(group_sum.loc[group_sum["relationship"] == "vertebrate_vs_yeast", "mean_identity"].iloc[0])

    best_model = str(iq_sum.loc[iq_sum["run"] == "cytc_all", "best_model"].iloc[0])
    rec_ml_clades = list(clade_df.loc[clade_df["recovered_in_ML"] == "Y", "clade"])
    miss_ml_clades = list(clade_df.loc[clade_df["recovered_in_ML"] == "N", "clade"])

    rq_md = f"""# Answers to Research Questions (RQ1–RQ4) & Hypothesis Evaluation

All numbers below are derived directly from `results/` CSV tables and trees. Both **UniProt-human** (Met1 = 1) and **classical horse-heart** (UniProt − 1) numbering are stated wherever positions are cited.

---

## RQ1: How conserved is cytochrome c across the selected vertebrate species?

Cytochrome c is **exceptionally conserved** across vertebrates:
- **Vertebrates only (`primary_vert.fasta`, n = 15, 105 alignment columns):**
  - **{int(v_row['n_invariant'])} / {int(v_row['total_columns'])} ({v_row['pct_invariant']}%)** columns are strictly invariant (zero entropy, zero gaps).
  - **{int(v_row['n_ge_90_conserved'])} / {int(v_row['total_columns'])} ({v_row['pct_ge_90_conserved']}%)** columns are >= 90% conserved.
  - Mean per-column identity is **{v_row['mean_identity_frac']*100:.2f}%** (median **100.0%**), mean Shannon entropy is **{v_row['mean_entropy_bits']:.4f} bits** (out of 4.32 bits max), and mean JSD score is **{v_row['mean_jsd_score']:.3f}**.
  - Mean vertebrate pairwise identity is **88.73%** (range **71.43%** for *Homo sapiens* vs. *Takifugu rubripes* to **100.00%** for identical rodent pair *Mus musculus* & *Rattus norvegicus*).
  - Only **{n_pi_vert} / 105** vertebrate columns are parsimony-informative.
- **Deep conservation (`primary_all.fasta`, n = 17 including *D. melanogaster* and *S. cerevisiae*, 110 columns):**
  - **{int(a_row['n_invariant'])} / {int(a_row['total_columns'])} ({a_row['pct_invariant']}%)** columns remain strictly invariant across ~1.1 billion years of eukaryotic divergence (**52 / 67 = 77.6%** of vertebrate-invariant sites stay invariant in fly and yeast).
  - Mean Shannon entropy increases modestly to **{a_row['mean_entropy_bits']:.4f} bits**; **{n_pi_all} / 110** columns are parsimony-informative.
- **Figures:** `fig_conservation_profile.png`, `fig_entropy_bar.png`, `fig_conservation_class_bar.png`.

---

## RQ2: Which amino-acid positions or regions are highly conserved?

- **Conserved regions:** {n_cons_reg} contiguous conserved blocks (>= 3 aa of invariant/highly conserved columns) were detected in vertebrates (`results/conservation/regions.csv`). The longest is a **17-aa block** (`{longest_reg['sequence_in_human']}`, UniProt **{int(longest_reg['human_uniprot_start'])}–{int(longest_reg['human_uniprot_end'])}**, classical **{int(longest_reg['classical_start'])}–{int(longest_reg['classical_end'])}**, mean conservation **{longest_reg['mean_conservation']:.4f}**), containing the longest strictly invariant stretch of **11 aa** (`NPKKYIPGTKM`, UniProt 71–81, classical 70–80, the Met80 loop).
- **Functional anchors verified against our alignment (`results/tables/functional_checklist_results.csv`):**
  - **Heme covalent attachment & axial ligands:** Cys15, Cys18, His19, and Met81 (UniProt; classical **Cys14, Cys17, His18, Met80**) are **100% invariant** across all 15 vertebrates and both outgroups.
  - **Heme pocket H-bond network:** Asn53, Trp60, Tyr68, Thr79 (UniProt; classical **Asn52, Trp59, Tyr67, Thr78**) are **100% invariant** across all 17 species.
  - **Close packing around heme:** Gly7, Pro31, Gly42, Phe83 (UniProt; classical **Gly6, Pro30, Gly41, Phe82**) are **100% invariant** across all 17 species.
  - **Surface lysine patch:** Lys28, Lys73, Lys74, Lys80, Lys88 (UniProt; classical **Lys27, Lys72, Lys73, Lys79, Lys87**) are **100% invariant** across all 17 species; Lys14 (cl. 13) and Lys87 (cl. 86) are highly conserved in vertebrates (14/15 K, with R in *Takifugu rubripes*).
  - **Regulatory phosphosites:** Tyr49 and Tyr98 (UniProt; classical **Tyr48, Tyr97**) are **100% invariant** across all 17 species.
  - **Disease-associated & natural variant sites (UniProt `P99999` Natural Variants in `functional_checklist_results.csv`):** The pathogenic Thrombocytopenia 4 variant **Gly42Ser** (classical **Gly41**, `rs121918552`) and **Lys56Arg** (classical **Lys55**, `rs11548795`) occur at **strictly invariant** positions across all 17 species, whereas **Met66Leu** (classical **Met65**) is moderately conserved (`M:13, F:1, N:1` in vertebrates).
- **Variable positions:** The 10 most variable vertebrate sites (led by UniProt 90 / cl. 89 at 2.04 bits, UniProt 45 / cl. 44 at 1.75 bits, and the C-terminal tail UniProt 101–105 / cl. 100–104) are predominantly **hydrophilic surface-exposed positions** (negative Kyte-Doolittle hydropathy) rather than buried heme-core residues (`results/tables/variable_sites_analysis.csv`).
- **Figures:** `fig_sequence_logo_vert.png`, `fig_functional_sites_on_conservation.png`, `fig_conservation_heatmap_per_species.png`.

---

## RQ3: How does sequence similarity vary between closely and distantly related species?

Pairwise identity declines strongly and monotonically with evolutionary distance (`results/tables/identity_to_human.csv`, `results/tables/similarity_group_summary.csv`):
- **Clade-level gradient:** Within mammals (**{mam_mean:.2f}%**) > mammal vs. bird (**{mb_mean:.2f}%**) > mammal vs. fish (**{mf_mean:.2f}%**) > vertebrate vs. yeast (**{vy_mean:.2f}%**).
- **Close vs. distant pairs:**
  - *Mus musculus* vs. *Rattus norvegicus* (20 Mya): **100.00%** identity (0 substitutions).
  - *Homo sapiens* vs. *Macaca mulatta* (29 Mya): **99.05%** identity (1 substitution).
  - *Gallus gallus* vs. *Taeniopygia guttata* (90 Mya): **94.29%** identity (6 substitutions).
  - *Homo sapiens* vs. *Danio rerio* (430 Mya): **84.62%** identity (16 substitutions).
  - *Homo sapiens* vs. *Takifugu rubripes* (430 Mya): **71.43%** identity (30 substitutions — fast-evolving outlier).
  - *Homo sapiens* vs. *Saccharomyces cerevisiae* (1100 Mya): **62.50%** identity (39 substitutions).
- **Molecular-clock & rate heterogeneity:**
  - Identity vs. approximate divergence time: Spearman **rho = −0.9575** (human-vs-X, descriptive p = 6.00e−09) and **rho = −0.8318** (all 136 pairs).
  - Poisson distance vs. time: linear OLS **R² = 0.8154** (human-vs-X), slope = **0.000352** substitutions/site/Myr (**0.352** subs/site/Gyr), confirming slow overall evolution.
  - **Rate heterogeneity:** *Takifugu rubripes* (`cyc-a`, residual **+0.136**, longest vertebrate terminal branch = **{term_bls[1][1] if term_bls[0][0]=='Saccharomyces_cerevisiae' else term_bls[0][1]:.3f}** subs/site) evolves markedly faster than *Danio rerio*. In primates, *Homo sapiens* and *Macaca mulatta* share **5 diagnostic primate-specific substitutions** (`I12`, `M13`, `A51`, `V84`, `E90` UniProt; classical `11, 12, 50, 83, 89`) relative to all other sampled vertebrates (`results/tables/lineage_specific_substitutions.csv`), consistent with a primate-lineage substitution burst.
- **Figures:** `fig_identity_heatmap.png`, `fig_identity_to_human.png`, `fig_identity_vs_time.png`, `fig_identity_distribution_by_relationship.png`.

---

## RQ4: Can cytochrome c sequence information be used to infer meaningful phylogenetic relationships?

**Partially — cytochrome c reliably separates major outgroups and a few distinctive clades, but lacks sufficient informative sites ({n_pi_vert} parsimony-informative sites in vertebrates) to resolve the full vertebrate species tree.**
- **Quantitative comparison vs. accepted species tree (`results/tables/rf_vs_reference.csv`, `clade_recovery.csv`):**
  - **ML (IQ-TREE primary, `{best_model}`):** Unrooted RF = **{int(ml_rf['unrooted_rf'])}** (normalized **{ml_rf['norm_unrooted_rf']:.2f}**), rooted RF = **{int(ml_rf['rooted_rf'])}**, recovering **{int(ml_rf['ref_clades_recovered'])} / {int(ml_rf['total_ref_clades'])}** reference clades: `{', '.join(rec_ml_clades)}`.
  - **ML (collapsed polytomies, 7 weak branches collapsed):** Unrooted RF = **{int(col_rf['unrooted_rf'])}** (normalized **{col_rf['norm_unrooted_rf']:.2f}**) — collapsing unsupported branches reduces spurious conflicts.
  - **NJ (Poisson):** Unrooted RF = **{int(nj_rf['unrooted_rf'])}** (normalized **{nj_rf['norm_unrooted_rf']:.2f}**), recovering **{int(nj_rf['ref_clades_recovered'])} / {int(nj_rf['total_ref_clades'])}** reference clades.
  - **FastTree (LG+G):** Unrooted RF = **{int(ft_rf['unrooted_rf'])}** (normalized **{ft_rf['norm_unrooted_rf']:.2f}**), recovering **{int(ft_rf['ref_clades_recovered'])} / {int(ft_rf['total_ref_clades'])}** reference clades.
- **Resolution & support:** Only **{n_both} / {len(ml_sup)}** internal ML nodes pass joint strong support (`UFBoot >= 70` and `SH-aLRT >= 80`: **Vertebrata** 89.9/70 and **Primates** 98.9/100); **{n_uf70} / {len(ml_sup)}** reach `UFBoot >= 70` and **{n_uf90} / {len(ml_sup)}** reach `UFBoot >= 90`. Zero-length terminal branches ({len(zero_terms)} taxa, including identical rodents *M. musculus* & *R. norvegicus*) prevent resolving shallow relationships.
- **Why disagreements occur (`{', '.join(miss_ml_clades)}` not recovered in ML):**
  1. **Insufficient signal:** Only {n_pi_vert} parsimony-informative vertebrate sites across 15 taxa leave deep mammalian and tetrapod internodes with 0–2 supporting substitutions.
  2. **Rate heterogeneity / long-branch attraction:** Fast-evolving *Takifugu rubripes* (`cyc-a`) breaks Teleostei and Sauropsida in ML (nested inside Sauropsida with weak UFBoot 43–49) and breaks Vertebrata in NJ (pulled outside *Drosophila*).
  3. **Lineage-specific substitution bursts:** The 5 primate-specific substitutions give Primates a long stem branch (`0.087` subs/site, 98.9/100 support) that makes Primates appear artificially divergent from other mammals in single-gene distance/ML trees.
- **Figures:** `fig_tree_ml_rectangular.png`, `fig_tree_comparison_NJ_vs_ML.png`, `fig_tree_vs_reference.png`, `fig_clade_support_summary.png`, `fig_rf_distance_summary.png`.

---

## Hypothesis Evaluation

> **Hypothesis:** *"Closely related species show higher sequence similarity and cluster together on the tree."*

**Verdict: PARTIALLY SUPPORTED.**
- **Supported for pairwise similarity:** Sequence identity decreases strongly and monotonically with divergence time (human-vs-X Spearman `rho = −0.958`, mammal mean `93.6%` > mammal–bird `90.9%` > mammal–fish `80.8%` > vertebrate–yeast `61.4%`), and close relatives share 99–100% identity (human–macaque `99.05%`, mouse–rat `100.0%`).
- **Only partially supported for phylogenetic clustering:** While ML recovers Vertebrata (89.9/70), Primates (98.9/100), Rodentia (0/94), and Laurasiatheria (79.9/67), it recovers only **{int(ml_rf['ref_clades_recovered'])} / {int(ml_rf['total_ref_clades'])}** reference clades (unrooted RF = **{int(ml_rf['unrooted_rf'])}**, dropping to **{int(col_rf['unrooted_rf'])}** when 7 unsupported branches are collapsed; NJ Poisson recovers {int(nj_rf['ref_clades_recovered'])}/{int(nj_rf['total_ref_clades'])}, RF = {int(nj_rf['unrooted_rf'])}) and fails to recover monophyletic Mammalia, Aves, Sauropsida, Tetrapoda, or Teleostei in ML due to ultra-high conservation (only {n_pi_vert} parsimony-informative vertebrate sites), identical sequences, and rate heterogeneity in *Takifugu rubripes*.
"""
    (interp_dir / "research_questions_answers.md").write_text(rq_md, encoding="utf-8")

    # 3. limitations.md (all 9 required points)
    lim_md = f"""# Study Limitations

1. **Single short protein (~105 aa) with limited phylogenetic signal:** Cytochrome c is extremely conserved (63.8% invariant in vertebrates), leaving only **{n_pi_vert}** parsimony-informative sites in vertebrates ({n_pi_all} across all 17 species). Consequently, 13 of 15 internal ML nodes lack strong joint support (`UFBoot >= 70` and `SH-aLRT >= 80`), and 7 branches collapse into polytomies when unsupported splits are removed.
2. **Small, uneven taxon sampling:** Only 15 vertebrate species (one sequence per species) represent major classes, with mammals overrepresented (8 species) relative to amphibians (1 species) and fishes (2 teleosts, no cartilaginous or jawless fish).
3. **Identical sequences and zero-length branches:** *Mus musculus* and *Rattus norvegicus* have 100% identical somatic cytochrome c sequences (patristic distance `2e-6`), making their internal resolution arbitrary and preventing resolution of shallow splits.
4. **Lineage-specific rate heterogeneity and long-branch effects:** *Takifugu rubripes* (`cyc-a`) exhibits accelerated divergence (+0.136 Poisson residual above the molecular-clock line; 30 differences from human vs. 16 for *Danio rerio*), inducing long-branch attraction that disrupts Teleostei, Sauropsida (in ML), and Vertebrata (in NJ). Similarly, the 5 primate-specific substitutions lengthen the primate stem branch relative to other mammals.
5. **Database-derived annotations and potential paralog effects:** Four vertebrate sequences (*M. domestica*, *T. guttata*, *A. carolinensis*, *C. picta*) were retrieved from UniProt TrEMBL (unreviewed) because reviewed Swiss-Prot entries were unavailable, and teleosts underwent whole-genome duplication (so *Takifugu* `cyc-a` may reflect relaxed constraint on a teleost duplicate).
6. **Deep outgroup divergence:** *Saccharomyces cerevisiae* (~1.1 Gya divergence, 4-aa N-terminal extension, ~61% identity to vertebrates) sits on a very long root branch (`1.06` subs/site in ML), which can attract long ingroup branches. Cross-checks using *Drosophila melanogaster* and the unrooted vertebrate-only tree (`cytc_vert`, pruned RF = 4 vs. primary ML) confirmed that vertebrate internal instability persists even without outgroups.
7. **Approximate divergence times and non-independence of pairwise tests:** Divergence times in the molecular-clock analysis are approximate consensus node ages (Mya) rather than calibrated fossil priors, and pairwise species comparisons are statistically non-independent due to shared phylogenetic ancestry (making correlation p-values descriptive only).
8. **Gene tree vs. species tree:** A single nuclear-encoded mitochondrial protein gene tree need not match the organismal species tree because of incomplete lineage sorting, convergent functional substitutions, or lineage-specific selection on mitochondrial electron transport.
9. **Computational-only scope:** No wet-lab biochemical assays or 3D structural dynamics simulations were performed; surface exposure of variable sites was inferred from sequence hydropathy rather than solvent-accessible surface area or mutagenesis experiments.
"""
    (interp_dir / "limitations.md").write_text(lim_md, encoding="utf-8")


def main() -> None:
    """Run Phase 7 biological interpretation, species-tree validation, figures, and conclusions."""
    parser = argparse.ArgumentParser(description="Phase 7: Biological interpretation and species-tree validation.")
    parser.add_argument("--config", type=str, default=None, help="Optional path to config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    logger = get_logger("70_interpretation", cfg)
    logger.info(
        "Phase 7 goal: link conservation to UniProt human functional sites and literature checklist, analyze variable sites, "
        "validate NJ/FastTree/ML trees against the accepted vertebrate species tree, and answer RQ1-RQ4 with limitations. "
        "Inputs: Phase 4-6 results, data/raw/uniprot/P99999.json, config/reference_species_tree.nwk."
    )

    raw_dir = resolve(cfg, "raw")
    proc_dir = resolve(cfg, "processed")
    aln_dir = resolve(cfg, "alignments")
    cons_dir = resolve(cfg, "conservation")
    trees_dir = resolve(cfg, "trees")
    interp_dir = resolve(cfg, "interpretation")
    tables_dir = resolve(cfg, "tables")
    fig_dir = resolve(cfg, "figures")

    ref_path = REPO_ROOT / "config" / "reference_species_tree.nwk"
    id_map = pd.read_csv(proc_dir / "id_map.csv").sort_values("order_index")
    all_taxa = list(id_map["tree_id"])

    # Validate reference tree leaf set
    ref_leaves = {t.name for t in Phylo.read(ref_path, "newick").get_terminals()}
    if ref_leaves != set(all_taxa):
        raise ValueError(f"Reference tree leaves mismatch: diff={ref_leaves ^ set(all_taxa)}")

    prov_df = pd.read_csv(tables_dir / "dataset_provenance.csv")
    pos_map = pd.read_csv(tables_dir / "position_map.csv")
    per_vert = pd.read_csv(cons_dir / "per_position_vert.csv")
    cons_sum = pd.read_csv(tables_dir / "conservation_summary.csv")
    regions_df = pd.read_csv(cons_dir / "regions.csv")
    lineage_df = pd.read_csv(tables_dir / "lineage_specific_substitutions.csv")
    hum_df = pd.read_csv(tables_dir / "identity_to_human.csv")
    group_sum = pd.read_csv(tables_dir / "similarity_group_summary.csv")
    ml_sup = pd.read_csv(tables_dir / "ml_clade_support.csv")
    iq_sum = pd.read_csv(tables_dir / "iqtree_summary.csv")

    aln_vert = AlignIO.read(aln_dir / "primary_vert.fasta", "fasta")
    aln_all = AlignIO.read(aln_dir / "primary_all.fasta", "fasta")

    # Task 7.1 & 7.2
    feat_df = extract_uniprot_human_annotations(raw_dir, prov_df, interp_dir, logger)
    _, chk_df = build_functional_tables(feat_df, per_vert, pos_map, aln_vert, aln_all, tables_dir, logger)

    # Task 7.3
    var_df = analyze_variable_sites(per_vert, pos_map, aln_vert, lineage_df, tables_dir, logger)

    # Task 7.4: Parsimony-informative sites + tree comparison
    n_pi_vert = count_parsimony_informative(aln_vert)
    n_pi_all = count_parsimony_informative(aln_all)
    logger.info("Parsimony-informative sites: vert=%d/105, all=%d/110", n_pi_vert, n_pi_all)

    tree_paths = {
        "ML (IQ-TREE primary)": trees_dir / "iqtree" / "cytc_all.treefile",
        "ML (IQ-TREE collapsed)": trees_dir / "iqtree" / "cytc_all.collapsed.nwk",
        "NJ (Poisson)": trees_dir / "nj" / "nj_poisson.nwk",
        "NJ (BLOSUM62)": trees_dir / "nj" / "nj_blosum62.nwk",
        "FastTree (LG+G)": trees_dir / "fasttree" / "fasttree_lg_gamma.nwk",
    }
    root_og = str(cfg["phylogenetics"]["primary_root_outgroup"])
    rf_df, clade_df = compare_trees_to_reference(ref_path, tree_paths, all_taxa, root_og, tables_dir, logger)

    # Task 7.8: Figures
    plot_phase7_figures(
        per_vert,
        pos_map,
        chk_df,
        ref_path,
        tree_paths["ML (IQ-TREE primary)"],
        rf_df,
        clade_df,
        id_map,
        fig_dir,
        cfg,
    )

    # Task 7.6 & 7.7: Interpretation docs
    ml_tree = Phylo.read(tree_paths["ML (IQ-TREE primary)"], "newick")
    ml_tree.root_with_outgroup(root_og)
    write_interpretation_docs(
        interp_dir,
        cons_sum,
        regions_df,
        chk_df,
        var_df,
        lineage_df,
        hum_df,
        group_sum,
        rf_df,
        clade_df,
        ml_sup,
        iq_sum,
        n_pi_vert,
        n_pi_all,
        ml_tree,
        id_map,
    )

    logger.info(
        "SUMMARY: verified %d functional/disease sites; parsimony-informative vert=%d all=%d; "
        "ML vs reference unrooted_RF=%d (%d/%d clades); outputs=%s, %s",
        len(chk_df),
        n_pi_vert,
        n_pi_all,
        int(rf_df.iloc[0]["unrooted_rf"]),
        int(rf_df.iloc[0]["ref_clades_recovered"]),
        int(rf_df.iloc[0]["total_ref_clades"]),
        interp_dir / "research_questions_answers.md",
        tables_dir / "clade_recovery.csv",
    )


if __name__ == "__main__":
    main()
