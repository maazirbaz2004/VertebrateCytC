"""Phase 6: Phylogenetic reconstruction (NJ, FastTree, IQ-TREE), sensitivity analysis, and visualization."""
from __future__ import annotations

import argparse
import copy
import math
import random
import re
import shlex
import subprocess
import sys
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
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Clade, Tree
from Bio.Phylo.Consensus import bootstrap_trees, get_support, majority_consensus
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceMatrix, DistanceTreeConstructor

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cytc.utils import REPO_ROOT, get_logger, load_config, resolve  # noqa: E402

GROUP_TO_TREE_CLADE: dict[str, str] = {
    "Primate": "Primate",
    "Rodent": "Rodent",
    "Carnivore": "OtherMammal",
    "Artiodactyl": "OtherMammal",
    "Perissodactyl": "OtherMammal",
    "Marsupial": "OtherMammal",
    "Bird": "Bird",
    "Reptile": "Reptile",
    "Amphibian": "Amphibian",
    "Fish": "Fish",
    "Insect": "Insect",
    "Fungus": "Fungus",
}

TREE_CLADE_ORDER = [
    "Primate", "Rodent", "OtherMammal", "Bird", "Reptile", "Amphibian", "Fish", "Insect", "Fungus"
]
TREE_CLADE_COLORS = {
    "Primate": "#1b9e77",
    "Rodent": "#d95f02",
    "OtherMammal": "#7570b3",
    "Bird": "#e7298a",
    "Reptile": "#66a61e",
    "Amphibian": "#e6ab02",
    "Fish": "#a6761d",
    "Insect": "#666666",
    "Fungus": "#1f78b4",
}


class PoissonDistanceCalculator(DistanceCalculator):
    """Custom Biopython DistanceCalculator computing gap-free pairwise Poisson distance (-ln(1-p))."""

    def __init__(self) -> None:
        super().__init__("identity")

    def get_distance(self, msa: MultipleSeqAlignment) -> DistanceMatrix:
        """Compute lower-triangular DistanceMatrix of Poisson distances for an MSA."""
        names = [r.id for r in msa]
        seqs = [str(r.seq) for r in msa]
        n = len(names)
        matrix: list[list[float]] = []
        for i in range(n):
            row: list[float] = []
            si = seqs[i]
            for j in range(i + 1):
                if i == j:
                    row.append(0.0)
                    continue
                sj = seqs[j]
                n_comp = 0
                n_ident = 0
                for a, b in zip(si, sj):
                    if a == "-" or b == "-":
                        continue
                    n_comp += 1
                    if a == b:
                        n_ident += 1
                p = 1.0 - (n_ident / n_comp) if n_comp > 0 else 0.9999
                if p <= 0.0:
                    dist = 0.0
                elif p >= 0.9999:
                    dist = -math.log(0.0001)
                else:
                    dist = -math.log(1.0 - p)
                row.append(dist)
            matrix.append(row)
        return DistanceMatrix(names, matrix)


def df_to_markdown(df: pd.DataFrame) -> str:
    """Format a DataFrame as a GitHub-flavored Markdown table without tabulate."""
    headers = [str(c) for c in df.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in df.itertuples(index=False):
        lines.append("| " + " | ".join(str(v) for v in row) + " |")
    return "\n".join(lines)


def run_cmd(cmd: list[str], logger: Any, stdin_text: str | None = None) -> subprocess.CompletedProcess[str]:
    """Run an external command without shell=True and fail loudly on error."""
    logger.info("CMD: %s", " ".join(cmd))
    return subprocess.run(cmd, input=stdin_text, check=True, capture_output=True, text=True)


def clamp_negative_branches(tree: Tree) -> int:
    """Clamp negative branch lengths to 0.0 and return count clamped."""
    n_clamped = 0
    for clade in tree.find_clades():
        if clade.branch_length is not None and clade.branch_length < 0.0:
            clade.branch_length = 0.0
            n_clamped += 1
    return n_clamped


def strip_inner_names(tree: Tree) -> None:
    """Clear auto-generated Biopython 'InnerN' names so Newick exports clean numeric confidences."""
    for clade in tree.get_nonterminals():
        if clade.name and clade.name.startswith("Inner"):
            clade.name = None


def normalize_confidence_pct(tree: Tree) -> None:
    """Ensure clade.confidence values are on the 0-100 percentage scale."""
    confs = [c.confidence for c in tree.get_nonterminals() if c.confidence is not None]
    if confs and max(confs) <= 1.0:
        for c in tree.get_nonterminals():
            if c.confidence is not None:
                c.confidence = round(float(c.confidence) * 100.0, 1)
    else:
        for c in tree.get_nonterminals():
            if c.confidence is not None:
                c.confidence = round(float(c.confidence), 1)


def is_monophyletic_pair(tree: Tree, taxa: list[str]) -> bool:
    """Check whether a set of taxa forms a monophyletic clade in a rooted tree."""
    target = set(taxa)
    for clade in tree.get_nonterminals():
        leaves = {t.name for t in clade.get_terminals()}
        if leaves == target:
            return True
    return False


def run_nj_analysis(
    aln_all: MultipleSeqAlignment,
    nj_dir: Path,
    cfg: dict[str, Any],
    logger: Any,
) -> tuple[Tree, Tree, int]:
    """Task 6.1: Build Poisson NJ (primary, with bootstrap) and BLOSUM62 NJ trees."""
    seed = int(cfg["project"]["random_seed"])
    n_boot = int(cfg["phylogenetics"]["bootstrap_replicates_nj"])
    root_outgroup = str(cfg["phylogenetics"]["primary_root_outgroup"])
    outgroup_pair = list(cfg["phylogenetics"]["outgroups"])

    # Primary NJ: Poisson distance
    poiss_calc = PoissonDistanceCalculator()
    poiss_ctor = DistanceTreeConstructor(poiss_calc, "nj")
    nj_poiss = poiss_ctor.build_tree(aln_all)
    n_clamped_p = clamp_negative_branches(nj_poiss)

    random.seed(seed)
    logger.info("Running NJ Poisson bootstrap: %d replicates (seed=%d)", n_boot, seed)
    boot_list = list(bootstrap_trees(aln_all, times=n_boot, tree_constructor=poiss_ctor))
    for bt in boot_list:
        clamp_negative_branches(bt)

    nj_poiss_sup = get_support(nj_poiss, boot_list)
    normalize_confidence_pct(nj_poiss_sup)
    strip_inner_names(nj_poiss_sup)
    nj_poiss_sup.root_with_outgroup(root_outgroup)
    nj_poiss_sup.ladderize()

    outgroups_mono = is_monophyletic_pair(nj_poiss_sup, outgroup_pair)
    logger.info(
        "NJ Poisson: clamped %d negative branch(es); outgroup pair %s monophyletic=%s",
        n_clamped_p,
        outgroup_pair,
        outgroups_mono,
    )

    cons_tree = majority_consensus(boot_list, cutoff=0.5)
    normalize_confidence_pct(cons_tree)
    strip_inner_names(cons_tree)
    cons_tree.root_with_outgroup(root_outgroup)
    cons_tree.ladderize()

    # Secondary NJ: BLOSUM62 distance
    blosum_calc = DistanceCalculator("blosum62")
    blosum_ctor = DistanceTreeConstructor(blosum_calc, "nj")
    nj_blosum = blosum_ctor.build_tree(aln_all)
    n_clamped_b = clamp_negative_branches(nj_blosum)
    strip_inner_names(nj_blosum)
    nj_blosum.root_with_outgroup(root_outgroup)
    nj_blosum.ladderize()
    logger.info("NJ BLOSUM62: clamped %d negative branch(es).", n_clamped_b)

    Phylo.write(nj_poiss_sup, nj_dir / "nj_poisson.nwk", "newick")
    Phylo.write(nj_blosum, nj_dir / "nj_blosum62.nwk", "newick")
    Phylo.write(cons_tree, nj_dir / "nj_bootstrap_consensus.nwk", "newick")
    return nj_poiss_sup, nj_blosum, n_clamped_p


def run_fasttree(
    aln_path: Path, ft_dir: Path, cfg: dict[str, Any], logger: Any
) -> Path:
    """Task 6.2: Run FastTree approximate ML cross-check."""
    ft_bin = cfg["tools"]["fasttree"]
    ft_flags = shlex.split(str(cfg["phylogenetics"]["fasttree_model_flags"]))
    out_nwk = ft_dir / "fasttree_lg_gamma.nwk"

    fasta_text = aln_path.read_text(encoding="utf-8")
    res = run_cmd([ft_bin, *ft_flags], logger, stdin_text=fasta_text)
    out_nwk.write_text(res.stdout.strip() + "\n", encoding="utf-8")

    # Parse version and log-likelihood from stderr
    ll_match = re.search(r"Gamma20LogLk\s+([-\d.]+)", res.stderr)
    ver_line = next((ln.strip() for ln in res.stderr.splitlines() if "FastTree" in ln), "FastTree")
    logger.info(
        "FastTree complete: %s | Gamma20LogLk=%s | output=%s",
        ver_line,
        ll_match.group(1) if ll_match else "n/a",
        out_nwk,
    )
    return out_nwk


def run_iqtree_job(
    aln_path: Path,
    prefix: Path,
    cfg: dict[str, Any],
    logger: Any,
    model: str | None = None,
    seed: int | None = None,
    outgroup: str | None = None,
) -> Path:
    """Run a single IQ-TREE job with -keep-ident, UFBoot, SH-aLRT, and -bnni."""
    iq_bin = cfg["tools"]["iqtree"]
    ph = cfg["phylogenetics"]
    mod = model or str(ph["model_selection"])
    sd = int(seed if seed is not None else cfg["project"]["random_seed"])
    n_threads = str(cfg["project"]["n_threads"])

    cmd = [
        iq_bin,
        "-s", str(aln_path),
        "-m", mod,
        "-B", str(ph["ufboot_replicates"]),
        "-alrt", str(ph["alrt_replicates"]),
        "-bnni",
        "-T", "AUTO",
        "--threads-max", n_threads,
        "-seed", str(sd),
        "-keep-ident",
        "-redo",
        "--prefix", str(prefix),
    ]
    if outgroup:
        cmd.extend(["-o", outgroup])

    run_cmd(cmd, logger)
    return Path(f"{prefix}.treefile")


def parse_iqtree_report(iqtree_file: Path, run_name: str, aln_name: str) -> dict[str, Any]:
    """Extract best-fit model, log-likelihood, BIC, tree length, and identical-sequence notes from .iqtree."""
    txt = iqtree_file.read_text(encoding="utf-8")

    def _find(pattern: str, default: str = "") -> str:
        m = re.search(pattern, txt)
        return m.group(1).strip() if m else default

    best_model = _find(r"Best-fit model according to BIC:\s*(.+)") or _find(r"Model of substitution:\s*(.+)")
    ll = _find(r"Log-likelihood of the tree:\s*([-\d.]+)")
    n_params = _find(r"Number of free parameters.*?:\s*(\d+)")
    bic = _find(r"Bayesian information criterion \(BIC\) score:\s*([-\d.]+)")
    tree_len = _find(r"Total tree length \(sum of branch lengths\):\s*([\d.]+)")

    ident_lines = [ln.strip() for ln in txt.splitlines() if "identical" in ln.lower()]
    ident_note = "; ".join(ident_lines) if ident_lines else "none noted in .iqtree"

    return {
        "run": run_name,
        "alignment": aln_name,
        "best_model": best_model,
        "log_likelihood": float(ll) if ll else np.nan,
        "n_parameters": int(n_params) if n_params else np.nan,
        "bic": float(bic) if bic else np.nan,
        "tree_length": float(tree_len) if tree_len else np.nan,
        "identical_seq_warnings": ident_note,
    }


def parse_support_label(label: str | None) -> tuple[float | None, float | None]:
    """Parse IQ-TREE 'SH-aLRT/UFBoot' node label into (sh_alrt, ufboot)."""
    if not label:
        return None, None
    if "/" in label:
        parts = label.split("/")
        try:
            return float(parts[0]), float(parts[1])
        except ValueError:
            return None, None
    try:
        return None, float(label)
    except ValueError:
        return None, None


def extract_ml_clade_support(
    treefile: Path, root_outgroup: str, support_thresh: float
) -> tuple[Tree, pd.DataFrame]:
    """Read IQ-TREE treefile, root cleanly on outgroup, and build ml_clade_support table."""
    tree = Phylo.read(treefile, "newick")
    tree.root_with_outgroup(root_outgroup)
    tree.ladderize()

    rows: list[dict[str, Any]] = []
    node_idx = 0
    for clade in tree.get_nonterminals():
        if clade is tree.root:
            continue
        node_idx += 1
        label = clade.name if clade.name else (str(clade.confidence) if clade.confidence is not None else None)
        sh_alrt, ufboot = parse_support_label(label)
        leaves = sorted(t.name for t in clade.get_terminals())
        supported = bool(
            ufboot is not None
            and sh_alrt is not None
            and ufboot >= support_thresh
            and sh_alrt >= 80.0
        )
        rows.append(
            {
                "node_id": f"N{node_idx:02d}",
                "n_leaves": len(leaves),
                "leaves": ";".join(leaves),
                "sh_alrt": sh_alrt,
                "ufboot": ufboot,
                "branch_length": round(float(clade.branch_length or 0.0), 6),
                "supported": supported,
            }
        )
    return tree, pd.DataFrame(rows)


def collapse_unsupported_branches(
    treefile: Path, out_path: Path, root_outgroup: str, min_len: float = 1e-6, min_ufboot: float = 50.0
) -> int:
    """Collapse internal branches with length < min_len or UFBoot < min_ufboot into polytomies."""
    tree = Phylo.read(treefile, "newick")
    tree.root_with_outgroup(root_outgroup)
    n_collapsed = 0

    def _collapse_recursive(clade: Clade) -> None:
        nonlocal n_collapsed
        # First recurse on children
        for child in list(clade.clades):
            if not child.is_terminal():
                _collapse_recursive(child)
        new_children: list[Clade] = []
        for child in clade.clades:
            if child.is_terminal():
                new_children.append(child)
                continue
            bl = float(child.branch_length or 0.0)
            label = child.name if child.name else (str(child.confidence) if child.confidence is not None else None)
            _, uf = parse_support_label(label)
            if bl < min_len or (uf is not None and uf < min_ufboot):
                n_collapsed += 1
                for grand in child.clades:
                    grand.branch_length = float(grand.branch_length or 0.0) + bl
                    new_children.append(grand)
            else:
                new_children.append(child)
        clade.clades = new_children

    _collapse_recursive(tree.root)
    tree.ladderize()
    Phylo.write(tree, out_path, "newick")
    return n_collapsed


def compute_sensitivity_rf(
    primary_path: Path,
    variants: list[tuple[str, Path, bool]],
    all_taxa: list[str],
    vert_taxa: list[str],
    outgroups: list[str],
) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Compute unrooted RF distances against primary ML tree and list stable vs unstable clades."""
    tns_all = dendropy.TaxonNamespace(all_taxa)
    prim_all = dendropy.Tree.get(path=str(primary_path), schema="newick", taxon_namespace=tns_all, preserve_underscores=True)
    prim_all.encode_bipartitions()

    tns_vert = dendropy.TaxonNamespace(vert_taxa)
    prim_vert = dendropy.Tree.get(path=str(primary_path), schema="newick", taxon_namespace=tns_all, preserve_underscores=True)
    prim_vert.prune_taxa_with_labels(outgroups)
    prim_vert_re = dendropy.Tree.get(data=prim_vert.as_string(schema="newick"), schema="newick", taxon_namespace=tns_vert, preserve_underscores=True)
    prim_vert_re.encode_bipartitions()

    def _clade_splits(dtree: dendropy.Tree, total_taxa: set[str]) -> set[frozenset[str]]:
        splits: set[frozenset[str]] = set()
        for edge in dtree.postorder_edge_iter():
            if edge.bipartition is None:
                continue
            leafset = frozenset(n.taxon.label for n in edge.head_node.leaf_iter() if n.taxon is not None)
            if 1 < len(leafset) < len(total_taxa) - 1:
                comp = frozenset(total_taxa - leafset)
                canon = leafset if (len(leafset) < len(comp) or (len(leafset) == len(comp) and sorted(leafset) < sorted(comp))) else comp
                splits.add(canon)
        return splits

    all_set = set(all_taxa)
    prim_splits = _clade_splits(prim_all, all_set)
    split_counts: dict[frozenset[str], int] = {s: 1 for s in prim_splits}
    n_full_trees = 1

    rows: list[dict[str, Any]] = []
    for name, path, is_vert_only in variants:
        if is_vert_only:
            t_var = dendropy.Tree.get(path=str(path), schema="newick", taxon_namespace=tns_vert, preserve_underscores=True)
            t_var.encode_bipartitions()
            rf = int(treecompare.symmetric_difference(prim_vert_re, t_var))
            max_rf = 2 * (len(vert_taxa) - 3)
            n_taxa = len(vert_taxa)
            note = "compared on 15-vertebrate pruned namespace"
        else:
            t_var = dendropy.Tree.get(path=str(path), schema="newick", taxon_namespace=tns_all, preserve_underscores=True)
            t_var.encode_bipartitions()
            rf = int(treecompare.symmetric_difference(prim_all, t_var))
            max_rf = 2 * (len(all_taxa) - 3)
            n_taxa = len(all_taxa)
            v_splits = _clade_splits(t_var, all_set)
            n_full_trees += 1
            for s in v_splits:
                split_counts[s] = split_counts.get(s, 0) + 1
            note = "identical unrooted topology" if rf == 0 else f"{rf // 2} differing split(s)"

        rows.append(
            {
                "variant": name,
                "n_taxa": n_taxa,
                "rf_distance_to_primary": rf,
                "max_rf": max_rf,
                "normalized_rf": round(rf / max_rf, 4) if max_rf > 0 else 0.0,
                "identical_topology": rf == 0,
                "notes": note,
            }
        )

    stable = [
        "{" + ", ".join(sorted(s)) + "}"
        for s, cnt in sorted(split_counts.items(), key=lambda kv: (-kv[1], len(kv[0])))
        if cnt == n_full_trees
    ]
    unstable = [
        f"({cnt}/{n_full_trees}) " + "{" + ", ".join(sorted(s)) + "}"
        for s, cnt in sorted(split_counts.items(), key=lambda kv: (-kv[1], len(kv[0])))
        if cnt < n_full_trees and s in prim_splits
    ]
    return pd.DataFrame(rows), stable, unstable


def draw_phylogram(
    ax: plt.Axes,
    tree: Tree,
    sp_to_clade: dict[str, str],
    support_mode: str,
    title: str,
    font_size: int,
    show_legend: bool = True,
) -> list[str]:
    """Draw a clean rectangular phylogram with coloured tips, support values (>=50), and scale bar."""
    terminals = list(tree.get_terminals())
    y_coords: dict[Clade, float] = {t: float(i) for i, t in enumerate(terminals)}
    x_coords: dict[Clade, float] = {}

    def _compute_x(clade: Clade, x: float) -> None:
        x_coords[clade] = x
        for ch in clade.clades:
            _compute_x(ch, x + float(ch.branch_length or 0.0))

    def _compute_y(clade: Clade) -> float:
        if clade in y_coords:
            return y_coords[clade]
        ys = [_compute_y(ch) for ch in clade.clades]
        y_coords[clade] = 0.5 * (min(ys) + max(ys))
        return y_coords[clade]

    _compute_x(tree.root, 0.0)
    _compute_y(tree.root)

    max_x = max(x_coords.values()) if x_coords else 1.0
    pad_x = max_x * 0.015

    # Draw branches
    for clade in tree.find_clades():
        x0, y0 = x_coords[clade], y_coords[clade]
        if not clade.is_terminal() and clade.clades:
            ys = [y_coords[ch] for ch in clade.clades]
            ax.plot([x0, x0], [min(ys), max(ys)], color="#333333", lw=1.3, zorder=2)
            for ch in clade.clades:
                ax.plot([x0, x_coords[ch]], [y_coords[ch], y_coords[ch]], color="#333333", lw=1.3, zorder=2)

    # Tip markers & labels
    ordered_tip_names: list[str] = []
    for t in terminals:
        ordered_tip_names.append(t.name)
        cl = sp_to_clade[t.name]
        col = TREE_CLADE_COLORS[cl]
        ax.scatter(x_coords[t], y_coords[t], color=col, edgecolor="black", lw=0.6, s=55, zorder=4)
        ax.text(
            x_coords[t] + pad_x * 1.4,
            y_coords[t],
            t.name.replace("_", " "),
            va="center",
            ha="left",
            fontsize=font_size - 0.5,
            color="#111111",
        )

    # Internal node support labels (hide < 50)
    for clade in tree.get_nonterminals():
        if clade is tree.root:
            continue
        lbl_str = ""
        if support_mode == "iqtree":
            raw = clade.name if clade.name else (str(clade.confidence) if clade.confidence is not None else None)
            sh, uf = parse_support_label(raw)
            if uf is not None and uf >= 50.0:
                lbl_str = f"{sh:.0f}/{uf:.0f}" if sh is not None else f"{uf:.0f}"
        elif support_mode == "nj":
            if clade.confidence is not None and float(clade.confidence) >= 50.0:
                lbl_str = f"{float(clade.confidence):.0f}"
        if lbl_str:
            ax.text(
                x_coords[clade] - pad_x * 0.6,
                y_coords[clade] - 0.22,
                lbl_str,
                va="bottom",
                ha="right",
                fontsize=max(7, font_size - 2.5),
                color="#005580",
                fontweight="bold",
            )

    # Scale bar
    scale_len = 0.1 if max_x >= 0.25 else 0.05
    y_bar = len(terminals) - 0.2
    ax.plot([0.0, scale_len], [y_bar, y_bar], color="black", lw=2.0)
    ax.text(scale_len / 2, y_bar + 0.35, f"{scale_len:.2f} subs/site", ha="center", va="top", fontsize=font_size - 1.5)

    ax.set_xlim(-max_x * 0.08, max_x * 1.48)
    ax.set_ylim(len(terminals) + 0.6, -0.8)
    ax.axis("off")
    ax.set_title(title, loc="left", fontweight="bold", fontsize=font_size + 1)

    if show_legend:
        present = [c for c in TREE_CLADE_ORDER if c in set(sp_to_clade.values())]
        handles = [mpatches.Patch(facecolor=TREE_CLADE_COLORS[c], edgecolor="black", lw=0.5, label=c) for c in present]
        ax.legend(
            handles=handles,
            title="Clade",
            loc="lower left",
            bbox_to_anchor=(0.0, 0.14),
            frameon=True,
            fontsize=font_size - 2,
            ncol=2,
        )

    return ordered_tip_names


def save_fig(fig: plt.Figure, stem: Path, cfg: dict[str, Any]) -> None:
    """Save figure in all configured formats and close it."""
    dpi = int(cfg["reporting"]["dpi"])
    for fmt in cfg["reporting"]["figure_formats"]:
        fig.savefig(stem.with_suffix(f".{fmt}"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_all_tree_figures(
    ml_tree: Tree,
    nj_tree: Tree,
    id_map: pd.DataFrame,
    hum_df: pd.DataFrame,
    best_model: str,
    fig_dir: Path,
    cfg: dict[str, Any],
    logger: Any,
) -> None:
    """Task 6.6: Render ML, NJ, side-by-side comparison, and ML+identity-strip tree figures."""
    fs = int(cfg["reporting"]["font_size"])
    plt.rcParams.update({"font.size": fs})
    sp_to_clade = {r.tree_id: GROUP_TO_TREE_CLADE[r.group] for r in id_map.itertuples()}

    # 1. ML rectangular
    fig, ax = plt.subplots(figsize=(10.5, 7.2))
    draw_phylogram(
        ax,
        ml_tree,
        sp_to_clade,
        "iqtree",
        f"Maximum Likelihood Phylogeny (IQ-TREE, Model: {best_model})\nNode labels: SH-aLRT / UFBoot % (values with UFBoot < 50 hidden)",
        fs,
        show_legend=True,
    )
    save_fig(fig, fig_dir / "fig_tree_ml_rectangular", cfg)

    # 2. NJ rectangular
    fig, ax = plt.subplots(figsize=(10.5, 7.2))
    draw_phylogram(
        ax,
        nj_tree,
        sp_to_clade,
        "nj",
        f"Neighbor-Joining Phylogeny (Poisson Distance, {cfg['phylogenetics']['bootstrap_replicates_nj']} Bootstrap Replicates)\nNode labels: Bootstrap % (values < 50 hidden)",
        fs,
        show_legend=True,
    )
    save_fig(fig, fig_dir / "fig_tree_nj_rectangular", cfg)

    # 3. Side-by-side NJ vs ML
    fig, axes = plt.subplots(1, 2, figsize=(17.5, 7.4))
    draw_phylogram(axes[0], nj_tree, sp_to_clade, "nj", "A. Neighbor-Joining (Poisson; Bootstrap %)", fs, show_legend=False)
    draw_phylogram(axes[1], ml_tree, sp_to_clade, "iqtree", f"B. Maximum Likelihood (IQ-TREE {best_model}; SH-aLRT/UFBoot %)", fs, show_legend=True)
    fig.suptitle("Comparison of Neighbor-Joining and Maximum Likelihood Cytochrome c Trees", y=0.99, fontsize=fs + 3, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_fig(fig, fig_dir / "fig_tree_comparison_NJ_vs_ML", cfg)

    # 4. ML tree + identity-to-human heat strip
    fig = plt.figure(figsize=(13.5, 7.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.2, 1.1], wspace=0.08)
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[0, 1])

    tip_order = draw_phylogram(
        ax_tree,
        ml_tree,
        sp_to_clade,
        "iqtree",
        f"ML Phylogeny (IQ-TREE {best_model}) & Identity to Human",
        fs,
        show_legend=True,
    )
    id_lookup = dict(zip(hum_df["species"], hum_df["identity_pct"]))
    vals = [id_lookup[sp] for sp in tip_order]
    y = np.arange(len(tip_order))
    colors = [TREE_CLADE_COLORS[sp_to_clade[sp]] for sp in tip_order]
    bars = ax_bar.barh(y, vals, color=colors, edgecolor="black", lw=0.5, height=0.68)
    for bar, v in zip(bars, vals):
        ax_bar.text(v + 1.0, bar.get_y() + bar.get_height() / 2, f"{v:.1f}%", va="center", fontsize=fs - 1.5)
    ax_bar.set_ylim(len(tip_order) + 0.6, -0.8)
    ax_bar.set_xlim(0, 118)
    ax_bar.set_yticks([])
    ax_bar.set_xlabel("Identity to Human (%)")
    ax_bar.set_title("Identity to Homo sapiens", fontsize=fs + 0.5, fontweight="bold")
    ax_bar.grid(axis="x", linestyle=":", alpha=0.4)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    save_fig(fig, fig_dir / "fig_tree_with_identity_heatmap", cfg)

    logger.info("Saved 4 tree figures (PNG + SVG) to %s (radial unrooted skipped in favor of clean rectangular + comparison layouts).", fig_dir)


def export_itol_files(
    ml_tree: Tree, id_map: pd.DataFrame, species_df: pd.DataFrame, itol_dir: Path, logger: Any
) -> None:
    """Task 6.7: Export rooted Newick (UFBoot as labels), colorstrip, common-name labels, and README."""
    itol_tree = copy.deepcopy(ml_tree)
    for clade in itol_tree.get_nonterminals():
        if clade is itol_tree.root:
            clade.name = None
            clade.confidence = None
            continue
        raw = clade.name if clade.name else (str(clade.confidence) if clade.confidence is not None else None)
        _, uf = parse_support_label(raw)
        clade.name = f"{uf:.0f}" if uf is not None else None
        clade.confidence = None
    Phylo.write(itol_tree, itol_dir / "cytc_ml_rooted.nwk", "newick")

    # Colorstrip
    cs_lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR COMMA",
        "DATASET_LABEL,Clade",
        "COLOR,#ff0000",
        "LEGEND_TITLE,Clade",
        "LEGEND_SHAPES," + ",".join(["1"] * len(TREE_CLADE_ORDER)),
        "LEGEND_COLORS," + ",".join(TREE_CLADE_COLORS[c] for c in TREE_CLADE_ORDER),
        "LEGEND_LABELS," + ",".join(TREE_CLADE_ORDER),
        "DATA",
    ]
    for r in id_map.itertuples():
        cl = GROUP_TO_TREE_CLADE[r.group]
        cs_lines.append(f"{r.tree_id},{TREE_CLADE_COLORS[cl]},{cl}")
    (itol_dir / "itol_colorstrip_clade.txt").write_text("\n".join(cs_lines) + "\n", encoding="utf-8")

    # Common-name labels
    sp_to_common = dict(zip(species_df["species"], species_df["common_name"]))
    lbl_lines = ["LABELS", "SEPARATOR TAB", "DATA"]
    for r in id_map.itertuples():
        pretty = r.tree_id.replace("_", " ")
        common = sp_to_common.get(r.tree_id, "")
        lbl_lines.append(f"{r.tree_id}\t{pretty} ({common})")
    (itol_dir / "itol_labels_common_names.txt").write_text("\n".join(lbl_lines) + "\n", encoding="utf-8")

    readme = (
        "# iTOL Export Instructions\n\n"
        "1. Open https://itol.embl.de/upload.cgi and upload `cytc_ml_rooted.nwk`.\n"
        "2. Drag and drop `itol_colorstrip_clade.txt` and `itol_labels_common_names.txt` onto the tree canvas.\n"
        "3. In Control Panel -> Advanced, enable Bootstrap / metadata display as text (values represent IQ-TREE UFBoot %).\n"
        "4. Verify rooting on *Saccharomyces cerevisiae*, then export as SVG/PDF via the Export tab.\n"
    )
    (itol_dir / "README_itol.md").write_text(readme, encoding="utf-8")
    logger.info("Exported iTOL files to %s", itol_dir)


def write_phylogeny_summary(
    out_md: Path,
    iq_summary_df: pd.DataFrame,
    clade_sup_df: pd.DataFrame,
    sens_df: pd.DataFrame,
    stable: list[str],
    unstable: list[str],
    n_clamped_nj: int,
    n_collapsed_ml: int,
    ident_check_note: str,
    cfg: dict[str, Any],
) -> None:
    """Write comprehensive Phase 6 summary markdown."""
    n_same = int(sens_df["identical_topology"].sum())
    n_sup = int(clade_sup_df["supported"].sum())
    lines = [
        "# Phase 6 Phylogenetic Reconstruction Summary (RQ4)",
        "",
        "## 1. Methods & Random Seeds",
        f"- Random seed: `{cfg['project']['random_seed']}` (plus sensitivity seeds `1` and `2` for IQ-TREE).",
        f"- NJ: Biopython Neighbor-Joining with Poisson distance ({cfg['phylogenetics']['bootstrap_replicates_nj']} bootstrap replicates; {n_clamped_nj} negative branch(es) clamped to 0) and BLOSUM62 cross-check.",
        f"- FastTree: approximate ML with `{cfg['phylogenetics']['fasttree_model_flags']}` (SH-like local supports, 0–1).",
        f"- IQ-TREE (primary): ModelFinder Plus (`{cfg['phylogenetics']['model_selection']}`), `{cfg['phylogenetics']['ufboot_replicates']}` UFBoot + `{cfg['phylogenetics']['alrt_replicates']}` SH-aLRT + `-bnni` + `-keep-ident`, rooted on `{cfg['phylogenetics']['primary_root_outgroup']}`.",
        f"- Polytomy-collapsed ML tree (`cytc_all.collapsed.nwk`): **{n_collapsed_ml}** internal branch(es) collapsed (length < 1e-6 or UFBoot < 50).",
        f"- Identical-sequence check: {ident_check_note}",
        "",
        "## 2. IQ-TREE Model & Likelihood Summary",
        df_to_markdown(iq_summary_df),
        "",
        f"## 3. Primary ML Clade Support ({n_sup}/{len(clade_sup_df)} internal nodes supported at UFBoot >= {cfg['phylogenetics']['support_threshold']} & SH-aLRT >= 80)",
        df_to_markdown(clade_sup_df),
        "",
        f"## 4. Sensitivity Analysis (Unrooted Robinson-Foulds Distance vs. Primary ML Tree; {n_same}/{len(sens_df)} identical topology)",
        df_to_markdown(sens_df),
        "",
        "## 5. Stable vs. Unstable Clades Across All 17-Taxon Variants",
        "### Stable clades (recovered in 100% of 17-taxon variants)",
        *([f"- `{s}`" for s in stable] if stable else ["- _None_"]),
        "",
        "### Unstable primary-tree clades (absent in >= 1 sensitivity variant)",
        *([f"- `{u}`" for u in unstable] if unstable else ["- _None_"]),
        "",
    ]
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run Phase 6 phylogenetic analyses, sensitivity comparisons, visualizations, and iTOL exports."""
    parser = argparse.ArgumentParser(description="Phase 6: Phylogenetic tree construction.")
    parser.add_argument("--config", type=str, default=None, help="Optional path to config.yaml")
    parser.add_argument("--reuse-trees", action="store_true", help="Reuse existing NJ/FastTree/IQ-TREE outputs if present")
    args = parser.parse_args()

    cfg = load_config(args.config)
    logger = get_logger("60_phylogenetics", cfg)
    logger.info(
        "Phase 6 goal: reconstruct NJ, FastTree, and IQ-TREE phylogenies with bootstrapping, sensitivity variants, "
        "polytomy collapsing, RF comparisons, figures, and iTOL exports. "
        "Inputs: results/alignments/primary_all.fasta, primary_vert.fasta, primary_all_trimmed.fasta, all_muscle.fasta."
    )

    aln_dir = resolve(cfg, "alignments")
    trees_dir = resolve(cfg, "trees")
    tables_dir = resolve(cfg, "tables")
    fig_dir = resolve(cfg, "figures")
    proc_dir = resolve(cfg, "processed")
    species_path = resolve(cfg, "species_table")

    nj_dir = trees_dir / "nj"
    ft_dir = trees_dir / "fasttree"
    iq_dir = trees_dir / "iqtree"
    itol_dir = trees_dir / "itol"
    for d in (nj_dir, ft_dir, iq_dir, itol_dir):
        d.mkdir(parents=True, exist_ok=True)

    id_map = pd.read_csv(proc_dir / "id_map.csv").sort_values("order_index")
    species_df = pd.read_csv(species_path)
    hum_df = pd.read_csv(tables_dir / "identity_to_human.csv")
    ident_df = pd.read_csv(tables_dir / "identical_sequence_groups.csv")

    all_taxa = list(id_map["tree_id"])
    vert_taxa = list(id_map.loc[id_map["role"] == "ingroup", "tree_id"])
    outgroups = list(cfg["phylogenetics"]["outgroups"])
    root_outgroup = str(cfg["phylogenetics"]["primary_root_outgroup"])
    support_thresh = float(cfg["phylogenetics"]["support_threshold"])

    aln_all_path = aln_dir / "primary_all.fasta"
    aln_vert_path = aln_dir / "primary_vert.fasta"
    aln_muscle_path = aln_dir / "all_muscle.fasta"
    aln_trim_path = aln_dir / "primary_all_trimmed.fasta"
    aln_all = AlignIO.read(aln_all_path, "fasta")

    # Task 6.1: NJ
    if args.reuse_trees and (nj_dir / "nj_poisson.nwk").exists() and (nj_dir / "nj_blosum62.nwk").exists():
        nj_poiss = Phylo.read(nj_dir / "nj_poisson.nwk", "newick")
        nj_blosum = Phylo.read(nj_dir / "nj_blosum62.nwk", "newick")
        n_clamped_nj = 0
        logger.info("Reusing existing NJ trees in %s", nj_dir)
    else:
        nj_poiss, nj_blosum, n_clamped_nj = run_nj_analysis(aln_all, nj_dir, cfg, logger)

    # Task 6.2: FastTree
    ft_path = ft_dir / "fasttree_lg_gamma.nwk"
    if not (args.reuse_trees and ft_path.exists()):
        ft_path = run_fasttree(aln_all_path, ft_dir, cfg, logger)

    # Task 6.3 & 6.4: IQ-TREE primary + sensitivity variants
    iq_jobs = [
        ("cytc_all", aln_all_path, iq_dir / "cytc_all", None, None, root_outgroup),
        ("cytc_vert", aln_vert_path, iq_dir / "cytc_vert", None, None, None),
        ("cytc_all_lgg4", aln_all_path, iq_dir / "cytc_all_lgg4", "LG+G4", None, root_outgroup),
        ("cytc_all_muscle", aln_muscle_path, iq_dir / "cytc_all_muscle", None, None, root_outgroup),
        ("cytc_all_trimmed", aln_trim_path, iq_dir / "cytc_all_trimmed", None, None, root_outgroup),
        ("cytc_all_seed1", aln_all_path, iq_dir / "cytc_all_seed1", None, 1, root_outgroup),
        ("cytc_all_seed2", aln_all_path, iq_dir / "cytc_all_seed2", None, 2, root_outgroup),
    ]
    iq_reports: list[dict[str, Any]] = []
    iq_treefiles: dict[str, Path] = {}
    for run_name, aln_p, pref, mod, sd, og in iq_jobs:
        tf = Path(f"{pref}.treefile")
        if not (args.reuse_trees and tf.exists() and Path(f"{pref}.iqtree").exists()):
            tf = run_iqtree_job(aln_p, pref, cfg, logger, model=mod, seed=sd, outgroup=og)
        iq_treefiles[run_name] = tf
        iq_reports.append(parse_iqtree_report(Path(f"{pref}.iqtree"), run_name, aln_p.name))

    iq_summary_df = pd.DataFrame(iq_reports)
    iq_summary_df.to_csv(tables_dir / "iqtree_summary.csv", index=False)

    # Primary ML clade support
    ml_tree, clade_sup_df = extract_ml_clade_support(iq_treefiles["cytc_all"], root_outgroup, support_thresh)
    clade_sup_df.to_csv(tables_dir / "ml_clade_support.csv", index=False)

    # Task 6.5: Identical-sequence check & polytomy collapse
    ident_notes: list[str] = []
    for r in ident_df.itertuples():
        members = str(r.members).split(";")
        if len(members) == 2:
            d = ml_tree.distance(members[0], members[1])
            ident_notes.append(f"{members[0]} & {members[1]} patristic distance = {d:.6f} (resolution arbitrary)")
    ident_check_note = "; ".join(ident_notes) if ident_notes else "No identical-sequence groups."
    logger.info("Identical-sequence check: %s", ident_check_note)

    collapsed_path = iq_dir / "cytc_all.collapsed.nwk"
    n_collapsed_ml = collapse_unsupported_branches(iq_treefiles["cytc_all"], collapsed_path, root_outgroup)
    logger.info("Polytomy-collapsed ML tree written (%d branches collapsed): %s", n_collapsed_ml, collapsed_path)

    # Task 6.4: Sensitivity RF table
    variants = [
        ("nj_poisson", nj_dir / "nj_poisson.nwk", False),
        ("nj_blosum62", nj_dir / "nj_blosum62.nwk", False),
        ("fasttree_lg_gamma", ft_path, False),
        ("iqtree_vert_only", iq_treefiles["cytc_vert"], True),
        ("iqtree_fixed_lgg4", iq_treefiles["cytc_all_lgg4"], False),
        ("iqtree_muscle", iq_treefiles["cytc_all_muscle"], False),
        ("iqtree_trimmed", iq_treefiles["cytc_all_trimmed"], False),
        ("iqtree_seed1", iq_treefiles["cytc_all_seed1"], False),
        ("iqtree_seed2", iq_treefiles["cytc_all_seed2"], False),
    ]
    sens_df, stable, unstable = compute_sensitivity_rf(
        iq_treefiles["cytc_all"], variants, all_taxa, vert_taxa, outgroups
    )
    sens_path = trees_dir / "sensitivity_summary.csv"
    sens_df.to_csv(sens_path, index=False)

    # Task 6.6: Tree figures
    best_model = str(iq_summary_df.loc[iq_summary_df["run"] == "cytc_all", "best_model"].iloc[0])
    plot_all_tree_figures(ml_tree, nj_poiss, id_map, hum_df, best_model, fig_dir, cfg, logger)

    # Task 6.7: iTOL exports
    export_itol_files(ml_tree, id_map, species_df, itol_dir, logger)

    # Summary markdown
    write_phylogeny_summary(
        trees_dir / "phylogeny_summary.md",
        iq_summary_df,
        clade_sup_df,
        sens_df,
        stable,
        unstable,
        n_clamped_nj,
        n_collapsed_ml,
        ident_check_note,
        cfg,
    )

    logger.info(
        "SUMMARY: ML best_model=%s (supported nodes=%d/%d, collapsed=%d); NJ clamped=%d; "
        "sensitivity variants=%d (identical=%d); outputs=%s, %s",
        best_model,
        int(clade_sup_df["supported"].sum()),
        len(clade_sup_df),
        n_collapsed_ml,
        n_clamped_nj,
        len(sens_df),
        int(sens_df["identical_topology"].sum()),
        iq_treefiles["cytc_all"],
        sens_path,
    )


if __name__ == "__main__":
    main()
