"""Phase 5: Pairwise sequence identity, similarity, evolutionary distance, clock test, and heatmaps."""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Align, AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment, substitution_matrices
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from scipy.stats import linregress, pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cytc.utils import REPO_ROOT, get_logger, load_config, resolve  # noqa: E402

GROUP_TO_CLADE: dict[str, str] = {
    "Primate": "Mammal",
    "Rodent": "Mammal",
    "Carnivore": "Mammal",
    "Artiodactyl": "Mammal",
    "Perissodactyl": "Mammal",
    "Marsupial": "Mammal",
    "Bird": "Bird",
    "Reptile": "Reptile",
    "Amphibian": "Amphibian",
    "Fish": "Fish",
    "Insect": "Insect",
    "Fungus": "Fungus",
}

CLADE_ORDER = ["Mammal", "Bird", "Reptile", "Amphibian", "Fish", "Insect", "Fungus"]
CLADE_COLORS = dict(zip(CLADE_ORDER, sns.color_palette("colorblind", n_colors=len(CLADE_ORDER))))

# Approximate MRCA divergence times (Mya) consistent with config/species.csv and TimeTree consensus
# (Labelled approximate; pairs are non-independent due to shared phylogeny)
RELATIONSHIP_ORDER = [
    "within_primates",
    "within_rodents",
    "primate_vs_rodent",
    "within_other_mammals",
    "within_birds",
    "within_reptiles",
    "within_fish",
    "bird_vs_reptile",
    "mammal_vs_bird",
    "mammal_vs_reptile",
    "mammal_vs_amphibian",
    "bird_vs_amphibian",
    "reptile_vs_amphibian",
    "mammal_vs_fish",
    "bird_vs_fish",
    "reptile_vs_fish",
    "amphibian_vs_fish",
    "vertebrate_vs_drosophila",
    "vertebrate_vs_yeast",
    "drosophila_vs_yeast",
]


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


def classify_relationship(group1: str, group2: str) -> tuple[str, bool]:
    """Return (relationship_label, both_mammals) for a pair of taxonomic groups."""
    c1, c2 = GROUP_TO_CLADE[group1], GROUP_TO_CLADE[group2]
    both_mammals = c1 == "Mammal" and c2 == "Mammal"

    if c1 == "Insect" and c2 == "Fungus" or c1 == "Fungus" and c2 == "Insect":
        return "drosophila_vs_yeast", False
    if "Fungus" in (c1, c2):
        return "vertebrate_vs_yeast", False
    if "Insect" in (c1, c2):
        return "vertebrate_vs_drosophila", False

    if both_mammals:
        if group1 == "Primate" and group2 == "Primate":
            return "within_primates", True
        if group1 == "Rodent" and group2 == "Rodent":
            return "within_rodents", True
        if {group1, group2} == {"Primate", "Rodent"}:
            return "primate_vs_rodent", True
        return "within_other_mammals", True

    if c1 == c2:
        if c1 == "Bird":
            return "within_birds", False
        if c1 == "Reptile":
            return "within_reptiles", False
        if c1 == "Fish":
            return "within_fish", False
        return f"within_{c1.lower()}s", False

    pair_sorted = sorted([c1, c2], key=lambda c: CLADE_ORDER.index(c))
    return f"{pair_sorted[0].lower()}_vs_{pair_sorted[1].lower()}", False


def estimate_pair_mya(sp1: str, sp2: str, group1: str, group2: str) -> float:
    """Return approximate MRCA divergence time (Mya) consistent with species.csv and TimeTree."""
    if sp1 == sp2:
        return 0.0
    c1, c2 = GROUP_TO_CLADE[group1], GROUP_TO_CLADE[group2]
    clades = {c1, c2}
    groups = {group1, group2}

    if "Fungus" in clades:
        return 1100.0
    if "Insect" in clades:
        return 780.0
    if "Fish" in clades:
        return 230.0 if clades == {"Fish"} else 430.0
    if "Amphibian" in clades:
        return 352.0
    if "Mammal" in clades and clades != {"Mammal"}:
        return 320.0
    if clades <= {"Bird", "Reptile"}:
        if clades == {"Bird"}:
            return 90.0
        # Turtle (Chrysemys_picta) is Archelosauria (sister to birds ~255 Mya); Anole is squamate (~280 Mya)
        if "Anolis_carolinensis" in (sp1, sp2):
            return 280.0
        return 255.0
    # Both mammals
    if "Marsupial" in groups:
        return 160.0
    if groups == {"Primate"}:
        return 29.0
    if groups == {"Rodent"}:
        return 20.0
    if groups == {"Primate", "Rodent"}:
        return 90.0
    laurasiatheria = {"Carnivore", "Artiodactyl", "Perissodactyl"}
    if groups <= laurasiatheria:
        return 80.0
    return 94.0


def compute_msa_pairwise_metrics(
    aln: MultipleSeqAlignment,
    id_map: pd.DataFrame,
    cfg: dict[str, Any],
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Compute all 6 square symmetric matrices and tidy pairwise_long DataFrame from primary_all.fasta."""
    blosum = substitution_matrices.load(cfg["similarity"]["substitution_matrix"])
    ids = [r.id for r in aln]
    seqs = [str(r.seq) for r in aln]
    ungapped_lens = [len(s.replace("-", "")) for s in seqs]
    n = len(ids)

    sp_to_group = dict(zip(id_map["tree_id"], id_map["group"]))

    mat_id = np.full((n, n), 100.0)
    mat_sim = np.full((n, n), 100.0)
    mat_subs = np.zeros((n, n), dtype=int)
    mat_p = np.zeros((n, n), dtype=float)
    mat_poiss = np.zeros((n, n), dtype=float)
    mat_blosum = np.zeros((n, n), dtype=float)

    for i in range(n):
        mat_blosum[i, i] = sum(float(blosum[a, a]) for a in seqs[i] if a != "-")

    long_rows: list[dict[str, Any]] = []
    for i in range(n):
        si = seqs[i]
        for j in range(i + 1, n):
            sj = seqs[j]
            n_comp = 0
            n_ident = 0
            n_pos = 0
            b_score = 0.0
            for a, b in zip(si, sj):
                if a == "-" or b == "-":
                    continue
                n_comp += 1
                sc = float(blosum[a, b])
                b_score += sc
                if a == b:
                    n_ident += 1
                if sc > 0:
                    n_pos += 1

            n_sub = n_comp - n_ident
            id_pct = 100.0 * n_ident / n_comp if n_comp > 0 else 0.0
            id_overall = 100.0 * n_ident / min(ungapped_lens[i], ungapped_lens[j])
            sim_pct = 100.0 * n_pos / n_comp if n_comp > 0 else 0.0
            p_dist = 1.0 - (id_pct / 100.0)
            if p_dist <= 0.0:
                poiss = 0.0
            elif p_dist >= 1.0:
                poiss = float("inf")
            else:
                poiss = -math.log(1.0 - p_dist)

            mat_id[i, j] = mat_id[j, i] = round(id_pct, 4)
            mat_sim[i, j] = mat_sim[j, i] = round(sim_pct, 4)
            mat_subs[i, j] = mat_subs[j, i] = n_sub
            mat_p[i, j] = mat_p[j, i] = round(p_dist, 6)
            mat_poiss[i, j] = mat_poiss[j, i] = round(poiss, 6)
            mat_blosum[i, j] = mat_blosum[j, i] = round(b_score, 1)

            g1, g2 = sp_to_group[ids[i]], sp_to_group[ids[j]]
            rel, both_mam = classify_relationship(g1, g2)
            mya = estimate_pair_mya(ids[i], ids[j], g1, g2)
            long_rows.append(
                {
                    "sp1": ids[i],
                    "sp2": ids[j],
                    "group1": g1,
                    "group2": g2,
                    "clade1": GROUP_TO_CLADE[g1],
                    "clade2": GROUP_TO_CLADE[g2],
                    "both_mammals": both_mam,
                    "relationship": rel,
                    "identity_pct": round(id_pct, 4),
                    "identity_pct_overall": round(id_overall, 4),
                    "similarity_pct": round(sim_pct, 4),
                    "n_substitutions": n_sub,
                    "p_distance": round(p_dist, 6),
                    "poisson_distance": round(poiss, 6),
                    "blosum62_score": round(b_score, 1),
                    "n_compared": n_comp,
                    "mya_divergence": mya,
                }
            )

    matrices = {
        "identity_matrix": pd.DataFrame(mat_id, index=ids, columns=ids),
        "similarity_matrix": pd.DataFrame(mat_sim, index=ids, columns=ids),
        "n_substitutions_matrix": pd.DataFrame(mat_subs, index=ids, columns=ids),
        "p_distance_matrix": pd.DataFrame(mat_p, index=ids, columns=ids),
        "poisson_distance_matrix": pd.DataFrame(mat_poiss, index=ids, columns=ids),
        "blosum62_score_matrix": pd.DataFrame(mat_blosum, index=ids, columns=ids),
    }
    return matrices, pd.DataFrame(long_rows)


def parse_clustalo_distmat(path: Path) -> pd.DataFrame:
    """Parse Clustal Omega --percent-id --full distance matrix into a DataFrame."""
    lines = [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    n = int(lines[0])
    ids: list[str] = []
    vals: list[list[float]] = []
    for ln in lines[1 : n + 1]:
        parts = ln.split()
        ids.append(parts[0])
        vals.append([float(x) for x in parts[1 : n + 1]])
    return pd.DataFrame(vals, index=ids, columns=ids)


def run_cross_validation(
    curated_fasta: Path,
    long_df: pd.DataFrame,
    id_mat: pd.DataFrame,
    distmat_path: Path,
    cfg: dict[str, Any],
    out_md: Path,
    logger: Any,
) -> tuple[float, float]:
    """Cross-check MSA identities against global PairwiseAligner and Clustal Omega distmat."""
    recs = {r.id: str(r.seq) for r in SeqIO.parse(curated_fasta, "fasta")}
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load(cfg["similarity"]["substitution_matrix"])
    aligner.open_gap_score = float(cfg["similarity"]["gap_open"])
    aligner.extend_gap_score = float(cfg["similarity"]["gap_extend"])

    pw_ids: list[float] = []
    msa_ids: list[float] = []
    diffs_gt2: list[dict[str, Any]] = []

    for r in long_df.itertuples():
        aln_pw = aligner.align(recs[r.sp1], recs[r.sp2])[0]
        # Reconstruct aligned strings from coordinate blocks
        s1_blocks, s2_blocks = aln_pw.aligned
        n_comp = 0
        n_ident = 0
        for (a0, a1), (b0, b1) in zip(s1_blocks, s2_blocks):
            seg1 = recs[r.sp1][a0:a1]
            seg2 = recs[r.sp2][b0:b1]
            n_comp += len(seg1)
            n_ident += sum(1 for x, y in zip(seg1, seg2) if x == y)
        pw_pct = 100.0 * n_ident / n_comp if n_comp > 0 else 0.0
        pw_ids.append(pw_pct)
        msa_ids.append(float(r.identity_pct))
        diff = abs(float(r.identity_pct) - pw_pct)
        if diff > 2.0:
            diffs_gt2.append(
                {
                    "sp1": r.sp1,
                    "sp2": r.sp2,
                    "msa_identity_pct": round(float(r.identity_pct), 2),
                    "pairwise_identity_pct": round(pw_pct, 2),
                    "abs_diff_pct": round(diff, 2),
                }
            )

    r_pw, _ = pearsonr(msa_ids, pw_ids)
    rho_pw, _ = spearmanr(msa_ids, pw_ids)
    max_diff_pw = float(np.max(np.abs(np.array(msa_ids) - np.array(pw_ids))))
    mean_diff_pw = float(np.mean(np.abs(np.array(msa_ids) - np.array(pw_ids))))

    # Third check: Clustal Omega distmat
    co_note = "_Clustal Omega distmat not found._"
    if distmat_path.exists():
        co_df = parse_clustalo_distmat(distmat_path)
        co_aligned = co_df.loc[id_mat.index, id_mat.columns]
        tri = np.triu_indices(len(id_mat), k=1)
        co_vals = co_aligned.values[tri]
        msa_vals = id_mat.values[tri]
        r_co, _ = pearsonr(msa_vals, co_vals)
        max_diff_co = float(np.max(np.abs(msa_vals - co_vals)))
        co_note = (
            f"- **Clustal Omega `--percent-id` distmat (`all_clustalo.distmat`) vs. MSA identity matrix:** "
            f"Pearson r = **{r_co:.6f}**, max absolute difference = **{max_diff_co:.4f}** percentage points."
        )

    lines = [
        "# Phase 5 Cross-Validation of Pairwise Sequence Identities",
        "",
        "## 1. MSA-Derived Identity vs. Independent Global Pairwise Alignment (`Bio.Align.PairwiseAligner`)",
        f"- Settings: `mode=global`, `substitution_matrix={cfg['similarity']['substitution_matrix']}`, "
        f"`open_gap_score={cfg['similarity']['gap_open']}`, `extend_gap_score={cfg['similarity']['gap_extend']}`.",
        f"- Pairs compared: **{len(msa_ids)}**",
        f"- Pearson correlation: **r = {r_pw:.6f}**",
        f"- Spearman correlation: **rho = {rho_pw:.6f}**",
        f"- Mean absolute difference: **{mean_diff_pw:.4f}** percentage points",
        f"- Maximum absolute difference: **{max_diff_pw:.4f}** percentage points",
        "",
        "### Pairs differing by > 2.0 percentage points",
        df_to_markdown(pd.DataFrame(diffs_gt2)) if diffs_gt2 else "_None (all 136 pairs agree within <= 2.0 percentage points)._",
        "",
        "## 2. Third Cross-Check: Clustal Omega Distance Matrix",
        co_note,
        "",
    ]
    out_md.write_text("\n".join(lines), encoding="utf-8")
    logger.info(
        "Cross-check OK: MSA vs PairwiseAligner r=%.6f, max_diff=%.4f%% (pairs >2%%: %d)",
        r_pw,
        max_diff_pw,
        len(diffs_gt2),
    )
    return r_pw, max_diff_pw


def write_metric_definitions(out_md: Path, cfg: dict[str, Any]) -> None:
    """Write precise definitions of all pairwise metrics to metric_definitions.md."""
    text = f"""# Phase 5 Pairwise Similarity Metric Definitions

**Headline metrics:** `identity_pct` (gap-free pairwise aligned columns) and `similarity_pct` (BLOSUM62 positive substitutions over gap-free pairwise aligned columns).

| Metric | Definition |
| --- | --- |
| `n_compared` | Number of alignment columns where **both** sequences contain a standard amino acid (gap-gap and residue-gap columns excluded). |
| `identity_pct` | `100 * identical / n_compared` (headline percent identity). |
| `identity_pct_overall` | `100 * identical / min(len_i, len_j)` using ungapped sequence lengths (penalizes terminal extensions/indels). |
| `similarity_pct` | `100 * (compared columns with {cfg['similarity']['substitution_matrix']} score > 0) / n_compared` (includes identical residues and biochemically conservative substitutions). |
| `p_distance` | `1 - identity_pct / 100` (observed proportion of differing sites). |
| `poisson_distance` | `-ln(1 - p_distance)` (Poisson correction for multiple substitutions; `0` when `p = 0`; guarded if `p >= 1`). |
| `n_substitutions` | Number of non-identical compared columns (`n_compared - identical`). |
| `blosum62_score` | Sum of `{cfg['similarity']['substitution_matrix']}` substitution scores across all `n_compared` columns. |
"""
    out_md.write_text(text, encoding="utf-8")


def build_group_summary(long_df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-relationship summary statistics (plus mammal_vs_mammal_all superset)."""
    rows: list[dict[str, Any]] = []
    present_rels = [r for r in RELATIONSHIP_ORDER if r in set(long_df["relationship"])]
    subsets: list[tuple[str, pd.DataFrame]] = [
        (rel, long_df[long_df["relationship"] == rel]) for rel in present_rels
    ]
    subsets.append(("mammal_vs_mammal_all", long_df[long_df["both_mammals"]]))

    for label, sub in subsets:
        idx_min = sub["identity_pct"].idxmin()
        idx_max = sub["identity_pct"].idxmax()
        r_min = sub.loc[idx_min]
        r_max = sub.loc[idx_max]
        rows.append(
            {
                "relationship": label,
                "n_pairs": int(len(sub)),
                "mean_identity": round(float(sub["identity_pct"].mean()), 2),
                "median_identity": round(float(sub["identity_pct"].median()), 2),
                "min_identity": round(float(sub["identity_pct"].min()), 2),
                "max_identity": round(float(sub["identity_pct"].max()), 2),
                "sd_identity": round(float(sub["identity_pct"].std(ddof=1)), 2) if len(sub) > 1 else 0.0,
                "min_pair": f"{r_min['sp1']} vs {r_min['sp2']}",
                "max_pair": f"{r_max['sp1']} vs {r_max['sp2']}",
            }
        )
    return pd.DataFrame(rows)


def build_identity_to_human(
    matrices: dict[str, pd.DataFrame], id_map: pd.DataFrame
) -> pd.DataFrame:
    """Build the human-centric identity/similarity/distance table sorted by identity descending."""
    hum = "Homo_sapiens"
    rows: list[dict[str, Any]] = []
    for r in id_map.itertuples():
        sp = r.tree_id
        grp = r.group
        rows.append(
            {
                "species": sp,
                "group": grp,
                "clade": GROUP_TO_CLADE[grp],
                "identity_pct": round(float(matrices["identity_matrix"].loc[hum, sp]), 2),
                "similarity_pct": round(float(matrices["similarity_matrix"].loc[hum, sp]), 2),
                "n_substitutions": int(matrices["n_substitutions_matrix"].loc[hum, sp]),
                "poisson_distance": round(float(matrices["poisson_distance_matrix"].loc[hum, sp]), 4),
                "mya_divergence": estimate_pair_mya(hum, sp, "Primate", grp),
            }
        )
    df = pd.DataFrame(rows).sort_values(["identity_pct", "mya_divergence"], ascending=[False, True]).reset_index(drop=True)
    return df


def save_fig(fig: plt.Figure, stem: Path, cfg: dict[str, Any]) -> None:
    """Save matplotlib figure in all configured formats and close it."""
    dpi = int(cfg["reporting"]["dpi"])
    for fmt in cfg["reporting"]["figure_formats"]:
        fig.savefig(stem.with_suffix(f".{fmt}"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_annotated_heatmap(
    mat_df: pd.DataFrame,
    id_map: pd.DataFrame,
    title: str,
    cbar_label: str,
    stem: Path,
    cfg: dict[str, Any],
) -> None:
    """Plot square annotated heatmap in fixed biological order with a coloured clade side-bar."""
    fs = int(cfg["reporting"]["font_size"])
    plt.rcParams.update({"font.size": fs})
    labels = [s.replace("_", " ") for s in mat_df.index]
    sp_to_clade = {r.tree_id: GROUP_TO_CLADE[r.group] for r in id_map.itertuples()}
    clades = [sp_to_clade[s] for s in mat_df.index]

    fig = plt.figure(figsize=(12.2, 10.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.35, 10.0, 0.45], wspace=0.04)
    ax_bar = fig.add_subplot(gs[0, 0])
    ax_hm = fig.add_subplot(gs[0, 1])
    ax_cb = fig.add_subplot(gs[0, 2])

    # Clade side-bar
    for idx, cl in enumerate(clades):
        ax_bar.add_patch(mpatches.Rectangle((0, idx), 1, 1, facecolor=CLADE_COLORS[cl], edgecolor="white", lw=0.5))
    ax_bar.set_xlim(0, 1)
    ax_bar.set_ylim(len(clades), 0)
    ax_bar.set_xticks([])
    ax_bar.set_yticks(np.arange(len(labels)) + 0.5)
    ax_bar.set_yticklabels(labels, fontsize=fs - 1)
    ax_bar.set_ylabel("Species (Biological Order)")
    for spn in ax_bar.spines.values():
        spn.set_visible(False)

    sns.heatmap(
        mat_df,
        annot=True,
        fmt=".1f",
        cmap="cividis",
        ax=ax_hm,
        cbar_ax=ax_cb,
        cbar_kws={"label": cbar_label},
        linewidths=0.3,
        annot_kws={"size": 6.5},
        xticklabels=labels,
        yticklabels=False,
    )
    ax_hm.set_xticklabels(labels, rotation=45, ha="right", fontsize=fs - 1)
    ax_hm.set_title(title, pad=12, fontweight="bold")

    present_clades = [c for c in CLADE_ORDER if c in set(clades)]
    handles = [mpatches.Patch(facecolor=CLADE_COLORS[c], edgecolor="black", lw=0.5, label=c) for c in present_clades]
    ax_hm.legend(
        handles=handles,
        title="Clade",
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        ncol=len(present_clades),
        frameon=True,
        fontsize=fs - 1.5,
    )
    save_fig(fig, stem, cfg)


def plot_clustermap(
    id_mat: pd.DataFrame,
    poiss_mat: pd.DataFrame,
    id_map: pd.DataFrame,
    fig_dir: Path,
    cfg: dict[str, Any],
) -> None:
    """Plot seaborn clustermap (average linkage on Poisson distance) with clade side-colours."""
    fs = int(cfg["reporting"]["font_size"])
    plt.rcParams.update({"font.size": fs})
    sp_to_clade = {r.tree_id: GROUP_TO_CLADE[r.group] for r in id_map.itertuples()}
    clade_series = pd.Series([CLADE_COLORS[sp_to_clade[s]] for s in id_mat.index], index=id_mat.index, name="Clade")

    condensed = squareform(poiss_mat.values, checks=False)
    Z = linkage(condensed, method="average")

    g = sns.clustermap(
        id_mat,
        row_linkage=Z,
        col_linkage=Z,
        row_colors=clade_series,
        col_colors=clade_series,
        cmap="cividis",
        annot=True,
        fmt=".1f",
        annot_kws={"size": 6.5},
        linewidths=0.3,
        figsize=(11.5, 10.5),
        cbar_kws={"label": "Pairwise Identity (%)"},
    )
    g.ax_heatmap.set_xticklabels([t.get_text().replace("_", " ") for t in g.ax_heatmap.get_xticklabels()], rotation=45, ha="right", fontsize=fs - 1)
    g.ax_heatmap.set_yticklabels([t.get_text().replace("_", " ") for t in g.ax_heatmap.get_yticklabels()], rotation=0, fontsize=fs - 1)
    g.fig.suptitle("Clustered Pairwise Identity Heatmap (Average Linkage on Poisson Distance)", y=1.01, fontsize=fs + 2, fontweight="bold")

    present_clades = [c for c in CLADE_ORDER if c in set(sp_to_clade.values())]
    handles = [mpatches.Patch(facecolor=CLADE_COLORS[c], edgecolor="black", lw=0.5, label=c) for c in present_clades]
    g.fig.legend(
        handles=handles,
        title="Clade",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.04),
        ncol=len(present_clades),
        frameon=True,
        fontsize=fs - 1.5,
    )

    dpi = int(cfg["reporting"]["dpi"])
    for fmt in cfg["reporting"]["figure_formats"]:
        g.fig.savefig((fig_dir / "fig_identity_clustermap").with_suffix(f".{fmt}"), dpi=dpi, bbox_inches="tight")
    plt.close(g.fig)


def plot_identity_to_human(hum_df: pd.DataFrame, fig_dir: Path, cfg: dict[str, Any]) -> None:
    """Horizontal bar chart of pairwise identity to Homo sapiens, annotated with Mya."""
    fs = int(cfg["reporting"]["font_size"])
    sub = hum_df[hum_df["species"] != "Homo_sapiens"].copy().reset_index(drop=True)
    y = np.arange(len(sub))
    colors = [CLADE_COLORS[c] for c in sub["clade"]]

    fig, ax = plt.subplots(figsize=(10.5, 6.4))
    bars = ax.barh(y, sub["identity_pct"], color=colors, edgecolor="black", linewidth=0.6)
    for bar, r in zip(bars, sub.itertuples()):
        ax.text(
            r.identity_pct + 0.6,
            bar.get_y() + bar.get_height() / 2,
            f"{r.identity_pct:.1f}% ({int(r.mya_divergence)} Mya, {r.n_substitutions} diff)",
            va="center",
            fontsize=fs - 1.5,
        )

    ax.set_yticks(y)
    ax.set_yticklabels([s.replace("_", " ") for s in sub["species"]])
    ax.invert_yaxis()
    ax.set_xlim(0, 128)
    ax.set_xlabel("Pairwise amino-acid identity to Homo sapiens (%)")
    ax.set_title("Cytochrome c Sequence Identity to Human Ordered by Similarity")
    ax.grid(axis="x", linestyle=":", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    present_clades = [c for c in CLADE_ORDER if c in set(sub["clade"])]
    handles = [mpatches.Patch(facecolor=CLADE_COLORS[c], edgecolor="black", lw=0.5, label=c) for c in present_clades]
    ax.legend(handles=handles, title="Clade", loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=True)
    save_fig(fig, fig_dir / "fig_identity_to_human", cfg)


def plot_identity_vs_time(
    long_df: pd.DataFrame,
    hum_df: pd.DataFrame,
    clock_stats: dict[str, Any],
    fig_dir: Path,
    cfg: dict[str, Any],
) -> None:
    """Two-panel scatter: Identity % and Poisson distance vs divergence time (Mya) with clock regression."""
    fs = int(cfg["reporting"]["font_size"])
    hum_sub = hum_df[hum_df["species"] != "Homo_sapiens"].copy()

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))

    # Panel A: Identity vs Mya
    ax = axes[0]
    ax.scatter(long_df["mya_divergence"], long_df["identity_pct"], color="#BBBBBB", alpha=0.55, s=28, label="All 136 pairs")
    for cl in CLADE_ORDER:
        csub = hum_sub[hum_sub["clade"] == cl]
        if not csub.empty:
            ax.scatter(csub["mya_divergence"], csub["identity_pct"], color=CLADE_COLORS[cl], edgecolor="black", lw=0.6, s=65, label=f"Human vs {cl}", zorder=4)
    ax.set_xlabel("Approximate divergence time (Mya)")
    ax.set_ylabel("Pairwise sequence identity (%)")
    ax.set_title(
        f"A. Sequence Identity vs. Divergence Time\n"
        f"(All pairs Spearman rho = {clock_stats['all_id_rho']:.3f}; Human-vs-X rho = {clock_stats['hum_id_rho']:.3f})"
    )
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel B: Poisson distance vs Mya (molecular clock)
    ax = axes[1]
    ax.scatter(long_df["mya_divergence"], long_df["poisson_distance"], color="#BBBBBB", alpha=0.55, s=28, label="All 136 pairs")
    for cl in CLADE_ORDER:
        csub = hum_sub[hum_sub["clade"] == cl]
        if not csub.empty:
            ax.scatter(csub["mya_divergence"], csub["poisson_distance"], color=CLADE_COLORS[cl], edgecolor="black", lw=0.6, s=65, zorder=4)

    x_line = np.linspace(0, 1150, 100)
    ax.plot(
        x_line,
        clock_stats["hum_intercept"] + clock_stats["hum_slope"] * x_line,
        color="#D55E00",
        linestyle="--",
        linewidth=1.8,
        label=f"Human-vs-X OLS (R²={clock_stats['hum_r2']:.3f}, slope={clock_stats['hum_slope']*1000:.2f}/Gyr)",
    )

    # Label top 3 outliers in Panel B
    for idx, r in enumerate(clock_stats["top_outliers"].itertuples()):
        dy = 0.028 if r.residual > 0 else -0.036
        dx = -25 if r.mya_divergence > 700 else 20
        ha = "right" if r.mya_divergence > 700 else "left"
        ax.annotate(
            f"{r.species.replace('_', ' ')} ({r.residual:+.3f})",
            xy=(r.mya_divergence, r.poisson_distance),
            xytext=(r.mya_divergence + dx, r.poisson_distance + dy),
            fontsize=fs - 1.5,
            ha=ha,
            arrowprops={"arrowstyle": "->", "color": "black", "lw": 0.7},
        )

    ax.set_xlabel("Approximate divergence time (Mya)")
    ax.set_ylabel("Poisson evolutionary distance (-ln(1 - p))")
    ax.set_title(
        f"B. Molecular-Clock Test (Poisson Distance vs. Time)\n"
        f"(All pairs R² = {clock_stats['all_r2']:.3f}; Human-vs-X R² = {clock_stats['hum_r2']:.3f})"
    )
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper left", frameon=True, fontsize=fs - 1.5)
    axes[0].legend(loc="lower left", frameon=True, fontsize=fs - 2)

    fig.tight_layout()
    save_fig(fig, fig_dir / "fig_identity_vs_time", cfg)


def plot_distribution_by_relationship(long_df: pd.DataFrame, fig_dir: Path, cfg: dict[str, Any]) -> None:
    """Box + strip plot of pairwise identity grouped by taxonomic relationship (close to distant)."""
    fs = int(cfg["reporting"]["font_size"])
    present = [r for r in RELATIONSHIP_ORDER if r in set(long_df["relationship"])]

    fig, ax = plt.subplots(figsize=(12.5, 5.8))
    sns.boxplot(
        data=long_df,
        x="relationship",
        y="identity_pct",
        order=present,
        color="#E8F4FA",
        linewidth=1.0,
        fliersize=0,
        ax=ax,
    )
    sns.stripplot(
        data=long_df,
        x="relationship",
        y="identity_pct",
        order=present,
        color="#0072B2",
        size=5,
        alpha=0.8,
        jitter=0.18,
        ax=ax,
    )
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels([r.replace("_", " ") for r in present], rotation=40, ha="right", fontsize=fs - 1.5)
    ax.set_xlabel("Pairwise Taxonomic Relationship (Ordered Close -> Distant)")
    ax.set_ylabel("Pairwise sequence identity (%)")
    ax.set_title("Distribution of Cytochrome c Pairwise Identity by Evolutionary Relationship")
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save_fig(fig, fig_dir / "fig_identity_distribution_by_relationship", cfg)


def write_summary_md(
    out_md: Path,
    long_df: pd.DataFrame,
    group_sum: pd.DataFrame,
    hum_df: pd.DataFrame,
    clock_stats: dict[str, Any],
) -> None:
    """Write RQ3 headline findings, key pair comparisons, and molecular-clock statistics."""
    ident_pairs = long_df[long_df["identity_pct"] >= 100.0 - 1e-6]
    most_sim = long_df.sort_values("identity_pct", ascending=False).iloc[0]
    vert_pairs = long_df[(long_df["clade1"].isin(["Mammal", "Bird", "Reptile", "Amphibian", "Fish"])) &
                         (long_df["clade2"].isin(["Mammal", "Bird", "Reptile", "Amphibian", "Fish"]))]
    least_vert = vert_pairs.sort_values("identity_pct", ascending=True).iloc[0]

    def _pair_lookup(a: str, b: str) -> pd.Series:
        m = long_df[((long_df["sp1"] == a) & (long_df["sp2"] == b)) | ((long_df["sp1"] == b) & (long_df["sp2"] == a))]
        return m.iloc[0]

    key_pairs = [
        ("Human vs Macaque (close primate)", "Homo_sapiens", "Macaca_mulatta"),
        ("Mouse vs Rat (close rodent)", "Mus_musculus", "Rattus_norvegicus"),
        ("Chicken vs Zebra finch (close bird)", "Gallus_gallus", "Taeniopygia_guttata"),
        ("Zebrafish vs Fugu (teleost fish)", "Danio_rerio", "Takifugu_rubripes"),
        ("Human vs Zebrafish (distant vert)", "Homo_sapiens", "Danio_rerio"),
        ("Human vs Fugu (distant vert)", "Homo_sapiens", "Takifugu_rubripes"),
        ("Human vs Yeast (outgroup)", "Homo_sapiens", "Saccharomyces_cerevisiae"),
    ]
    kp_rows = []
    for label, a, b in key_pairs:
        r = _pair_lookup(a, b)
        kp_rows.append(
            {
                "comparison": label,
                "identity_pct": round(float(r["identity_pct"]), 2),
                "similarity_pct": round(float(r["similarity_pct"]), 2),
                "n_substitutions": int(r["n_substitutions"]),
                "poisson_distance": round(float(r["poisson_distance"]), 4),
                "mya": int(r["mya_divergence"]),
            }
        )

    mam_mean = float(long_df.loc[long_df["both_mammals"], "identity_pct"].mean())
    mb_mean = float(long_df.loc[long_df["relationship"] == "mammal_vs_bird", "identity_pct"].mean())
    mf_mean = float(long_df.loc[long_df["relationship"] == "mammal_vs_fish", "identity_pct"].mean())
    vy_mean = float(long_df.loc[long_df["relationship"] == "vertebrate_vs_yeast", "identity_pct"].mean())

    lines = [
        "# Phase 5 Pairwise Similarity & Molecular-Clock Summary (RQ3)",
        "",
        "## 1. Key Extremes & Identical Pairs",
        f"- **Most similar pair overall:** `{most_sim['sp1']}` vs `{most_sim['sp2']}` = **{most_sim['identity_pct']:.2f}%** identity ({int(most_sim['n_substitutions'])} substitutions).",
        f"- **Identical pairs (100.0% identity):** {len(ident_pairs)} pair(s) — "
        + (", ".join(f"`{r.sp1}` & `{r.sp2}`" for r in ident_pairs.itertuples()) if not ident_pairs.empty else "none") + ".",
        f"- **Least similar vertebrate pair:** `{least_vert['sp1']}` vs `{least_vert['sp2']}` = **{least_vert['identity_pct']:.2f}%** identity ({int(least_vert['n_substitutions'])} substitutions).",
        "",
        "## 2. Gradient Across Evolutionary Distance",
        f"- Mean identity **within mammals**: **{mam_mean:.2f}%**",
        f"- Mean identity **mammal vs. bird**: **{mb_mean:.2f}%**",
        f"- Mean identity **mammal vs. fish**: **{mf_mean:.2f}%**",
        f"- Mean identity **vertebrate vs. yeast**: **{vy_mean:.2f}%**",
        "",
        "## 3. Close vs. Distant Pair Highlights",
        df_to_markdown(pd.DataFrame(kp_rows)),
        "",
        "## 4. Identity to *Homo sapiens*",
        df_to_markdown(hum_df),
        "",
        "## 5. Identity vs. Divergence Time & Molecular-Clock Test",
        "Divergence times are approximate MRCA node ages (Mya) consistent with `config/species.csv` and TimeTree consensus estimates. "
        "**Caveat:** pairwise comparisons are **not statistically independent** due to shared phylogenetic ancestry, so p-values are descriptive only.",
        "",
        f"- **All 136 pairs — Identity vs. Time:** Spearman rho = **{clock_stats['all_id_rho']:.4f}** (p = {clock_stats['all_id_p']:.2e})",
        f"- **All 136 pairs — Poisson Distance vs. Time:** Spearman rho = **{clock_stats['all_poiss_rho']:.4f}** (p = {clock_stats['all_poiss_p']:.2e}); "
        f"Pearson R² = **{clock_stats['all_r2']:.4f}**, slope = **{clock_stats['all_slope']:.6f}** subs/site/Myr (**{clock_stats['all_slope']*1000:.3f}** subs/site/Gyr)",
        f"- **Human-vs-X (n = 16) — Identity vs. Time:** Spearman rho = **{clock_stats['hum_id_rho']:.4f}** (p = {clock_stats['hum_id_p']:.2e})",
        f"- **Human-vs-X (n = 16) — Poisson Distance vs. Time:** Spearman rho = **{clock_stats['hum_poiss_rho']:.4f}** (p = {clock_stats['hum_poiss_p']:.2e}); "
        f"Pearson R² = **{clock_stats['hum_r2']:.4f}**, slope = **{clock_stats['hum_slope']:.6f}** subs/site/Myr (**{clock_stats['hum_slope']*1000:.3f}** subs/site/Gyr)",
        "",
        "### Top 3 Molecular-Clock Outliers (Human-vs-X Residuals)",
        df_to_markdown(clock_stats["top_outliers"][["species", "clade", "mya_divergence", "poisson_distance", "predicted_distance", "residual"]]),
        "",
        "## 6. Group-Level Summary",
        df_to_markdown(group_sum),
        "",
    ]
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run Phase 5 pairwise similarity analysis, cross-validation, clock test, and figures."""
    parser = argparse.ArgumentParser(description="Phase 5: Pairwise similarity and molecular clock.")
    parser.add_argument("--config", type=str, default=None, help="Optional path to config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    logger = get_logger("50_similarity", cfg)
    logger.info(
        "Phase 5 goal: compute pairwise identity/similarity/distance matrices, cross-check against global pairwise "
        "alignment and Clustal Omega distmat, summarize by relationship, test identity-vs-time molecular clock, and plot figures. "
        "Inputs: results/alignments/primary_all.fasta, data/processed/cytc_curated_treeid.fasta, data/processed/id_map.csv."
    )

    aln_dir = resolve(cfg, "alignments")
    sim_dir = resolve(cfg, "similarity")
    tables_dir = resolve(cfg, "tables")
    fig_dir = resolve(cfg, "figures")
    proc_dir = resolve(cfg, "processed")

    aln_all = AlignIO.read(aln_dir / "primary_all.fasta", "fasta")
    id_map = pd.read_csv(proc_dir / "id_map.csv").sort_values("order_index")

    # Task 5.1: metric definitions
    write_metric_definitions(sim_dir / "metric_definitions.md", cfg)

    # Task 5.2 & 5.3: matrices and pairwise_long
    matrices, long_df = compute_msa_pairwise_metrics(aln_all, id_map, cfg)
    for name, df in matrices.items():
        df.to_csv(sim_dir / f"{name}.csv")
    long_df.to_csv(sim_dir / "pairwise_long.csv", index=False)

    # Cross-validation
    r_pw, max_diff_pw = run_cross_validation(
        proc_dir / "cytc_curated_treeid.fasta",
        long_df,
        matrices["identity_matrix"],
        aln_dir / "all_clustalo.distmat",
        cfg,
        sim_dir / "crosscheck.md",
        logger,
    )

    # Task 5.4: group summary and human-centric table
    group_sum = build_group_summary(long_df)
    group_sum.to_csv(tables_dir / "similarity_group_summary.csv", index=False)

    hum_df = build_identity_to_human(matrices, id_map)
    hum_df.to_csv(tables_dir / "identity_to_human.csv", index=False)

    # Task 5.5: molecular-clock correlation & regression
    all_id_rho, all_id_p = spearmanr(long_df["mya_divergence"], long_df["identity_pct"])
    all_poiss_rho, all_poiss_p = spearmanr(long_df["mya_divergence"], long_df["poisson_distance"])
    all_reg = linregress(long_df["mya_divergence"], long_df["poisson_distance"])

    hum_sub = hum_df[hum_df["species"] != "Homo_sapiens"].copy()
    hum_id_rho, hum_id_p = spearmanr(hum_sub["mya_divergence"], hum_sub["identity_pct"])
    hum_poiss_rho, hum_poiss_p = spearmanr(hum_sub["mya_divergence"], hum_sub["poisson_distance"])
    hum_reg = linregress(hum_sub["mya_divergence"], hum_sub["poisson_distance"])

    hum_sub["predicted_distance"] = (hum_reg.intercept + hum_reg.slope * hum_sub["mya_divergence"]).round(4)
    hum_sub["residual"] = (hum_sub["poisson_distance"] - hum_sub["predicted_distance"]).round(4)
    top_outliers = hum_sub.reindex(hum_sub["residual"].abs().sort_values(ascending=False).index).head(3)

    clock_stats = {
        "all_id_rho": float(all_id_rho),
        "all_id_p": float(all_id_p),
        "all_poiss_rho": float(all_poiss_rho),
        "all_poiss_p": float(all_poiss_p),
        "all_slope": float(all_reg.slope),
        "all_intercept": float(all_reg.intercept),
        "all_r2": float(all_reg.rvalue ** 2),
        "hum_id_rho": float(hum_id_rho),
        "hum_id_p": float(hum_id_p),
        "hum_poiss_rho": float(hum_poiss_rho),
        "hum_poiss_p": float(hum_poiss_p),
        "hum_slope": float(hum_reg.slope),
        "hum_intercept": float(hum_reg.intercept),
        "hum_r2": float(hum_reg.rvalue ** 2),
        "top_outliers": top_outliers,
    }
    logger.info(
        "Clock test: all pairs id rho=%.3f, poiss R^2=%.3f; human-vs-X id rho=%.3f, poiss R^2=%.3f, top outlier=%s (res=%+.3f)",
        all_id_rho,
        clock_stats["all_r2"],
        hum_id_rho,
        clock_stats["hum_r2"],
        top_outliers.iloc[0]["species"],
        top_outliers.iloc[0]["residual"],
    )

    write_summary_md(sim_dir / "summary.md", long_df, group_sum, hum_df, clock_stats)

    # Task 5.5b: Figures
    plot_annotated_heatmap(
        matrices["identity_matrix"],
        id_map,
        "Cytochrome c Pairwise Amino-Acid Identity (%) — Fixed Biological Order",
        "Pairwise Identity (%)",
        fig_dir / "fig_identity_heatmap",
        cfg,
    )
    plot_annotated_heatmap(
        matrices["similarity_matrix"],
        id_map,
        "Cytochrome c Pairwise BLOSUM62 Similarity (%) — Fixed Biological Order",
        "BLOSUM62 Similarity (%)",
        fig_dir / "fig_similarity_heatmap",
        cfg,
    )
    plot_clustermap(matrices["identity_matrix"], matrices["poisson_distance_matrix"], id_map, fig_dir, cfg)
    plot_identity_to_human(hum_df, fig_dir, cfg)
    plot_identity_vs_time(long_df, hum_df, clock_stats, fig_dir, cfg)
    plot_distribution_by_relationship(long_df, fig_dir, cfg)

    logger.info(
        "SUMMARY: computed 6 matrices & %d pairwise comparisons; crosscheck r=%.5f (max_diff=%.2f%%); "
        "human-vs-X clock R^2=%.3f; outputs=%s, %s",
        len(long_df),
        r_pw,
        max_diff_pw,
        clock_stats["hum_r2"],
        sim_dir,
        tables_dir / "identity_to_human.csv",
    )


if __name__ == "__main__":
    main()
