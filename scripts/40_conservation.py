"""Phase 4: Per-position conservation analysis, region detection, logos, and lineage substitutions."""
from __future__ import annotations

import argparse
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import logomaker
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from matplotlib.colors import ListedColormap

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cytc.utils import REPO_ROOT, get_logger, load_config, resolve  # noqa: E402

AA_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")

# BLOSUM62 background amino-acid frequencies (Capra & Singh 2007 / Henikoff & Henikoff 1992)
BLOSUM62_BG_RAW: dict[str, float] = {
    "A": 0.074, "R": 0.052, "N": 0.045, "D": 0.054, "C": 0.025,
    "Q": 0.034, "E": 0.054, "G": 0.074, "H": 0.026, "I": 0.068,
    "L": 0.099, "K": 0.058, "M": 0.025, "F": 0.047, "P": 0.039,
    "S": 0.057, "T": 0.051, "W": 0.013, "Y": 0.032, "V": 0.073,
}
_bg_total = sum(BLOSUM62_BG_RAW.values())
BLOSUM62_BG: dict[str, float] = {k: v / _bg_total for k, v in BLOSUM62_BG_RAW.items()}

# Physico-chemical classes (G, P, C treated as structurally distinct special singletons)
PHYSCHEM_MAP: dict[str, str] = {}
for _cls, _chars in (
    ("aliphatic", "AVLIM"),
    ("aromatic", "FWY"),
    ("positive", "KRH"),
    ("negative", "DE"),
    ("polar", "STNQ"),
    ("special_G", "G"),
    ("special_P", "P"),
    ("special_C", "C"),
):
    for _c in _chars:
        PHYSCHEM_MAP[_c] = _cls

CLASS_ORDER = ["invariant", "highly_conserved", "moderately_conserved", "variable", "gap_rich"]
CLASS_COLORS = {
    "invariant": "#0072B2",
    "highly_conserved": "#009E73",
    "moderately_conserved": "#E69F00",
    "variable": "#D55E00",
    "gap_rich": "#999999",
}


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


def build_position_map(aln_all: MultipleSeqAlignment, logger: Any) -> pd.DataFrame:
    """Map primary_all alignment columns to human UniProt, classical, and yeast positions."""
    hum_seq = str(next(r for r in aln_all if r.id == "Homo_sapiens").seq)
    yeast_seq = str(next(r for r in aln_all if r.id == "Saccharomyces_cerevisiae").seq)

    rows: list[dict[str, Any]] = []
    h_pos = 0
    y_pos = 0
    for col_idx in range(1, aln_all.get_alignment_length() + 1):
        h_ch = hum_seq[col_idx - 1]
        y_ch = yeast_seq[col_idx - 1]
        if h_ch != "-":
            h_pos += 1
            h_u: int | None = h_pos
            h_c: int | None = h_pos - 1 if h_pos > 1 else None
        else:
            h_u, h_c = None, None
        if y_ch != "-":
            y_pos += 1
            y_u: int | None = y_pos
        else:
            y_u = None

        rows.append(
            {
                "aln_col": col_idx,
                "human_uniprot_pos": h_u,
                "classical_pos": h_c,
                "human_residue": h_ch,
                "yeast_pos": y_u,
            }
        )

    pm = pd.DataFrame(rows)
    pm["human_uniprot_pos"] = pm["human_uniprot_pos"].astype("Int64")
    pm["classical_pos"] = pm["classical_pos"].astype("Int64")
    pm["yeast_pos"] = pm["yeast_pos"].astype("Int64")

    # Verify classical Cys14/Cys17/His18 <-> UniProt 15/18/19 and Met80 <-> UniProt 81
    checks = [(14, 15, "C"), (17, 18, "C"), (18, 19, "H"), (80, 81, "M")]
    for cl_pos, u_pos, expected_aa in checks:
        sub = pm[pm["classical_pos"] == cl_pos]
        if len(sub) != 1:
            raise ValueError(f"Coordinate mapping failed for classical pos {cl_pos}: {len(sub)} rows")
        r = sub.iloc[0]
        if int(r["human_uniprot_pos"]) != u_pos or str(r["human_residue"]) != expected_aa:
            raise ValueError(
                f"Anchor mismatch at classical {cl_pos}: got UniProt={r['human_uniprot_pos']}, "
                f"res={r['human_residue']} (expected UniProt={u_pos}, res={expected_aa})"
            )

    logger.info(
        "Position map verified (%d cols): Cys14/Cys17/His18 -> UniProt 15/18/19, Met80 -> UniProt 81.",
        len(pm),
    )
    return pm


def compute_henikoff_weights(aln: MultipleSeqAlignment) -> np.ndarray:
    """Compute Henikoff & Henikoff (1994) position-based sequence weights (sum to 1)."""
    n_seq = len(aln)
    aln_len = aln.get_alignment_length()
    seqs = [str(r.seq) for r in aln]
    weights = np.zeros(n_seq, dtype=float)

    for c in range(aln_len):
        col = [s[c] for s in seqs]
        non_gap = [ch for ch in col if ch != "-"]
        if not non_gap:
            continue
        counts = Counter(non_gap)
        r_types = len(counts)
        for i, ch in enumerate(col):
            if ch != "-":
                weights[i] += 1.0 / (r_types * counts[ch])

    return weights / weights.sum() if weights.sum() > 0 else np.full(n_seq, 1.0 / n_seq)


def jensen_shannon_divergence(p_dist: dict[str, float], q_dist: dict[str, float]) -> float:
    """Compute JSD(p || q) in bits (range [0, 1]) over the 20 standard amino acids."""
    jsd = 0.0
    for aa in AA_ALPHABET:
        p = p_dist.get(aa, 0.0)
        q = q_dist[aa]
        m = 0.5 * (p + q)
        if p > 0.0:
            jsd += 0.5 * p * math.log2(p / m)
        if q > 0.0:
            jsd += 0.5 * q * math.log2(q / m)
    return max(0.0, min(1.0, jsd))


def classify_column(
    identity_frac: float, n_gaps: int, gap_fraction: float, cfg: dict[str, Any]
) -> str:
    """Classify a single alignment column using unsmoothed config thresholds."""
    cc = cfg["conservation"]
    if gap_fraction > float(cc["max_gap_fraction_per_column"]):
        return "gap_rich"
    if identity_frac >= float(cc["invariant_threshold"]) and n_gaps == 0:
        return "invariant"
    if identity_frac >= float(cc["high_threshold"]):
        return "highly_conserved"
    if identity_frac >= float(cc["moderate_threshold"]):
        return "moderately_conserved"
    return "variable"


def compute_per_position(
    aln: MultipleSeqAlignment,
    aln_col_mapping: list[int],
    cfg: dict[str, Any],
    include_vert_col: bool = False,
    weights: np.ndarray | None = None,
) -> pd.DataFrame:
    """Compute unweighted (or weighted) per-column conservation metrics."""
    n_seq = len(aln)
    aln_len = aln.get_alignment_length()
    seqs = [str(r.seq) for r in aln]
    max_gap_frac = float(cfg["conservation"]["max_gap_fraction_per_column"])
    log2_20 = math.log2(20.0)

    if weights is None:
        w = np.full(n_seq, 1.0 / n_seq, dtype=float)
    else:
        w = weights

    rows: list[dict[str, Any]] = []
    for c in range(aln_len):
        col = [s[c] for s in seqs]
        n_gaps = sum(1 for ch in col if ch == "-")
        gap_frac = n_gaps / n_seq

        aa_weights: dict[str, float] = {}
        total_w = 0.0
        non_gap_chars: list[str] = []
        for i, ch in enumerate(col):
            if ch != "-":
                non_gap_chars.append(ch)
                aa_weights[ch] = aa_weights.get(ch, 0.0) + float(w[i])
                total_w += float(w[i])

        if total_w == 0.0:
            consensus = "-"
            id_frac = 0.0
            n_distinct = 0
            entropy = 0.0
            norm_ent = 0.0
            cons_score = 0.0
            jsd_score = 0.0
            physchem_cons = 0.0
        else:
            p_dist = {aa: val / total_w for aa, val in aa_weights.items()}
            # Deterministic tie-breaking: highest frequency, then alphabetical
            top_aa, top_freq = sorted(p_dist.items(), key=lambda kv: (-kv[1], kv[0]))[0]
            consensus = "-" if gap_frac > max_gap_frac else top_aa
            id_frac = top_freq
            n_distinct = len(p_dist)
            entropy = -sum(p * math.log2(p) for p in p_dist.values() if p > 0.0)
            entropy = max(0.0, entropy)
            norm_ent = entropy / log2_20
            cons_score = 1.0 - norm_ent
            jsd_score = (1.0 - gap_frac) * jensen_shannon_divergence(p_dist, BLOSUM62_BG)
            top_cls = PHYSCHEM_MAP[top_aa]
            physchem_cons = sum(p for aa, p in p_dist.items() if PHYSCHEM_MAP[aa] == top_cls)

        col_class = classify_column(id_frac, n_gaps, gap_frac, cfg)
        row: dict[str, Any] = {"aln_col": aln_col_mapping[c]}
        if include_vert_col:
            row["vert_aln_col"] = c + 1
        row.update(
            {
                "n_seqs": n_seq,
                "n_gaps": n_gaps,
                "gap_fraction": round(gap_frac, 4),
                "consensus": consensus,
                "identity_frac": round(id_frac, 6),
                "n_distinct_residues": n_distinct,
                "shannon_entropy_bits": round(entropy, 6),
                "norm_entropy": round(norm_ent, 6),
                "conservation_score": round(cons_score, 6),
                "jsd_score": round(jsd_score, 6),
                "physchem_conservation": round(physchem_cons, 6),
                "class": col_class,
            }
        )
        rows.append(row)

    df = pd.DataFrame(rows)
    win = int(cfg["conservation"]["smoothing_window"])
    df["smoothed_conservation"] = (
        df["conservation_score"].rolling(window=win, center=True, min_periods=1).mean().round(6)
    )
    return df


def detect_regions(
    per_pos_vert: pd.DataFrame, pos_map: pd.DataFrame, cfg: dict[str, Any]
) -> pd.DataFrame:
    """Detect strict contiguous conserved (>= min_region_length) and variable (>= 2) regions."""
    min_len = int(cfg["conservation"]["min_region_length"])
    merged = per_pos_vert.merge(pos_map, on="aln_col", how="left").sort_values("aln_col").reset_index(drop=True)

    regions: list[dict[str, Any]] = []
    for rtype, predicate, threshold in (
        ("conserved", lambda c: c in ("invariant", "highly_conserved"), min_len),
        ("variable", lambda c: c == "variable", 2),
    ):
        start_i: int | None = None
        for i, row in merged.iterrows():
            if predicate(row["class"]):
                if start_i is None:
                    start_i = i
            else:
                if start_i is not None:
                    if i - start_i >= threshold:
                        regions.append(_build_region_row(merged.iloc[start_i:i], rtype))
                    start_i = None
        if start_i is not None and len(merged) - start_i >= threshold:
            regions.append(_build_region_row(merged.iloc[start_i:], rtype))

    cols = [
        "region_type", "aln_start", "aln_end", "human_uniprot_start", "human_uniprot_end",
        "classical_start", "classical_end", "length", "mean_conservation", "sequence_in_human",
    ]
    return pd.DataFrame(regions, columns=cols)


def _build_region_row(sub: pd.DataFrame, region_type: str) -> dict[str, Any]:
    """Summarize a contiguous slice of merged per-position rows into a region record."""
    first, last = sub.iloc[0], sub.iloc[-1]
    return {
        "region_type": region_type,
        "aln_start": int(first["aln_col"]),
        "aln_end": int(last["aln_col"]),
        "human_uniprot_start": int(first["human_uniprot_pos"]) if pd.notna(first["human_uniprot_pos"]) else None,
        "human_uniprot_end": int(last["human_uniprot_pos"]) if pd.notna(last["human_uniprot_pos"]) else None,
        "classical_start": int(first["classical_pos"]) if pd.notna(first["classical_pos"]) else None,
        "classical_end": int(last["classical_pos"]) if pd.notna(last["classical_pos"]) else None,
        "length": int(len(sub)),
        "mean_conservation": round(float(sub["conservation_score"].mean()), 4),
        "sequence_in_human": "".join(str(x) for x in sub["human_residue"]),
    }


def find_invariant_stretches(df_merged: pd.DataFrame) -> list[dict[str, Any]]:
    """Return all maximal contiguous runs of invariant columns."""
    stretches: list[dict[str, Any]] = []
    start_i: int | None = None
    for i, row in df_merged.iterrows():
        if row["class"] == "invariant":
            if start_i is None:
                start_i = i
        else:
            if start_i is not None:
                stretches.append(_build_region_row(df_merged.iloc[start_i:i], "invariant"))
                start_i = None
    if start_i is not None:
        stretches.append(_build_region_row(df_merged.iloc[start_i:], "invariant"))
    return stretches


def find_lineage_specific_substitutions(
    aln_vert: MultipleSeqAlignment,
    vert_col_to_all_col: list[int],
    pos_map: pd.DataFrame,
    id_map: pd.DataFrame,
) -> pd.DataFrame:
    """Identify diagnostic substitutions shared by all members of a multi-species group and absent elsewhere."""
    sp_to_group = dict(zip(id_map["tree_id"], id_map["group"]))
    pm_idx = pos_map.set_index("aln_col")
    ids = [r.id for r in aln_vert]
    seqs = [str(r.seq) for r in aln_vert]

    group_members: dict[str, list[int]] = {}
    for i, sp in enumerate(ids):
        g = sp_to_group[sp]
        group_members.setdefault(g, []).append(i)
    multi_groups = {g: idxs for g, idxs in group_members.items() if len(idxs) >= 2}

    rows: list[dict[str, Any]] = []
    for c in range(aln_vert.get_alignment_length()):
        col = [s[c] for s in seqs]
        if len(set(col)) == 1:
            continue
        aln_col = vert_col_to_all_col[c]
        pm_row = pm_idx.loc[aln_col]
        for group, idxs in multi_groups.items():
            g_res = {col[i] for i in idxs}
            if len(g_res) != 1:
                continue
            diag = next(iter(g_res))
            if diag == "-":
                continue
            other_res = sorted({col[i] for i in range(len(ids)) if i not in idxs and col[i] != "-"})
            if diag not in other_res:
                rows.append(
                    {
                        "aln_col": aln_col,
                        "human_uniprot_pos": int(pm_row["human_uniprot_pos"]) if pd.notna(pm_row["human_uniprot_pos"]) else None,
                        "classical_pos": int(pm_row["classical_pos"]) if pd.notna(pm_row["classical_pos"]) else None,
                        "group": group,
                        "diagnostic_residue": diag,
                        "other_residues_seen": ",".join(other_res),
                    }
                )

    cols = ["aln_col", "human_uniprot_pos", "classical_pos", "group", "diagnostic_residue", "other_residues_seen"]
    return pd.DataFrame(rows, columns=cols)


def save_fig(fig: plt.Figure, stem: Path, cfg: dict[str, Any]) -> None:
    """Save figure in all configured formats and close it."""
    dpi = int(cfg["reporting"]["dpi"])
    for fmt in cfg["reporting"]["figure_formats"]:
        fig.savefig(stem.with_suffix(f".{fmt}"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_conservation_profile(
    vert_m: pd.DataFrame, all_m: pd.DataFrame, regions: pd.DataFrame, fig_dir: Path, cfg: dict[str, Any]
) -> None:
    """Plot 2-panel conservation profile (Vertebrates vs All-17) over human UniProt position."""
    fs = int(cfg["reporting"]["font_size"])
    plt.rcParams.update({"font.size": fs})
    fig, axes = plt.subplots(2, 1, figsize=(13, 7.6), sharex=True, sharey=True)

    cons_regions = regions[regions["region_type"] == "conserved"]
    all_hum = all_m[all_m["human_uniprot_pos"].notna()].copy()

    panels = [
        (axes[0], vert_m, "A. Vertebrates Only (n = 15, primary conservation basis)"),
        (axes[1], all_hum, "B. All Species Including Outgroups (n = 17, deep conservation)"),
    ]
    for ax, df, title in panels:
        x = df["human_uniprot_pos"].astype(int).to_numpy()
        y_raw = df["conservation_score"].to_numpy()
        y_smooth = df["smoothed_conservation"].to_numpy()

        for r in cons_regions.itertuples():
            ax.axvspan(r.human_uniprot_start - 0.5, r.human_uniprot_end + 0.5, color="#009E73", alpha=0.14, lw=0)

        ax.plot(x, y_raw, color="#56B4E9", linewidth=1.0, alpha=0.75, label="Raw conservation (1 - norm H)")
        ax.plot(x, y_smooth, color="#0072B2", linewidth=2.2, label=f"Smoothed (window={cfg['conservation']['smoothing_window']})")

        # Mark heme motif (UniProt 15, 18, 19) and Met81 (classical Met80)
        for pos, lbl, col in ((15, "C15", "#D55E00"), (18, "C18", "#D55E00"), (19, "H19", "#CC79A7"), (81, "M81 (Met80)", "#000000")):
            ax.axvline(pos, color=col, linestyle="--", linewidth=0.9, alpha=0.75)
        ax.text(17.0, 0.08, "Heme CxxCH\n(UniProt 15-19;\nclassical 14-18)", ha="center", fontsize=fs - 1.5,
                bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.85, "edgecolor": "#D55E00"})
        ax.text(81.0, 0.08, "Met81\n(classical Met80)", ha="center", fontsize=fs - 1.5,
                bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.85, "edgecolor": "black"})

        ax.set_ylabel("Conservation score")
        ax.set_ylim(-0.02, 1.08)
        ax.set_xlim(1, 105)
        ax.set_title(title, loc="left", fontweight="bold")
        ax.grid(True, linestyle=":", alpha=0.4)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Secondary top x-axis on Panel A for classical numbering
    secax = axes[0].secondary_xaxis("top", functions=(lambda u: u - 1, lambda c: c + 1))
    secax.set_xlabel("Classical (horse-heart) residue number (UniProt - 1)")
    secax.set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])

    axes[1].set_xlabel("Human UniProt residue position (Met1 = 1)")
    axes[1].set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 105])

    patch = mpatches.Patch(color="#009E73", alpha=0.22, label="Detected conserved region (vert)")
    handles, labels = axes[0].get_legend_handles_labels()
    handles.append(patch)
    axes[0].legend(handles=handles, loc="lower center", frameon=True, fontsize=fs - 1)

    fig.suptitle("Cytochrome c Per-Position Conservation Profile", y=1.02, fontsize=fs + 3, fontweight="bold")
    save_fig(fig, fig_dir / "fig_conservation_profile", cfg)


def plot_entropy_bar(vert_m: pd.DataFrame, fig_dir: Path, cfg: dict[str, Any]) -> None:
    """Plot per-column Shannon entropy (bits) coloured by conservation class (vertebrates)."""
    fs = int(cfg["reporting"]["font_size"])
    fig, ax = plt.subplots(figsize=(13.5, 4.5))
    x = vert_m["human_uniprot_pos"].astype(int).to_numpy()
    y = vert_m["shannon_entropy_bits"].to_numpy()
    colors = [CLASS_COLORS[c] for c in vert_m["class"]]

    ax.bar(x, y, color=colors, width=0.85, edgecolor="none")
    ax.set_xlim(0.5, 105.5)
    ax.set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 105])
    ax.set_xlabel("Human UniProt residue position (Met1 = 1)")
    ax.set_ylabel("Shannon entropy (bits)")
    ax.set_title("Per-Position Shannon Entropy in Vertebrate Cytochrome c (n = 15)")
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    present_classes = [c for c in CLASS_ORDER if c in set(vert_m["class"])]
    handles = [mpatches.Patch(color=CLASS_COLORS[c], label=c.replace("_", " ").title()) for c in present_classes]
    ax.legend(handles=handles, title="Conservation Class", loc="upper left", frameon=True)
    save_fig(fig, fig_dir / "fig_entropy_bar", cfg)


def plot_sequence_logo(
    aln: MultipleSeqAlignment,
    pos_labels: list[str],
    title: str,
    stem: Path,
    cfg: dict[str, Any],
) -> None:
    """Generate a 3-row split logomaker information-content sequence logo (pseudocount=0)."""
    fs = int(cfg["reporting"]["font_size"])
    seqs = [str(r.seq) for r in aln]
    info_df = logomaker.alignment_to_matrix(
        seqs, to_type="information", characters_to_ignore="-", pseudocount=0
    )
    n_cols = len(info_df)
    n_rows = 3
    chunk = math.ceil(n_cols / n_rows)

    fig, axes = plt.subplots(n_rows, 1, figsize=(13, 6.8), sharey=True)
    for row_idx in range(n_rows):
        start = row_idx * chunk
        end = min((row_idx + 1) * chunk, n_cols)
        sub_df = info_df.iloc[start:end].reset_index(drop=True)
        ax = axes[row_idx]
        logomaker.Logo(sub_df, color_scheme="chemistry", ax=ax)
        ax.set_ylim(0, math.log2(20.0) + 0.1)
        ax.set_ylabel("Bits")
        ticks = list(range(len(sub_df)))
        tick_lbls = pos_labels[start:end]
        # Show every position or every 2nd if dense
        ax.set_xticks(ticks)
        ax.set_xticklabels(tick_lbls, fontsize=max(7, fs - 2.5), rotation=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[-1].set_xlabel("Human UniProt position (outgroup-only N-terminal columns prefixed with 'N')")
    fig.suptitle(title, y=0.99, fontsize=fs + 2.5, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_fig(fig, stem, cfg)


def plot_class_bar(vert_df: pd.DataFrame, all_df: pd.DataFrame, fig_dir: Path, cfg: dict[str, Any]) -> None:
    """Grouped bar chart of conservation class counts (Vertebrates vs All-17)."""
    fs = int(cfg["reporting"]["font_size"])
    v_counts = [int((vert_df["class"] == c).sum()) for c in CLASS_ORDER]
    a_counts = [int((all_df["class"] == c).sum()) for c in CLASS_ORDER]

    x = np.arange(len(CLASS_ORDER))
    w = 0.36
    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    b1 = ax.bar(x - w / 2, v_counts, width=w, color="#0072B2", edgecolor="black", linewidth=0.6, label=f"Vertebrates (n=15, {len(vert_df)} cols)")
    b2 = ax.bar(x + w / 2, a_counts, width=w, color="#E69F00", edgecolor="black", linewidth=0.6, label=f"All Species (n=17, {len(all_df)} cols)")

    for bars, total in ((b1, len(vert_df)), (b2, len(all_df))):
        for bar in bars:
            h = bar.get_height()
            pct = 100.0 * h / total
            ax.text(bar.get_x() + bar.get_width() / 2, h + 1.0, f"{h}\n({pct:.0f}%)", ha="center", va="bottom", fontsize=fs - 1.5)

    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("_", "\n").title() for c in CLASS_ORDER])
    ax.set_ylabel("Number of alignment columns")
    ax.set_ylim(0, max(v_counts + a_counts) + 14)
    ax.set_title("Distribution of Alignment Columns by Conservation Class")
    ax.legend(frameon=True)
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save_fig(fig, fig_dir / "fig_conservation_class_bar", cfg)


def plot_species_consensus_heatmap(
    aln_all: MultipleSeqAlignment,
    vert_m: pd.DataFrame,
    fig_dir: Path,
    cfg: dict[str, Any],
) -> None:
    """Binary heatmap of species x non-invariant vertebrate columns (match vs mismatch to vertebrate consensus)."""
    fs = int(cfg["reporting"]["font_size"])
    var_cols = vert_m[vert_m["class"] != "invariant"].copy()
    col_indices_0 = [int(c) - 1 for c in var_cols["aln_col"]]
    cons_residues = list(var_cols["consensus"])
    col_labels = [f"{int(u)} ({c})" for u, c in zip(var_cols["human_uniprot_pos"], cons_residues)]

    n_seq = len(aln_all)
    mat = np.zeros((n_seq, len(col_indices_0)), dtype=int)
    annots: list[list[str]] = []
    for i, rec in enumerate(aln_all):
        s = str(rec.seq)
        row_ann: list[str] = []
        for j, (c0, cons) in enumerate(zip(col_indices_0, cons_residues)):
            ch = s[c0]
            if ch == "-":
                mat[i, j] = 2
                row_ann.append("-")
            elif ch == cons:
                mat[i, j] = 1
                row_ann.append("·")
            else:
                mat[i, j] = 0
                row_ann.append(ch)
        annots.append(row_ann)

    cmap = ListedColormap(["#D55E00", "#E8F4FA", "#CCCCCC"])
    fig, ax = plt.subplots(figsize=(max(12, 0.36 * len(col_indices_0)), 6.8))
    ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=2, interpolation="nearest")

    for i in range(n_seq):
        for j in range(len(col_indices_0)):
            txt = annots[i][j]
            color = "white" if mat[i, j] == 0 else "#333333"
            weight = "bold" if mat[i, j] == 0 else "normal"
            ax.text(j, i, txt, ha="center", va="center", color=color, fontsize=max(7, fs - 2.5), fontweight=weight)

    ax.set_yticks(np.arange(n_seq))
    ax.set_yticklabels([r.id.replace("_", " ") for r in aln_all])
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=90, fontsize=max(7.5, fs - 2))
    ax.set_xlabel("Human UniProt position (vertebrate consensus residue)")
    ax.set_ylabel("Species (Biological Order)")
    ax.set_title("Substitutions Relative to Vertebrate Consensus Across Non-Invariant Vertebrate Positions")

    handles = [
        mpatches.Patch(facecolor="#E8F4FA", edgecolor="black", linewidth=0.5, label="Matches vertebrate consensus (·)"),
        mpatches.Patch(facecolor="#D55E00", edgecolor="black", linewidth=0.5, label="Substitution (residue shown)"),
        mpatches.Patch(facecolor="#CCCCCC", edgecolor="black", linewidth=0.5, label="Gap (-)"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, frameon=True)
    save_fig(fig, fig_dir / "fig_conservation_heatmap_per_species", cfg)


def summarize_dataset(
    name: str, df_m: pd.DataFrame
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Compute headline summary stats and invariant stretches for a dataset."""
    total = len(df_m)
    gapfree = int((df_m["n_gaps"] == 0).sum())
    counts = {c: int((df_m["class"] == c).sum()) for c in CLASS_ORDER}
    ge90 = counts["invariant"] + counts["highly_conserved"]
    stretches = find_invariant_stretches(df_m)
    longest = max(stretches, key=lambda s: s["length"]) if stretches else None

    summary = {
        "dataset": name,
        "total_columns": total,
        "gapfree_columns": gapfree,
        "n_invariant": counts["invariant"],
        "pct_invariant": round(100.0 * counts["invariant"] / total, 2),
        "n_highly_conserved": counts["highly_conserved"],
        "n_ge_90_conserved": ge90,
        "pct_ge_90_conserved": round(100.0 * ge90 / total, 2),
        "n_moderately_conserved": counts["moderately_conserved"],
        "pct_moderately_conserved": round(100.0 * counts["moderately_conserved"] / total, 2),
        "n_variable": counts["variable"],
        "pct_variable": round(100.0 * counts["variable"] / total, 2),
        "n_gap_rich": counts["gap_rich"],
        "mean_identity_frac": round(float(df_m["identity_frac"].mean()), 4),
        "median_identity_frac": round(float(df_m["identity_frac"].median()), 4),
        "mean_entropy_bits": round(float(df_m["shannon_entropy_bits"].mean()), 4),
        "mean_conservation_score": round(float(df_m["conservation_score"].mean()), 4),
        "mean_jsd_score": round(float(df_m["jsd_score"].mean()), 4),
        "n_invariant_stretches_ge2": sum(1 for s in stretches if s["length"] >= 2),
        "longest_invariant_stretch_len": longest["length"] if longest else 0,
        "longest_invariant_stretch_uniprot": (
            f"{longest['human_uniprot_start']}-{longest['human_uniprot_end']}" if longest else "none"
        ),
    }
    return summary, stretches


def residue_distribution_str(aln: MultipleSeqAlignment, col_0based: int) -> str:
    """Return formatted residue distribution string like 'K:9, R:4, Q:2'."""
    cnt = Counter(str(r.seq)[col_0based] for r in aln if str(r.seq)[col_0based] != "-")
    return ", ".join(f"{aa}:{n}" for aa, n in sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0])))


def write_summary_md(
    out_path: Path,
    summary_df: pd.DataFrame,
    vert_m: pd.DataFrame,
    all_m: pd.DataFrame,
    aln_vert: MultipleSeqAlignment,
    regions: pd.DataFrame,
    lineage_df: pd.DataFrame,
    henikoff_diffs: int,
    diffs_vs_human: list[tuple[str, int]],
) -> None:
    """Write human-readable Phase 4 summary (RQ1 & RQ2 answers) to results/conservation/summary.md."""
    # Top 10 most variable vertebrate positions
    top_var = vert_m.sort_values(["shannon_entropy_bits", "identity_frac"], ascending=[False, True]).head(10)
    var_rows = []
    for r in top_var.itertuples():
        dist = residue_distribution_str(aln_vert, int(r.vert_aln_col) - 1)
        var_rows.append(
            {
                "aln_col": int(r.aln_col),
                "human_uniprot_pos": int(r.human_uniprot_pos),
                "classical_pos": int(r.classical_pos) if pd.notna(r.classical_pos) else "-",
                "entropy_bits": round(float(r.shannon_entropy_bits), 3),
                "identity_frac": round(float(r.identity_frac), 3),
                "residues": dist,
            }
        )
    top_var_df = pd.DataFrame(var_rows)

    # Deep conservation: vertebrate invariant sites that stay invariant across all 17 species
    vert_inv_cols = set(vert_m.loc[vert_m["class"] == "invariant", "aln_col"])
    all_inv_df = all_m[all_m["class"] == "invariant"].copy()
    all_inv_cols = set(all_inv_df["aln_col"])
    deep_overlap = vert_inv_cols & all_inv_cols
    deep_pct = 100.0 * len(deep_overlap) / len(vert_inv_cols) if vert_inv_cols else 0.0

    deep_list = ", ".join(
        f"{r.human_residue}{int(r.human_uniprot_pos)} (cl.{int(r.classical_pos) if pd.notna(r.classical_pos) else 'M1'})"
        for r in all_inv_df.itertuples()
    )

    diff_df = pd.DataFrame(diffs_vs_human, columns=["species", "aa_differences_from_human"])

    lines = [
        "# Phase 4 Conservation Analysis Summary (RQ1 & RQ2)",
        "",
        "## 1. Headline Conservation Statistics (RQ1)",
        "- Sequence weights: **unweighted** (n = 15 vertebrates, n = 17 all; no pseudo-count smoothing applied).",
        f"- Henikoff sequence-weighting sensitivity check (vertebrates): **{henikoff_diffs}** / {len(vert_m)} columns changed conservation class under Henikoff weights.",
        "- JSD background: published BLOSUM62 amino-acid frequencies (Capra & Singh 2007, window = 0, multiplied by `1 - gap_fraction`).",
        "",
        df_to_markdown(summary_df),
        "",
        "## 2. Deep Conservation (Vertebrates -> All 17 Species)",
        f"- **{len(deep_overlap)} / {len(vert_inv_cols)} ({deep_pct:.1f}%)** of vertebrate-invariant positions remain strictly invariant when *D. melanogaster* and *S. cerevisiae* are included.",
        f"- Total strictly invariant positions across all 17 species: **{len(all_inv_df)} / {len(all_m)} ({100.0 * len(all_inv_df) / len(all_m):.1f}%)**.",
        f"- All-17 invariant residues (UniProt-human / classical): {deep_list}.",
        "",
        "## 3. Detected Conserved and Variable Regions (Vertebrates, RQ2)",
        "Conserved regions = maximal runs of >= 3 consecutive `invariant` or `highly_conserved` columns (strict, unmerged). "
        "Variable regions = maximal runs of >= 2 consecutive `variable` columns.",
        "",
        df_to_markdown(regions) if not regions.empty else "_No regions detected._",
        "",
        "## 4. Top 10 Most Variable Vertebrate Positions",
        df_to_markdown(top_var_df),
        "",
        "## 5. Diagnostic Lineage-Specific Substitutions (Multi-Species Groups)",
        df_to_markdown(lineage_df) if not lineage_df.empty else "_No strictly diagnostic multi-species substitutions found._",
        "",
        "## 6. Differences vs. *Homo sapiens* (Vertebrates Quick Check)",
        df_to_markdown(diff_df),
        "",
    ]
    out_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run Phase 4 conservation analysis, region detection, logos, and lineage substitution analysis."""
    parser = argparse.ArgumentParser(description="Phase 4: Conservation analysis and sequence logos.")
    parser.add_argument("--config", type=str, default=None, help="Optional path to config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    logger = get_logger("40_conservation", cfg)
    logger.info(
        "Phase 4 goal: map alignment coordinates, compute per-position conservation metrics (unweighted + Henikoff check), "
        "detect conserved/variable regions, identify lineage-specific substitutions, and render conservation figures/logos. "
        "Inputs: results/alignments/primary_vert.fasta, results/alignments/primary_all.fasta, data/processed/id_map.csv."
    )

    aln_dir = resolve(cfg, "alignments")
    cons_dir = resolve(cfg, "conservation")
    tables_dir = resolve(cfg, "tables")
    fig_dir = resolve(cfg, "figures")
    proc_dir = resolve(cfg, "processed")

    aln_all = AlignIO.read(aln_dir / "primary_all.fasta", "fasta")
    aln_vert = AlignIO.read(aln_dir / "primary_vert.fasta", "fasta")
    id_map = pd.read_csv(proc_dir / "id_map.csv").sort_values("order_index")

    # Task 4.1: Coordinate map
    pos_map = build_position_map(aln_all, logger)
    pos_map_path = tables_dir / "position_map.csv"
    pos_map.to_csv(pos_map_path, index=False)

    # Map each primary_vert column (1..105) to its canonical primary_all aln_col via Homo_sapiens UniProt position
    hum_vert = str(next(r for r in aln_vert if r.id == "Homo_sapiens").seq)
    u_to_all_col = {
        int(r.human_uniprot_pos): int(r.aln_col)
        for r in pos_map.itertuples()
        if pd.notna(r.human_uniprot_pos)
    }
    vert_col_to_all_col: list[int] = []
    u_counter = 0
    for ch in hum_vert:
        if ch != "-":
            u_counter += 1
            vert_col_to_all_col.append(u_to_all_col[u_counter])
        else:
            raise ValueError("Unexpected gap in Homo_sapiens inside primary_vert.fasta")

    all_col_list = list(range(1, aln_all.get_alignment_length() + 1))

    # Task 4.2: Per-position metrics
    per_vert = compute_per_position(aln_vert, vert_col_to_all_col, cfg, include_vert_col=True)
    per_all = compute_per_position(aln_all, all_col_list, cfg, include_vert_col=False)

    per_vert_path = cons_dir / "per_position_vert.csv"
    per_all_path = cons_dir / "per_position_all.csv"
    per_vert.to_csv(per_vert_path, index=False)
    per_all.to_csv(per_all_path, index=False)

    # Sensitivity check: Henikoff weights on vertebrates
    h_weights = compute_henikoff_weights(aln_vert)
    per_vert_weighted = compute_per_position(aln_vert, vert_col_to_all_col, cfg, include_vert_col=True, weights=h_weights)
    henikoff_diffs = int((per_vert["class"] != per_vert_weighted["class"]).sum())
    logger.info("Henikoff sequence-weighting sensitivity check: %d/%d vertebrate columns changed class.", henikoff_diffs, len(per_vert))

    vert_m = per_vert.merge(pos_map, on="aln_col", how="left")
    all_m = per_all.merge(pos_map, on="aln_col", how="left")

    # Quick check: differences vs human
    hum_s = str(next(r for r in aln_vert if r.id == "Homo_sapiens").seq)
    diffs_vs_human: list[tuple[str, int]] = []
    for r in aln_vert:
        s = str(r.seq)
        n_diff = sum(1 for a, b in zip(hum_s, s) if a != b)
        diffs_vs_human.append((r.id, n_diff))
    logger.info("Vertebrate differences vs human: %s", ", ".join(f"{sp}:{d}" for sp, d in diffs_vs_human))

    # Task 4.3: Regions
    regions_df = detect_regions(per_vert, pos_map, cfg)
    regions_path = cons_dir / "regions.csv"
    regions_df.to_csv(regions_path, index=False)

    # Task 4.4: Summary stats
    vert_sum, _ = summarize_dataset("vertebrates_only", vert_m)
    all_sum, _ = summarize_dataset("all_17_species", all_m)
    summary_df = pd.DataFrame([vert_sum, all_sum])
    summary_csv_path = tables_dir / "conservation_summary.csv"
    summary_df.to_csv(summary_csv_path, index=False)

    # Task 4.6: Lineage-specific substitutions
    lineage_df = find_lineage_specific_substitutions(aln_vert, vert_col_to_all_col, pos_map, id_map)
    lineage_path = tables_dir / "lineage_specific_substitutions.csv"
    lineage_df.to_csv(lineage_path, index=False)

    # Write summary.md
    summary_md_path = cons_dir / "summary.md"
    write_summary_md(
        summary_md_path, summary_df, vert_m, all_m, aln_vert, regions_df, lineage_df, henikoff_diffs, diffs_vs_human
    )

    # Task 4.5: Figures
    plot_conservation_profile(vert_m, all_m, regions_df, fig_dir, cfg)
    plot_entropy_bar(vert_m, fig_dir, cfg)

    vert_labels = [str(int(u)) for u in vert_m["human_uniprot_pos"]]
    plot_sequence_logo(
        aln_vert,
        vert_labels,
        "Cytochrome c Sequence Logo — Vertebrates Only (n = 15, Human UniProt Numbering)",
        fig_dir / "fig_sequence_logo_vert",
        cfg,
    )

    all_labels = [
        str(int(r.human_uniprot_pos)) if pd.notna(r.human_uniprot_pos) else f"N{int(r.aln_col)}"
        for r in all_m.itertuples()
    ]
    plot_sequence_logo(
        aln_all,
        all_labels,
        "Cytochrome c Sequence Logo — All 17 Species Including Outgroups (Human UniProt Numbering)",
        fig_dir / "fig_sequence_logo_all",
        cfg,
    )

    plot_class_bar(per_vert, per_all, fig_dir, cfg)
    plot_species_consensus_heatmap(aln_all, vert_m, fig_dir, cfg)

    logger.info(
        "SUMMARY: vert invariant=%d/%d (%.1f%%), >=90%%=%d (%.1f%%); all invariant=%d/%d (%.1f%%); "
        "regions=%d, diagnostic_subs=%d; outputs=%s, %s, %s",
        vert_sum["n_invariant"],
        vert_sum["total_columns"],
        vert_sum["pct_invariant"],
        vert_sum["n_ge_90_conserved"],
        vert_sum["pct_ge_90_conserved"],
        all_sum["n_invariant"],
        all_sum["total_columns"],
        all_sum["pct_invariant"],
        len(regions_df),
        len(lineage_df),
        per_vert_path,
        regions_path,
        summary_csv_path,
    )


if __name__ == "__main__":
    main()

