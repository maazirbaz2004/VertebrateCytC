"""Phase 2: QC, standardize, and curate cytochrome c sequences."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import Bio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cytc.utils import REPO_ROOT, get_logger, load_config, resolve  # noqa: E402

GROUP_ORDER = [
    "Primate",
    "Rodent",
    "Carnivore",
    "Artiodactyl",
    "Perissodactyl",
    "Marsupial",
    "Bird",
    "Reptile",
    "Amphibian",
    "Fish",
    "Insect",
    "Fungus",
]
TREE_ID_RE = re.compile(r"^[A-Za-z0-9_]+$")
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")


def clean_sequence(raw_seq: str) -> str:
    """Upper-case sequence and strip whitespace, stop symbols (*), and gaps."""
    return re.sub(r"[\s*\-.]", "", raw_seq.upper())


def sha256_file(path: Path) -> str:
    """Return hex SHA-256 digest of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def write_wrapped_fasta(records: list[SeqRecord], out_path: Path, width: int = 60) -> None:
    """Write FASTA records wrapped to a fixed column width (default 60)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
        for rec in records:
            fh.write(f">{rec.id}\n")
            seq_str = str(rec.seq)
            for i in range(0, len(seq_str), width):
                fh.write(seq_str[i : i + width] + "\n")


def load_and_validate_records(
    raw_fasta: Path,
    species_df: pd.DataFrame,
    prov_df: pd.DataFrame,
    cfg: dict[str, Any],
    logger: Any,
) -> list[dict[str, Any]]:
    """Parse raw FASTA, validate sequences, check duplicates, and sort biologically."""
    if not raw_fasta.exists():
        raise FileNotFoundError(f"Missing raw FASTA: {raw_fasta}")

    allowed = set(cfg["sequence_prep"]["allowed_alphabet"])
    max_ambig = float(cfg["sequence_prep"]["max_ambiguous_fraction"])
    motif_re = re.compile(cfg["data_collection"]["heme_motif_regex"])
    dc = cfg["data_collection"]

    species_meta = species_df.set_index("species").to_dict(orient="index")
    species_csv_rank = {sp: idx for idx, sp in enumerate(species_df["species"])}
    prov_acc = dict(zip(prov_df["species"], prov_df["accession"]))

    raw_records = list(SeqIO.parse(raw_fasta, "fasta"))
    if not raw_records:
        raise ValueError(f"No records found in {raw_fasta}")

    seen_species: set[str] = set()
    parsed: list[dict[str, Any]] = []

    for rec in raw_records:
        parts = rec.id.split("|")
        if len(parts) < 2:
            raise ValueError(f"Unexpected raw FASTA header format: {rec.id}")
        species, accession = parts[0].strip(), parts[1].strip()

        if species in seen_species:
            raise ValueError(f"Duplicate record within species: {species}")
        seen_species.add(species)

        if species not in species_meta:
            raise ValueError(f"Species in raw FASTA not found in species.csv: {species}")
        if not TREE_ID_RE.match(species):
            raise ValueError(f"Invalid tree_id characters in species name: {species}")
        if species in prov_acc and prov_acc[species] != accession:
            raise ValueError(
                f"Accession mismatch for {species}: FASTA={accession} vs provenance={prov_acc[species]}"
            )

        seq = clean_sequence(str(rec.seq))
        bad_chars = [ch for ch in seq if ch not in allowed]
        ambig_frac = len(bad_chars) / len(seq) if seq else 1.0
        if ambig_frac > max_ambig or bad_chars:
            raise ValueError(
                f"{species} ({accession}) has invalid residues {sorted(set(bad_chars))} "
                f"(fraction={ambig_frac:.4f} > {max_ambig})"
            )

        match = motif_re.search(seq)
        if not match:
            raise ValueError(f"{species} ({accession}) missing heme motif {motif_re.pattern}")
        motif_str = match.group(0)
        motif_start = match.start() + 1  # 1-based UniProt numbering

        meta = species_meta[species]
        role = str(meta["role"]).strip().lower()
        lo = int(dc["outgroup_len_min"] if role == "outgroup" else dc["vertebrate_len_min"])
        hi = int(dc["outgroup_len_max"] if role == "outgroup" else dc["vertebrate_len_max"])
        if not (lo <= len(seq) <= hi):
            raise ValueError(
                f"{species} ({accession}) length {len(seq)} outside [{lo}, {hi}] for {role}"
            )

        parsed.append(
            {
                "species": species,
                "accession": accession,
                "group": str(meta["group"]).strip(),
                "role": role,
                "taxon_id": int(meta["taxon_id"]),
                "common_name": str(meta["common_name"]).strip(),
                "sequence": seq,
                "length": len(seq),
                "heme_motif": motif_str,
                "heme_motif_start": motif_start,
            }
        )

    missing_species = set(species_df["species"]) - seen_species
    if missing_species:
        raise ValueError(f"Missing species from raw FASTA: {sorted(missing_species)}")

    group_rank = {g: i for i, g in enumerate(GROUP_ORDER)}
    parsed.sort(
        key=lambda r: (
            0 if r["role"] == "ingroup" else 1,
            group_rank.get(r["group"], len(GROUP_ORDER)),
            species_csv_rank[r["species"]],
        )
    )
    for idx, row in enumerate(parsed, start=1):
        row["order_index"] = idx
        classical = row["heme_motif_start"] - 1 if row["role"] == "ingroup" else None
        classical_note = f", classical={classical}" if classical is not None else ""
        logger.info(
            "Validated [%02d] %-24s %-10s len=%3d motif=%s (UniProt=%d%s)",
            idx,
            row["species"],
            row["accession"],
            row["length"],
            row["heme_motif"],
            row["heme_motif_start"],
            classical_note,
        )

    return parsed


def compute_qc_tables(
    records: list[dict[str, Any]], logger: Any
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute seq_qc, aa_composition (wide & long), length_summary, and identical_sequence_groups."""
    qc_rows: list[dict[str, Any]] = []
    aa_wide_rows: list[dict[str, Any]] = []
    aa_long_rows: list[dict[str, Any]] = []
    seq_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for rec in records:
        seq = rec["sequence"]
        seq_groups[seq].append(rec)
        pa = ProteinAnalysis(seq)
        counts = Counter(seq)
        length = rec["length"]
        n_lys = counts["K"]
        n_arg = counts["R"]
        n_his = counts["H"]
        n_cys = counts["C"]
        n_met = counts["M"]
        # Charged residues at physiological pH: K, R, D, E
        fraction_charged = (counts["K"] + counts["R"] + counts["D"] + counts["E"]) / length
        mw_kda = pa.molecular_weight() / 1000.0
        pi = pa.isoelectric_point()
        gravy = pa.gravy()

        qc_rows.append(
            {
                "species": rec["species"],
                "accession": rec["accession"],
                "group": rec["group"],
                "role": rec["role"],
                "length": length,
                "n_term_residue": seq[0],
                "c_term_residue": seq[-1],
                "heme_motif": rec["heme_motif"],
                "heme_motif_start": rec["heme_motif_start"],
                "mw_kda": round(mw_kda, 3),
                "isoelectric_point": round(pi, 3),
                "gravy": round(gravy, 4),
                "n_lys": n_lys,
                "n_arg": n_arg,
                "n_his": n_his,
                "n_cys": n_cys,
                "n_met": n_met,
                "fraction_charged": round(fraction_charged, 4),
            }
        )

        wide_row: dict[str, Any] = {
            "species": rec["species"],
            "group": rec["group"],
            "role": rec["role"],
        }
        for aa in AA_ORDER:
            pct = 100.0 * counts[aa] / length
            wide_row[aa] = round(pct, 3)
            aa_long_rows.append(
                {
                    "species": rec["species"],
                    "group": rec["group"],
                    "role": rec["role"],
                    "amino_acid": aa,
                    "count": counts[aa],
                    "percent": round(pct, 3),
                }
            )
        aa_wide_rows.append(wide_row)

    qc_df = pd.DataFrame(qc_rows)
    aa_wide_df = pd.DataFrame(aa_wide_rows)
    aa_long_df = pd.DataFrame(aa_long_rows)

    # Identical sequences across species
    ident_rows: list[dict[str, Any]] = []
    gid = 0
    for seq, members in seq_groups.items():
        if len(members) > 1:
            gid += 1
            ident_rows.append(
                {
                    "group_id": f"IDENT_{gid:02d}",
                    "n_members": len(members),
                    "length": len(seq),
                    "members": ";".join(m["species"] for m in members),
                    "accessions": ";".join(m["accession"] for m in members),
                }
            )
    ident_df = pd.DataFrame(
        ident_rows, columns=["group_id", "n_members", "length", "members", "accessions"]
    )
    if ident_rows:
        for r in ident_rows:
            logger.info("Identical sequence group %s (%d species): %s", r["group_id"], r["n_members"], r["members"])
    else:
        logger.info("No identical sequences across species.")

    # Length summary
    length_summary_df = build_length_summary(qc_df)

    # Biological expectation checks (warn, not fail)
    check_biological_expectations(qc_df, logger)

    return qc_df, aa_wide_df, aa_long_df, length_summary_df, ident_df


def build_length_summary(qc_df: pd.DataFrame) -> pd.DataFrame:
    """Build min/max/mean/median/mode length summary for all, vertebrates, and outgroups."""
    vert_df = qc_df[qc_df["role"] == "ingroup"]
    vert_mode = int(vert_df["length"].mode().iloc[0])

    rows: list[dict[str, Any]] = []
    subsets = [
        ("all", qc_df),
        ("ingroup_vertebrates", vert_df),
        ("outgroups", qc_df[qc_df["role"] == "outgroup"]),
    ]
    for label, sub in subsets:
        lens = sub["length"]
        mode_val = int(lens.mode().iloc[0])
        n_at_mode = int((lens == mode_val).sum())
        if label == "ingroup_vertebrates":
            non_modal = sub[sub["length"] != mode_val]
            diff_str = (
                ", ".join(f"{r.species} ({r.length} aa)" for r in non_modal.itertuples())
                if not non_modal.empty
                else "none"
            )
            note = f"{n_at_mode}/{len(sub)} vertebrates have modal length {mode_val} aa; non-modal: {diff_str}"
        elif label == "outgroups":
            diffs = [
                f"{r.species}: {r.length} aa ({r.length - vert_mode:+d} aa vs vertebrate mode {vert_mode})"
                for r in sub.itertuples()
            ]
            note = "; ".join(diffs)
        else:
            note = f"{len(sub)} total sequences ({len(vert_df)} vertebrates + {len(sub) - len(vert_df)} outgroups)"

        rows.append(
            {
                "subset": label,
                "n_species": int(len(sub)),
                "min_length": int(lens.min()),
                "max_length": int(lens.max()),
                "mean_length": round(float(lens.mean()), 2),
                "median_length": round(float(lens.median()), 2),
                "mode_length": mode_val,
                "n_at_mode": n_at_mode,
                "notes": note,
            }
        )
    return pd.DataFrame(rows)


def check_biological_expectations(qc_df: pd.DataFrame, logger: Any) -> None:
    """Log biological QC expectations for vertebrate cytochrome c and outgroups."""
    verts = qc_df[qc_df["role"] == "ingroup"]
    vert_mode = int(verts["length"].mode().iloc[0])

    low_pi = verts[verts["isoelectric_point"] < 8.5]
    if low_pi.empty:
        logger.info(
            "QC CHECK [pI]: All %d vertebrates are basic (pI range %.2f-%.2f).",
            len(verts),
            verts["isoelectric_point"].min(),
            verts["isoelectric_point"].max(),
        )
    else:
        for r in low_pi.itertuples():
            logger.warning("QC WARN [pI < 8.5]: %s has pI=%.2f", r.species, r.isoelectric_point)

    mw_outliers = verts[(verts["mw_kda"] < 11.5) | (verts["mw_kda"] > 12.0)]
    logger.info(
        "QC CHECK [MW]: Vertebrate MW range %.3f-%.3f kDa (expected ~11.5-12.0 kDa).",
        verts["mw_kda"].min(),
        verts["mw_kda"].max(),
    )
    for r in mw_outliers.itertuples():
        logger.warning("QC WARN [MW outside 11.5-12.0 kDa]: %s has MW=%.3f kDa", r.species, r.mw_kda)

    cys_outliers = verts[verts["n_cys"] != 2]
    if cys_outliers.empty:
        logger.info("QC CHECK [Cys]: All %d vertebrates have exactly 2 Cys (heme-binding).", len(verts))
    else:
        for r in cys_outliers.itertuples():
            logger.warning("QC WARN [Cys != 2]: %s has n_cys=%d", r.species, r.n_cys)

    for r in qc_df[qc_df["role"] == "outgroup"].itertuples():
        logger.info(
            "QC CHECK [Outgroup length]: %s = %d aa (%+d aa vs vertebrate mode %d aa; n_cys=%d, pI=%.2f, MW=%.3f kDa)",
            r.species,
            r.length,
            r.length - vert_mode,
            vert_mode,
            r.n_cys,
            r.isoelectric_point,
            r.mw_kda,
        )


def save_figure(fig: plt.Figure, stem: Path, cfg: dict[str, Any]) -> list[Path]:
    """Save a matplotlib figure in all configured formats and close it."""
    dpi = int(cfg["reporting"]["dpi"])
    formats = cfg["reporting"]["figure_formats"]
    saved: list[Path] = []
    for fmt in formats:
        out = stem.with_suffix(f".{fmt}")
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        saved.append(out)
    plt.close(fig)
    return saved


def make_qc_figures(
    qc_df: pd.DataFrame, aa_wide_df: pd.DataFrame, fig_dir: Path, cfg: dict[str, Any], logger: Any
) -> list[Path]:
    """Generate the three Phase 2 QC figures in PNG (300 dpi) and SVG."""
    font_size = int(cfg["reporting"]["font_size"])
    plt.rcParams.update({"font.size": font_size, "axes.titlesize": font_size + 2})

    present_groups = [g for g in GROUP_ORDER if g in set(qc_df["group"])]
    palette = sns.color_palette("colorblind", n_colors=len(present_groups))
    group_colors = dict(zip(present_groups, palette))

    saved_paths: list[Path] = []

    # 1. Horizontal bar chart of sequence length by species
    fig, ax = plt.subplots(figsize=(9, 6.5))
    y_pos = np.arange(len(qc_df))
    labels = [s.replace("_", " ") for s in qc_df["species"]]
    colors = [group_colors[g] for g in qc_df["group"]]
    bars = ax.barh(y_pos, qc_df["length"], color=colors, edgecolor="black", linewidth=0.6)
    for bar, role, length in zip(bars, qc_df["role"], qc_df["length"]):
        if role == "outgroup":
            bar.set_hatch("///")
        ax.text(length + 0.3, bar.get_y() + bar.get_height() / 2, f"{length}", va="center", fontsize=font_size - 1)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlim(0, int(qc_df["length"].max()) + 8)
    ax.set_xlabel("Sequence length (amino acids)")
    ax.set_title("Cytochrome c Sequence Length Across Species (Outgroups Hatched)")
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=group_colors[g], edgecolor="black",
                      hatch="///" if g in ("Insect", "Fungus") else "")
        for g in present_groups
    ]
    ax.legend(handles, present_groups, title="Taxonomic Group", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    saved_paths.extend(save_figure(fig, fig_dir / "fig_length_by_species", cfg))

    # 2. Amino-acid composition heatmap
    mat = aa_wide_df.set_index("species")[AA_ORDER]
    mat.index = [s.replace("_", " ") for s in mat.index]
    fig, ax = plt.subplots(figsize=(12.5, 7.2))
    sns.heatmap(
        mat,
        annot=True,
        fmt=".1f",
        cmap="viridis",
        cbar_kws={"label": "Amino Acid Composition (%)"},
        linewidths=0.3,
        annot_kws={"size": max(7, font_size - 2.5)},
        ax=ax,
    )
    ax.set_xlabel("Amino Acid")
    ax.set_ylabel("Species (Biological Order)")
    ax.set_title("Cytochrome c Amino-Acid Composition (%) by Species")
    saved_paths.extend(save_figure(fig, fig_dir / "fig_aa_composition_heatmap", cfg))

    # 3. Scatter plot of pI vs length
    fig, ax = plt.subplots(figsize=(10.2, 6.4))
    for group in present_groups:
        sub = qc_df[qc_df["group"] == group]
        is_out = (sub["role"] == "outgroup").all()
        ax.scatter(
            sub["length"],
            sub["isoelectric_point"],
            label=group,
            color=group_colors[group],
            marker="^" if is_out else "o",
            s=80,
            edgecolors="black",
            linewidths=0.7,
            zorder=4,
        )

    # Clean, non-overlapping label placement
    for length_val, grp in qc_df.groupby("length"):
        grp_sorted = grp.sort_values(["isoelectric_point", "species"])
        if len(grp_sorted) == 1:
            r = grp_sorted.iloc[0]
            dx = -0.18 if length_val >= 108 else -0.16
            ha = "right"
            ax.annotate(
                r["species"].replace("_", " "),
                xy=(r["length"], r["isoelectric_point"]),
                xytext=(r["length"] + dx, r["isoelectric_point"]),
                fontsize=max(8, font_size - 1.5),
                ha=ha,
                va="center",
                arrowprops={"arrowstyle": "-", "color": "gray", "lw": 0.6, "alpha": 0.7},
            )
        else:
            # Separate outlier high-pI point (e.g., Takifugu) from the dense cluster
            high_pts = grp_sorted[grp_sorted["isoelectric_point"] > 9.70]
            cluster_pts = grp_sorted[grp_sorted["isoelectric_point"] <= 9.70]
            for r in high_pts.itertuples():
                ax.annotate(
                    r.species.replace("_", " "),
                    xy=(r.length, r.isoelectric_point),
                    xytext=(r.length + 0.25, r.isoelectric_point),
                    fontsize=max(8, font_size - 1.5),
                    ha="left",
                    va="center",
                    arrowprops={"arrowstyle": "-", "color": "gray", "lw": 0.6, "alpha": 0.7},
                )
            left_pts = list(cluster_pts.iloc[0::2].itertuples())
            right_pts = list(cluster_pts.iloc[1::2].itertuples())
            left_ys = np.linspace(9.38, 9.74, len(left_pts)) if left_pts else []
            right_ys = np.linspace(9.40, 9.72, len(right_pts)) if right_pts else []
            for r, ty in zip(left_pts, left_ys):
                ax.annotate(
                    r.species.replace("_", " "),
                    xy=(r.length, r.isoelectric_point),
                    xytext=(r.length - 0.55, float(ty)),
                    fontsize=max(8, font_size - 1.5),
                    ha="right",
                    va="center",
                    arrowprops={"arrowstyle": "-", "color": "gray", "lw": 0.5, "alpha": 0.7},
                )
            for r, ty in zip(right_pts, right_ys):
                ax.annotate(
                    r.species.replace("_", " "),
                    xy=(r.length, r.isoelectric_point),
                    xytext=(r.length + 0.55, float(ty)),
                    fontsize=max(8, font_size - 1.5),
                    ha="left",
                    va="center",
                    arrowprops={"arrowstyle": "-", "color": "gray", "lw": 0.5, "alpha": 0.7},
                )

    ax.set_xlim(int(qc_df["length"].min()) - 1.8, int(qc_df["length"].max()) + 0.8)
    ax.set_ylim(9.33, 9.85)
    ax.set_xlabel("Sequence length (amino acids)")
    ax.set_ylabel("Theoretical isoelectric point (pI)")
    ax.set_title("Cytochrome c Isoelectric Point (pI) vs. Sequence Length\n(All 17 species strongly basic: pI 9.46–9.79, well above 8.5 threshold)")
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(title="Group", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, fontsize=font_size - 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    saved_paths.extend(save_figure(fig, fig_dir / "fig_pI_vs_length", cfg))

    logger.info("Saved %d QC figure files to %s", len(saved_paths), fig_dir)
    return saved_paths


def write_processed_outputs(
    records: list[dict[str, Any]],
    ident_df: pd.DataFrame,
    proc_dir: Path,
    logger: Any,
) -> dict[str, Path]:
    """Write curated FASTAs, id_map.csv, and dataset_manifest.json."""
    curated_recs = [
        SeqRecord(Seq(r["sequence"]), id=f"{r['species']}|{r['accession']}", description="")
        for r in records
    ]
    treeid_recs = [
        SeqRecord(Seq(r["sequence"]), id=r["species"], description="")
        for r in records
    ]
    vert_recs = [
        SeqRecord(Seq(r["sequence"]), id=r["species"], description="")
        for r in records
        if r["role"] == "ingroup"
    ]

    curated_path = proc_dir / "cytc_curated.fasta"
    treeid_path = proc_dir / "cytc_curated_treeid.fasta"
    vert_path = proc_dir / "cytc_vertebrates.fasta"
    id_map_path = proc_dir / "id_map.csv"
    manifest_path = proc_dir / "dataset_manifest.json"

    write_wrapped_fasta(curated_recs, curated_path, width=60)
    write_wrapped_fasta(treeid_recs, treeid_path, width=60)
    write_wrapped_fasta(vert_recs, vert_path, width=60)

    id_map_df = pd.DataFrame(
        [
            {
                "tree_id": r["species"],
                "accession": r["accession"],
                "group": r["group"],
                "role": r["role"],
                "taxon_id": r["taxon_id"],
                "order_index": r["order_index"],
            }
            for r in records
        ]
    )
    id_map_df.to_csv(id_map_path, index=False)

    manifest = {
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "counts": {
            "total": len(records),
            "ingroup_vertebrates": len(vert_recs),
            "outgroups": len(records) - len(vert_recs),
            "identical_sequence_groups": len(ident_df),
        },
        "files": {
            "cytc_curated.fasta": {
                "path": str(curated_path.relative_to(REPO_ROOT).as_posix()),
                "records": len(curated_recs),
                "sha256": sha256_file(curated_path),
            },
            "cytc_curated_treeid.fasta": {
                "path": str(treeid_path.relative_to(REPO_ROOT).as_posix()),
                "records": len(treeid_recs),
                "sha256": sha256_file(treeid_path),
            },
            "cytc_vertebrates.fasta": {
                "path": str(vert_path.relative_to(REPO_ROOT).as_posix()),
                "records": len(vert_recs),
                "sha256": sha256_file(vert_path),
            },
            "id_map.csv": {
                "path": str(id_map_path.relative_to(REPO_ROOT).as_posix()),
                "rows": len(id_map_df),
                "sha256": sha256_file(id_map_path),
            },
        },
        "tool_versions": {
            "python": platform.python_version(),
            "biopython": Bio.__version__,
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "matplotlib": matplotlib.__version__,
            "seaborn": sns.__version__,
        },
    }
    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)

    logger.info("Wrote curated FASTAs, id_map.csv, and dataset_manifest.json to %s", proc_dir)
    return {
        "curated": curated_path,
        "treeid": treeid_path,
        "vertebrates": vert_path,
        "id_map": id_map_path,
        "manifest": manifest_path,
    }


def main() -> None:
    """Run Phase 2 sequence preparation, QC, figure generation, and curation."""
    parser = argparse.ArgumentParser(description="Phase 2: QC and curate cytochrome c sequences.")
    parser.add_argument("--config", type=str, default=None, help="Optional path to config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    logger = get_logger("20_prepare_sequences", cfg)
    logger.info(
        "Phase 2 goal: standardize, validate, QC, and curate raw cytochrome c sequences into canonical processed datasets. "
        "Inputs: data/raw/cytc_raw.fasta, results/tables/dataset_provenance.csv, config/species.csv, config/config.yaml."
    )

    raw_dir = resolve(cfg, "raw")
    proc_dir = resolve(cfg, "processed")
    tables_dir = resolve(cfg, "tables")
    fig_dir = resolve(cfg, "figures")
    species_path = resolve(cfg, "species_table")

    raw_fasta = raw_dir / "cytc_raw.fasta"
    prov_path = tables_dir / "dataset_provenance.csv"

    species_df = pd.read_csv(species_path)
    prov_df = pd.read_csv(prov_path)

    records = load_and_validate_records(raw_fasta, species_df, prov_df, cfg, logger)
    qc_df, aa_wide_df, aa_long_df, length_summary_df, ident_df = compute_qc_tables(records, logger)

    # Write QC tables
    qc_path = tables_dir / "seq_qc.csv"
    aa_wide_path = tables_dir / "aa_composition.csv"
    aa_long_path = tables_dir / "aa_composition_long.csv"
    len_sum_path = tables_dir / "length_summary.csv"
    ident_path = tables_dir / "identical_sequence_groups.csv"

    qc_df.to_csv(qc_path, index=False)
    aa_wide_df.to_csv(aa_wide_path, index=False)
    aa_long_df.to_csv(aa_long_path, index=False)
    length_summary_df.to_csv(len_sum_path, index=False)
    ident_df.to_csv(ident_path, index=False)
    logger.info("Wrote QC tables to %s", tables_dir)

    # Generate QC figures
    make_qc_figures(qc_df, aa_wide_df, fig_dir, cfg, logger)

    # Write processed curated dataset
    out_paths = write_processed_outputs(records, ident_df, proc_dir, logger)

    n_vert = sum(1 for r in records if r["role"] == "ingroup")
    n_out = len(records) - n_vert
    logger.info(
        "SUMMARY: curated %d sequences (%d vertebrates, %d outgroups, %d identical-sequence groups); "
        "curated=%s, treeid=%s, vertebrates=%s, qc=%s",
        len(records),
        n_vert,
        n_out,
        len(ident_df),
        out_paths["curated"],
        out_paths["treeid"],
        out_paths["vertebrates"],
        qc_path,
    )


if __name__ == "__main__":
    main()

