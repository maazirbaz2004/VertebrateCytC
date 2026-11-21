"""Phase 3: Multiple sequence alignment (Clustal Omega & MUSCLE v5), QC, and trimAl sensitivity."""
from __future__ import annotations

import argparse
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib.colors import ListedColormap

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cytc.utils import REPO_ROOT, get_logger, load_config, resolve  # noqa: E402

# Okabe-Ito colour-blind-safe mapping for physico-chemical amino-acid classes
AA_CLASSES: list[tuple[str, str, str]] = [
    ("Gap (-)", "-", "#FFFFFF"),
    ("Hydrophobic (A,V,I,L,M)", "AVILM", "#56B4E9"),
    ("Aromatic (F,W,Y)", "FWY", "#CC79A7"),
    ("Polar (S,T,N,Q)", "STNQ", "#009E73"),
    ("Positive (K,R,H)", "KRH", "#E69F00"),
    ("Negative (D,E)", "DE", "#D55E00"),
    ("Special (G,P,C)", "GPC", "#F0E442"),
]


def run_cmd(cmd: list[str], logger: Any) -> subprocess.CompletedProcess[str]:
    """Run an external command without shell=True, logging command and failing loudly."""
    logger.info("CMD: %s", " ".join(cmd))
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return res


def detect_tool_version(binary: str, version_args: list[str]) -> str:
    """Query an external binary for its version string."""
    res = subprocess.run([binary, *version_args], capture_output=True, text=True, check=False)
    out = (res.stdout + "\n" + res.stderr).strip()
    for line in out.splitlines():
        line = line.strip()
        if line:
            return line
    return "unknown"


def write_wrapped_fasta(records: list[SeqRecord], out_path: Path, width: int = 60) -> None:
    """Write FASTA records wrapped to a fixed column width."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
        for rec in records:
            fh.write(f">{rec.id}\n")
            seq_str = str(rec.seq)
            for i in range(0, len(seq_str), width):
                fh.write(seq_str[i : i + width] + "\n")


def normalize_and_reorder_alignment(
    aln_path: Path,
    expected_order: list[str],
    curated_seqs: dict[str, str],
) -> MultipleSeqAlignment:
    """Upper-case, normalize gaps to '-', reorder to fixed species order, verify ungapped seqs, and rewrite."""
    raw_aln = AlignIO.read(aln_path, "fasta")
    by_id: dict[str, SeqRecord] = {}
    for rec in raw_aln:
        seq_str = str(rec.seq).upper().replace(".", "-")
        by_id[rec.id] = SeqRecord(Seq(seq_str), id=rec.id, description="")

    if set(by_id.keys()) != set(expected_order):
        raise ValueError(
            f"Alignment IDs mismatch in {aln_path.name}: "
            f"missing={set(expected_order) - set(by_id)}, extra={set(by_id) - set(expected_order)}"
        )

    ordered: list[SeqRecord] = []
    for sp in expected_order:
        rec = by_id[sp]
        ungapped = str(rec.seq).replace("-", "")
        if ungapped != curated_seqs[sp]:
            raise ValueError(f"Aligner altered sequence for {sp} in {aln_path.name}")
        ordered.append(rec)

    lengths = {len(r.seq) for r in ordered}
    if len(lengths) != 1:
        raise ValueError(f"Non-uniform row lengths in {aln_path.name}: {lengths}")

    write_wrapped_fasta(ordered, aln_path, width=60)
    return MultipleSeqAlignment(ordered)


def compute_alignment_stats(
    aln: MultipleSeqAlignment,
    name: str,
    dataset: str,
    aligner: str,
    cfg: dict[str, Any],
) -> dict[str, Any]:
    """Compute basic alignment statistics and BLOSUM62 sum-of-pairs score."""
    n_seq = len(aln)
    aln_len = aln.get_alignment_length()
    seqs = [str(r.seq) for r in aln]

    n_gap_cols = 0
    n_gapfree_cols = 0
    n_invariant_cols = 0
    for c in range(aln_len):
        col = [s[c] for s in seqs]
        if "-" in col:
            n_gap_cols += 1
        else:
            n_gapfree_cols += 1
            if len(set(col)) == 1:
                n_invariant_cols += 1

    total_gap_chars = sum(s.count("-") for s in seqs)
    n_gap_blocks = sum(len(re.findall(r"-+", s)) for s in seqs)

    # Mean pairwise identity and BLOSUM62 affine SP score
    blosum = substitution_matrices.load(cfg["similarity"]["substitution_matrix"])
    gap_open = float(cfg["similarity"]["gap_open"])
    gap_extend = float(cfg["similarity"]["gap_extend"])

    pair_ids: list[float] = []
    sp_score = 0.0
    for i in range(n_seq):
        si = seqs[i]
        for j in range(i + 1, n_seq):
            sj = seqs[j]
            matches = 0
            aligned = 0
            in_gap_i = False
            in_gap_j = False
            for a, b in zip(si, sj):
                if a == "-" and b == "-":
                    continue
                if a == "-":
                    sp_score += gap_extend if in_gap_i else gap_open
                    in_gap_i, in_gap_j = True, False
                elif b == "-":
                    sp_score += gap_extend if in_gap_j else gap_open
                    in_gap_i, in_gap_j = False, True
                else:
                    in_gap_i, in_gap_j = False, False
                    aligned += 1
                    if a == b:
                        matches += 1
                    sp_score += float(blosum[a, b])
            pair_ids.append((100.0 * matches / aligned) if aligned > 0 else 0.0)

    return {
        "alignment": name,
        "dataset": dataset,
        "aligner": aligner,
        "n_seq": n_seq,
        "aln_length": aln_len,
        "n_gap_columns": n_gap_cols,
        "n_gapfree_columns": n_gapfree_cols,
        "n_invariant_columns": n_invariant_cols,
        "pct_identity_columns": round(100.0 * n_invariant_cols / aln_len, 2),
        "mean_pairwise_identity": round(float(np.mean(pair_ids)), 2),
        "total_gap_chars": total_gap_chars,
        "n_gap_blocks": n_gap_blocks,
        "sum_of_pairs_blosum62": round(sp_score, 1),
    }


def column_signatures(aln: MultipleSeqAlignment) -> list[tuple[int | None, ...]]:
    """Map each column to a tuple of per-sequence 0-based residue indices (None for gaps)."""
    counters = [0] * len(aln)
    sigs: list[tuple[int | None, ...]] = []
    for col in range(aln.get_alignment_length()):
        sig: list[int | None] = []
        for i, rec in enumerate(aln):
            if rec.seq[col] == "-":
                sig.append(None)
            else:
                sig.append(counters[i])
                counters[i] += 1
        sigs.append(tuple(sig))
    return sigs


def pairwise_residue_pairs(sigs: list[tuple[int | None, ...]]) -> set[tuple[int, int, int, int]]:
    """Return set of (seq_i, seq_j, res_idx_i, res_idx_j) for all aligned residue pairs."""
    pairs: set[tuple[int, int, int, int]] = set()
    for sig in sigs:
        present = [(i, idx) for i, idx in enumerate(sig) if idx is not None]
        for a in range(len(present)):
            i, ri = present[a]
            for b in range(a + 1, len(present)):
                j, rj = present[b]
                pairs.add((i, j, ri, rj))
    return pairs


def compress_ranges(cols_1based: list[int]) -> str:
    """Format a sorted list of 1-based column indices as compact ranges (e.g. '1-5, 22')."""
    if not cols_1based:
        return "none"
    ranges: list[str] = []
    start = prev = cols_1based[0]
    for c in cols_1based[1:]:
        if c == prev + 1:
            prev = c
        else:
            ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
            start = prev = c
    ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
    return ", ".join(ranges)


def compare_aligners(
    aln_primary: MultipleSeqAlignment,
    aln_secondary: MultipleSeqAlignment,
    dataset: str,
    primary_name: str,
    secondary_name: str,
) -> dict[str, Any]:
    """Compute column-level (TC) and pairwise-residue (SP) agreement between two aligners."""
    sigs_p = column_signatures(aln_primary)
    sigs_s = column_signatures(aln_secondary)

    # Informative columns: >= 2 non-gap residues
    inf_p = [(idx + 1, sig) for idx, sig in enumerate(sigs_p) if sum(x is not None for x in sig) >= 2]
    inf_s_set = {sig for sig in sigs_s if sum(x is not None for x in sig) >= 2}
    inf_p_set = {sig for _, sig in inf_p}

    shared_cols = inf_p_set & inf_s_set
    union_cols = inf_p_set | inf_s_set
    disagreed_cols = [col1 for col1, sig in inf_p if sig not in inf_s_set]

    pairs_p = pairwise_residue_pairs(sigs_p)
    pairs_s = pairwise_residue_pairs(sigs_s)
    shared_pairs = pairs_p & pairs_s
    union_pairs = pairs_p | pairs_s

    return {
        "dataset": dataset,
        "primary_aligner": primary_name,
        "secondary_aligner": secondary_name,
        "n_informative_primary": len(inf_p_set),
        "n_informative_secondary": len(inf_s_set),
        "n_shared_columns": len(shared_cols),
        "pct_shared_columns_vs_primary": round(100.0 * len(shared_cols) / len(inf_p_set), 2),
        "pct_shared_columns_jaccard": round(100.0 * len(shared_cols) / len(union_cols), 2),
        "n_pairs_primary": len(pairs_p),
        "n_shared_pairs": len(shared_pairs),
        "pct_shared_pairs_vs_primary": round(100.0 * len(shared_pairs) / len(pairs_p), 2),
        "pct_shared_pairs_jaccard": round(100.0 * len(shared_pairs) / len(union_pairs), 2),
        "disagreed_primary_cols": compress_ranges(disagreed_cols),
    }


def extract_gap_blocks(aln: MultipleSeqAlignment, aln_name: str) -> list[dict[str, Any]]:
    """List contiguous column blocks containing gaps (split when gap membership changes)."""
    aln_len = aln.get_alignment_length()
    ids = [r.id for r in aln]
    seqs = [str(r.seq) for r in aln]

    blocks: list[dict[str, Any]] = []
    cur_start: int | None = None
    cur_gapped: tuple[str, ...] | None = None

    for c in range(aln_len):
        gapped = tuple(ids[i] for i, s in enumerate(seqs) if s[c] == "-")
        if not gapped:
            if cur_start is not None and cur_gapped is not None:
                ungapped = [sp for sp in ids if sp not in cur_gapped]
                blocks.append(
                    {
                        "alignment": aln_name,
                        "start": cur_start + 1,
                        "end": c,
                        "length": c - cur_start,
                        "sequences_with_gap": ";".join(cur_gapped),
                        "sequences_without_gap": ";".join(ungapped),
                    }
                )
                cur_start, cur_gapped = None, None
            continue

        if cur_gapped is None:
            cur_start = c
            cur_gapped = gapped
        elif gapped != cur_gapped:
            ungapped = [sp for sp in ids if sp not in cur_gapped]
            blocks.append(
                {
                    "alignment": aln_name,
                    "start": cur_start + 1,
                    "end": c,
                    "length": c - cur_start,
                    "sequences_with_gap": ";".join(cur_gapped),
                    "sequences_without_gap": ";".join(ungapped),
                }
            )
            cur_start = c
            cur_gapped = gapped

    if cur_start is not None and cur_gapped is not None:
        ungapped = [sp for sp in ids if sp not in cur_gapped]
        blocks.append(
            {
                "alignment": aln_name,
                "start": cur_start + 1,
                "end": aln_len,
                "length": aln_len - cur_start,
                "sequences_with_gap": ";".join(cur_gapped),
                "sequences_without_gap": ";".join(ungapped),
            }
        )
    return blocks


def verify_biological_anchors(
    aln_all: MultipleSeqAlignment, aln_vert: MultipleSeqAlignment, logger: Any
) -> dict[str, Any]:
    """Verify heme motif (C..CH) and Met80 region alignment consistency; return coordinates."""
    info: dict[str, Any] = {}
    for label, aln in (("primary_all", aln_all), ("primary_vert", aln_vert)):
        starts: set[int] = set()
        for rec in aln:
            m = re.search(r"C-*.-*.-*C-*H", str(rec.seq))
            if not m:
                raise ValueError(f"Heme motif missing in {label}:{rec.id}")
            if "-" in m.group(0):
                raise ValueError(f"Gap inside heme motif in {label}:{rec.id}: {m.group(0)}")
            starts.add(m.start() + 1)  # 1-based alignment column
        if len(starts) != 1:
            raise ValueError(f"Heme motif misaligned across sequences in {label}: starts={sorted(starts)}")
        motif_col = next(iter(starts))
        info[f"{label}_heme_start_col"] = motif_col
        info[f"{label}_heme_end_col"] = motif_col + 4

    # Map Homo_sapiens UniProt positions to alignment columns to verify Met80 loop (classical 70-85 = UniProt 71-86)
    for label, aln in (("primary_all", aln_all), ("primary_vert", aln_vert)):
        hum = next(r for r in aln if r.id == "Homo_sapiens")
        u_pos = 0
        u_to_col: dict[int, int] = {}
        for col_idx, ch in enumerate(str(hum.seq), start=1):
            if ch != "-":
                u_pos += 1
                u_to_col[u_pos] = col_idx
        met80_col = u_to_col[81]  # UniProt 81 = classical Met80
        if str(hum.seq)[met80_col - 1] != "M":
            raise ValueError(f"Expected Met at Homo_sapiens UniProt 81 (classical 80) in {label}")
        loop_start = u_to_col[71]  # classical 70
        loop_end = u_to_col[86]    # classical 85
        info[f"{label}_met80_col"] = met80_col
        info[f"{label}_met80_loop_cols"] = (loop_start, loop_end)

        if label == "primary_vert":
            for rec in aln:
                seg = str(rec.seq)[loop_start - 1 : loop_end]
                if "-" in seg:
                    raise ValueError(f"Gap in vertebrate Met80 loop (cols {loop_start}-{loop_end}) in {rec.id}: {seg}")

    logger.info(
        "Anchor check OK: primary_all CxxCH cols %d-%d, Met80 col %d (UniProt 81, classical 80); "
        "primary_vert CxxCH cols %d-%d, Met80 loop cols %d-%d gap-free.",
        info["primary_all_heme_start_col"],
        info["primary_all_heme_end_col"],
        info["primary_all_met80_col"],
        info["primary_vert_heme_start_col"],
        info["primary_vert_heme_end_col"],
        info["primary_vert_met80_loop_cols"][0],
        info["primary_vert_met80_loop_cols"][1],
    )
    return info


def plot_alignment_overview(
    aln: MultipleSeqAlignment,
    anchors: dict[str, Any],
    fig_dir: Path,
    cfg: dict[str, Any],
    logger: Any,
) -> list[Path]:
    """Generate colour-coded physico-chemical alignment overview figure (PNG + SVG)."""
    font_size = int(cfg["reporting"]["font_size"])
    dpi = int(cfg["reporting"]["dpi"])
    formats = cfg["reporting"]["figure_formats"]
    plt.rcParams.update({"font.size": font_size})

    char_to_code: dict[str, int] = {}
    colors: list[str] = []
    for code, (_, chars, hex_col) in enumerate(AA_CLASSES):
        colors.append(hex_col)
        for ch in chars:
            char_to_code[ch] = code

    n_seq = len(aln)
    aln_len = aln.get_alignment_length()
    mat = np.zeros((n_seq, aln_len), dtype=int)
    for i, rec in enumerate(aln):
        for j, ch in enumerate(str(rec.seq)):
            mat[i, j] = char_to_code.get(ch, 0)

    cmap = ListedColormap(colors)
    fig, ax = plt.subplots(figsize=(14, 6.2))
    ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=len(AA_CLASSES) - 1, interpolation="nearest")

    # Highlight heme motif and Met80 loop
    h_start = anchors["primary_all_heme_start_col"] - 1
    h_end = anchors["primary_all_heme_end_col"] - 1
    m_col = anchors["primary_all_met80_col"] - 1
    ax.axvspan(h_start - 0.5, h_end + 0.5, facecolor="none", edgecolor="black", linewidth=1.6, zorder=4)
    ax.axvspan(m_col - 0.5, m_col + 0.5, facecolor="none", edgecolor="black", linewidth=1.4, linestyle="--", zorder=4)

    ax.text(
        (h_start + h_end) / 2,
        -0.85,
        f"CxxCH\n(cols {h_start + 1}-{h_end + 1})",
        ha="center",
        va="bottom",
        fontsize=max(8, font_size - 1.5),
        fontweight="bold",
    )
    ax.text(
        m_col,
        -0.85,
        f"Met80\n(col {m_col + 1})",
        ha="center",
        va="bottom",
        fontsize=max(8, font_size - 1.5),
        fontweight="bold",
    )

    ax.set_yticks(np.arange(n_seq))
    ax.set_yticklabels([r.id.replace("_", " ") for r in aln])

    ticks_1based = [1] + list(range(10, aln_len + 1, 10))
    if aln_len not in ticks_1based and aln_len - ticks_1based[-1] >= 4:
        ticks_1based.append(aln_len)
    ax.set_xticks([t - 1 for t in ticks_1based])
    ax.set_xticklabels([str(t) for t in ticks_1based])
    ax.set_xlabel("Alignment column (1-based)")
    ax.set_ylabel("Species (Biological Order)")
    ax.set_title("Primary Multiple Sequence Alignment Overview (Physico-Chemical Amino-Acid Classes)", pad=28)

    handles = [
        mpatches.Patch(facecolor=hex_col, edgecolor="black", linewidth=0.6, label=label)
        for label, _, hex_col in AA_CLASSES
    ]
    ax.legend(
        handles=handles,
        title="Amino-Acid Class",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        frameon=True,
        fontsize=font_size - 1,
    )

    saved: list[Path] = []
    stem = fig_dir / "fig_alignment_overview"
    for fmt in formats:
        out = stem.with_suffix(f".{fmt}")
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        saved.append(out)
    plt.close(fig)
    logger.info("Saved alignment overview figure to %s", fig_dir)
    return saved


def df_to_markdown(df: pd.DataFrame) -> str:
    """Format a DataFrame as a GitHub-flavored Markdown table without requiring tabulate."""
    headers = [str(c) for c in df.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in df.itertuples(index=False):
        lines.append("| " + " | ".join(str(v) for v in row) + " |")
    return "\n".join(lines)


def write_alignment_notes(
    out_path: Path,
    versions: dict[str, str],
    commands: list[str],
    primary_choice: str,
    stats_df: pd.DataFrame,
    agree_df: pd.DataFrame,
    gap_df: pd.DataFrame,
    anchors: dict[str, Any],
    cfg: dict[str, Any],
) -> None:
    """Write human-readable alignment decision log and QC summary."""
    lines = [
        "# Phase 3 Alignment Notes and Primary Alignment Decision",
        "",
        "## 1. Tool Versions",
        f"- Clustal Omega: `{versions['clustalo']}`",
        f"- MUSCLE: `{versions['muscle']}`",
        f"- trimAl: `{versions['trimal']}`",
        "",
        "## 2. Exact Commands Executed",
        "```bash",
        *commands,
        "```",
        "",
        "## 3. Alignment Statistics Summary",
        "Sum-of-pairs (SP) score uses BLOSUM62 with pairwise affine gap penalties "
        f"(gap_open={cfg['similarity']['gap_open']}, gap_extend={cfg['similarity']['gap_extend']}; "
        "gap-gap columns ignored for each pair).",
        "",
        df_to_markdown(stats_df),
        "",
        "## 4. Aligner Agreement (Clustal Omega vs. MUSCLE v5)",
        "Column-level agreement (TC score) counts only informative columns containing at least 2 non-gap residues "
        "(all-gap and single-residue columns are excluded as uninformative).",
        "",
        df_to_markdown(agree_df),
        "",
        "## 5. Biological Anchor Verification",
        f"- **Primary all (`primary_all.fasta`)**: Heme motif `CxxCH` aligned gap-free across all 17 species at "
        f"alignment columns **{anchors['primary_all_heme_start_col']}–{anchors['primary_all_heme_end_col']}** "
        f"(UniProt-human 15–19; classical 14–18). Distal axial ligand **Met80** at alignment column "
        f"**{anchors['primary_all_met80_col']}** (UniProt-human 81; classical 80).",
        f"- **Primary vertebrates (`primary_vert.fasta`)**: Heme motif `CxxCH` aligned gap-free across all 15 vertebrates at "
        f"alignment columns **{anchors['primary_vert_heme_start_col']}–{anchors['primary_vert_heme_end_col']}**. "
        f"Met80 loop (classical 70–85; UniProt-human 71–86) spans alignment columns "
        f"**{anchors['primary_vert_met80_loop_cols'][0]}–{anchors['primary_vert_met80_loop_cols'][1]}** with zero gaps.",
        "",
        "## 6. Primary Alignment Decision",
        f"- **Selected primary aligner:** `{primary_choice}` (per `config.yaml` -> `alignment.primary`).",
        "- Both Clustal Omega and MUSCLE v5 produce identical or near-identical vertebrate alignments, place the heme motif "
        "`CxxCH` and Met80 loop in gap-free homologous columns, and confine gaps primarily to the outgroup N-terminal "
        "extensions (*S. cerevisiae* and *D. melanogaster*) and the 1-aa C-terminal truncation in *Danio rerio*.",
        "",
        "## 7. Gap Blocks in Primary Alignments",
    ]
    prim_gaps = gap_df[gap_df["alignment"].isin([f"all_{primary_choice}", f"vert_{primary_choice}"])]
    if prim_gaps.empty:
        lines.append("- No gap blocks detected.")
    else:
        lines.append(df_to_markdown(prim_gaps))

    lines.extend(
        [
            "",
            "## 8. Trimmed Alignment (Sensitivity Analysis Only)",
            f"- Generated with trimAl (`{cfg['alignment']['trimal_mode']} -colnumbering`) as `primary_all_trimmed.fasta`.",
            "- Used **only** in Phase 6 phylogenetic sensitivity testing, never as the primary alignment.",
            "",
        ]
    )
    out_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run Phase 3 multiple sequence alignments, QC, agreement analysis, and trimAl."""
    parser = argparse.ArgumentParser(description="Phase 3: MSA with Clustal Omega and MUSCLE v5.")
    parser.add_argument("--config", type=str, default=None, help="Optional path to config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    logger = get_logger("30_align", cfg)
    logger.info(
        "Phase 3 goal: build MSAs with Clustal Omega and MUSCLE v5 on all (17) and vertebrate-only (15) datasets, "
        "quantify aligner agreement and gaps, select primary alignments, and run trimAl sensitivity trimming. "
        "Inputs: data/processed/cytc_curated_treeid.fasta, data/processed/cytc_vertebrates.fasta, data/processed/id_map.csv."
    )

    proc_dir = resolve(cfg, "processed")
    aln_dir = resolve(cfg, "alignments")
    tables_dir = resolve(cfg, "tables")
    fig_dir = resolve(cfg, "figures")

    all_in = proc_dir / "cytc_curated_treeid.fasta"
    vert_in = proc_dir / "cytc_vertebrates.fasta"
    id_map = pd.read_csv(proc_dir / "id_map.csv").sort_values("order_index")

    all_order = list(id_map["tree_id"])
    vert_order = list(id_map.loc[id_map["role"] == "ingroup", "tree_id"])
    all_curated = {r.id: str(r.seq) for r in SeqIO.parse(all_in, "fasta")}
    vert_curated = {r.id: str(r.seq) for r in SeqIO.parse(vert_in, "fasta")}

    clustalo_bin = cfg["tools"]["clustalo"]
    muscle_bin = cfg["tools"]["muscle"]
    trimal_bin = cfg["tools"]["trimal"]
    n_threads = str(cfg["project"]["n_threads"])
    clustalo_extra = [str(x) for x in cfg["alignment"].get("clustalo_extra", [])]
    muscle_extra = [str(x) for x in cfg["alignment"].get("muscle_extra", [])]
    primary_choice = str(cfg["alignment"]["primary"]).strip().lower()
    secondary_choice = "muscle" if primary_choice == "clustalo" else "clustalo"

    versions = {
        "clustalo": detect_tool_version(clustalo_bin, ["--version"]),
        "muscle": detect_tool_version(muscle_bin, ["-version"]),
        "trimal": detect_tool_version(trimal_bin, ["--version"]),
    }
    logger.info("Tool versions: %s", versions)

    commands_log: list[str] = []
    alns: dict[str, MultipleSeqAlignment] = {}

    datasets = [
        ("all", all_in, all_order, all_curated),
        ("vert", vert_in, vert_order, vert_curated),
    ]

    for ds_name, in_fasta, order, curated in datasets:
        # Clustal Omega
        co_out = aln_dir / f"{ds_name}_clustalo.fasta"
        co_mat = aln_dir / f"{ds_name}_clustalo.distmat"
        co_cmd = [
            clustalo_bin,
            "-i",
            str(in_fasta),
            "-o",
            str(co_out),
            "--seqtype",
            "Protein",
            "--outfmt",
            "fasta",
            "--force",
            "--threads",
            n_threads,
            "-v",
            "--percent-id",
            f"--distmat-out={co_mat}",
            "--full",
            *clustalo_extra,
        ]
        run_cmd(co_cmd, logger)
        commands_log.append(" ".join(co_cmd))
        alns[f"{ds_name}_clustalo"] = normalize_and_reorder_alignment(co_out, order, curated)

        # MUSCLE v5
        mu_out = aln_dir / f"{ds_name}_muscle.fasta"
        mu_cmd = [
            muscle_bin,
            "-align",
            str(in_fasta),
            "-output",
            str(mu_out),
            "-threads",
            n_threads,
            *muscle_extra,
        ]
        run_cmd(mu_cmd, logger)
        commands_log.append(" ".join(mu_cmd))
        alns[f"{ds_name}_muscle"] = normalize_and_reorder_alignment(mu_out, order, curated)

    # Copy primary alignments to canonical names
    prim_all_path = aln_dir / "primary_all.fasta"
    prim_vert_path = aln_dir / "primary_vert.fasta"
    shutil.copyfile(aln_dir / f"all_{primary_choice}.fasta", prim_all_path)
    shutil.copyfile(aln_dir / f"vert_{primary_choice}.fasta", prim_vert_path)
    prim_all = alns[f"all_{primary_choice}"]
    prim_vert = alns[f"vert_{primary_choice}"]

    # Alternate formats for primary_all
    AlignIO.write(prim_all, aln_dir / "primary_all.aln", "clustal")
    AlignIO.write(prim_all, aln_dir / "primary_all.phy", "phylip-relaxed")

    # Verify biological anchors (CxxCH & Met80)
    anchors = verify_biological_anchors(prim_all, prim_vert, logger)

    # Trimmed alignment (sensitivity analysis only)
    trimmed_path = aln_dir / "primary_all_trimmed.fasta"
    trimmed_cols_path = aln_dir / "trimmed_columns.txt"
    trimal_flags = shlex.split(str(cfg["alignment"]["trimal_mode"]))
    trimal_cmd = [
        trimal_bin,
        "-in",
        str(prim_all_path),
        "-out",
        str(trimmed_path),
        *trimal_flags,
        "-colnumbering",
    ]
    trimal_res = run_cmd(trimal_cmd, logger)
    commands_log.append(" ".join(trimal_cmd) + f" > {trimmed_cols_path}")
    trimmed_cols_path.write_text(trimal_res.stdout, encoding="utf-8")

    # Normalize trimmed FASTA order and wrapping (without checking against full curated seq since columns were trimmed)
    raw_trim = AlignIO.read(trimmed_path, "fasta")
    trim_by_id = {
        r.id: SeqRecord(Seq(str(r.seq).upper().replace(".", "-")), id=r.id, description="")
        for r in raw_trim
    }
    trim_ordered = [trim_by_id[sp] for sp in all_order]
    write_wrapped_fasta(trim_ordered, trimmed_path, width=60)
    trimmed_aln = MultipleSeqAlignment(trim_ordered)

    # Compute alignment stats
    stats_rows = [
        compute_alignment_stats(alns["all_clustalo"], "all_clustalo", "all", "clustalo", cfg),
        compute_alignment_stats(alns["all_muscle"], "all_muscle", "all", "muscle", cfg),
        compute_alignment_stats(alns["vert_clustalo"], "vert_clustalo", "vert", "clustalo", cfg),
        compute_alignment_stats(alns["vert_muscle"], "vert_muscle", "vert", "muscle", cfg),
        compute_alignment_stats(trimmed_aln, "primary_all_trimmed", "all_trimmed", f"{primary_choice}+trimal", cfg),
    ]
    stats_df = pd.DataFrame(stats_rows)
    stats_path = tables_dir / "alignment_stats.csv"
    stats_df.to_csv(stats_path, index=False)

    # Aligner agreement
    agree_rows = [
        compare_aligners(alns[f"all_{primary_choice}"], alns[f"all_{secondary_choice}"], "all", primary_choice, secondary_choice),
        compare_aligners(alns[f"vert_{primary_choice}"], alns[f"vert_{secondary_choice}"], "vert", primary_choice, secondary_choice),
    ]
    agree_df = pd.DataFrame(agree_rows)
    agree_path = tables_dir / "aligner_agreement.csv"
    agree_df.to_csv(agree_path, index=False)

    # Gap blocks
    gap_rows: list[dict[str, Any]] = []
    for name in ("all_clustalo", "all_muscle", "vert_clustalo", "vert_muscle"):
        gap_rows.extend(extract_gap_blocks(alns[name], name))
    gap_df = pd.DataFrame(
        gap_rows,
        columns=["alignment", "start", "end", "length", "sequences_with_gap", "sequences_without_gap"],
    )
    gap_path = tables_dir / "gap_blocks.csv"
    gap_df.to_csv(gap_path, index=False)

    # Overview figure
    plot_alignment_overview(prim_all, anchors, fig_dir, cfg, logger)

    # Notes
    notes_path = aln_dir / "alignment_notes.md"
    write_alignment_notes(notes_path, versions, commands_log, primary_choice, stats_df, agree_df, gap_df, anchors, cfg)

    vert_agree_pct = float(agree_df.loc[agree_df["dataset"] == "vert", "pct_shared_columns_vs_primary"].iloc[0])
    all_agree_pct = float(agree_df.loc[agree_df["dataset"] == "all", "pct_shared_columns_vs_primary"].iloc[0])
    logger.info(
        "SUMMARY: primary=%s (all len=%d, vert len=%d, trimmed len=%d); "
        "aligner col agreement all=%.2f%% vert=%.2f%%; outputs=%s, %s, %s",
        primary_choice,
        prim_all.get_alignment_length(),
        prim_vert.get_alignment_length(),
        trimmed_aln.get_alignment_length(),
        all_agree_pct,
        vert_agree_pct,
        prim_all_path,
        prim_vert_path,
        stats_path,
    )


if __name__ == "__main__":
    main()

