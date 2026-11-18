# PHASE 4: Conservation Analysis

**Read `CLAUDE.md` first. Prerequisite: Phase 3 is done.**

## 0. Objective

Quantify, for every alignment column, how conserved cytochrome c is; classify positions (invariant / highly conserved / moderately conserved / variable); detect contiguous conserved regions; map positions to UniProt-human and classical numbering; and produce conservation plots and a sequence logo.

**Serves:** Objective 3; **RQ1** (how conserved?) and **RQ2** (which positions/regions?).

## 1. Prerequisites

- `results/alignments/primary_all.fasta` and `primary_vert.fasta` pass Phase 3 tests.
- `data/processed/id_map.csv` exists.

## 2. Inputs

- `results/alignments/primary_vert.fasta` (**primary basis for conservation**, vertebrates only; this avoids the outgroup's N-terminal extension inflating gaps).
- `results/alignments/primary_all.fasta` (used for a second, "deep conservation" view including Drosophila and yeast).
- `config.yaml` -> `conservation` thresholds.

## 3. Tasks

### Task 4.1: Coordinate mapping (build this first, since everything else depends on it)

Create `results/tables/position_map.csv`, one row per alignment column of `primary_all.fasta`:

| Column | Meaning |
|--------|---------|
| `aln_col` | 1-based alignment column |
| `human_uniprot_pos` | position in *Homo sapiens* sequence (Met1 = 1), blank if human has a gap there |
| `classical_pos` | `human_uniprot_pos - 1` (classical numbering with Met removed), blank for Met1 or human gaps |
| `human_residue` | residue in human |
| `yeast_pos` | position in the yeast sequence (UniProt numbering), blank if gap |

Verify the convention: the human heme-binding Cys residues should map to classical positions 14 and 17 and His to 18, and in UniProt numbering to 15, 18 and 19. If not, log an error and investigate before continuing (this validates the whole mapping).

### Task 4.2: Per-column metrics (`scripts/40_conservation.py`)

Compute for **each alignment** (vertebrates-only and all) and write to `results/conservation/per_position_{vert,all}.csv`:

| Metric | Definition |
|--------|-----------|
| `n_seqs`, `n_gaps`, `gap_fraction` | counts per column |
| `consensus` | most frequent residue (ignoring gaps); `-` if column is gap-dominated (> `max_gap_fraction_per_column`) |
| `identity_frac` | frequency of the consensus residue **among non-gap sequences** |
| `n_distinct_residues` | number of distinct residues among non-gap sequences |
| `shannon_entropy_bits` | `H = -sum p_i * log2(p_i)` over non-gap residues (0 = invariant, max about 4.32 bits) |
| `norm_entropy` | `H / log2(20)` |
| `conservation_score` | `1 - norm_entropy` |
| `jsd_score` | Jensen-Shannon divergence to the BLOSUM62 background distribution (Capra and Singh 2007, window = 0, gap penalty = multiply by `(1 - gap_fraction)`). Implement from the published definition; document the background frequencies used. |
| `physchem_conservation` | fraction of non-gap residues in the consensus residue's physico-chemical class (classes: aliphatic AVLIM; aromatic FWY; positive KRH; negative DE; polar STNQ; special G, P, C) |
| `class` | `invariant` (identity = 1.0, gap-free), `highly_conserved` (>= `high_threshold`), `moderately_conserved` (>= `moderate_threshold`), `variable` (below), `gap_rich` |

Implementation notes:
- Treat the sequence set as unweighted (n is small and the tree is not star-like, but **state this**). Optionally also compute a sequence-weighted version (Henikoff weights) as a sensitivity check and log whether classification changes.
- With only 15 vertebrates and 17 sequences overall, entropy estimates are coarse. **Do not apply pseudo-count smoothing silently.** If you add one, document it.
- Smoothing for plots only: moving average of `conservation_score` with window `smoothing_window` (centered, `min_periods=1`). **Classification uses unsmoothed values.**

### Task 4.3: Region detection

Detect **conserved regions** and **variable regions** from the vertebrate alignment:
- A conserved region = a maximal run of >= `min_region_length` consecutive columns that are `invariant` or `highly_conserved`.
- A variable region = a maximal run of >= 2 consecutive `variable` columns (list them even if short).
- Output `results/conservation/regions.csv`: `region_type, aln_start, aln_end, human_uniprot_start, human_uniprot_end, classical_start, classical_end, length, mean_conservation, sequence_in_human`.
- **Merge regions separated by a single non-conserved column** into one "conserved block with a gap" only if you document the rule; otherwise keep strict.

### Task 4.4: Summary statistics (these directly answer RQ1)

Write `results/tables/conservation_summary.csv` and a short text block into `results/conservation/summary.md`:

For vertebrates only and all-17:
- total columns; gap-free columns;
- number and percentage of **invariant** columns; of **>= 90% conserved**; of moderately conserved; of variable;
- mean and median `identity_frac`;
- mean entropy;
- the **N most variable positions** (top 10) with their residue distributions (for example `K:9, R:4, Q:2`);
- number of **substitution-free** stretches (the longest invariant stretch and its location);
- proportion of vertebrate invariant sites that remain invariant when the two outgroups are added (**deep conservation**). Report and list the positions invariant across all 17 sequences.

Also compute **how many positions differ between each vertebrate and human** (full table appears in Phase 5; here only a quick log).

### Task 4.5: Figures

1. **`fig_conservation_profile.png/.svg`**: line plot of `conservation_score` (raw as thin line, smoothed as thick line) over **human UniProt position** (x-axis); shaded bands for detected conserved regions; markers on the heme motif (Cys15/Cys18/His19 UniProt) and Met81 (UniProt numbering of classical Met80); secondary axis or annotation for classical numbering. Use vertebrate alignment. Panel B (below): same for all-17 alignment.
2. **`fig_entropy_bar.png/.svg`**: per-column Shannon entropy as bars, coloured by class.
3. **`fig_sequence_logo_vert.png/.svg`** and **`fig_sequence_logo_all.png/.svg`**: logomaker information-content logos (`logomaker.alignment_to_matrix(seqs, to_type="information", characters_to_ignore="-")`), split into 3 stacked rows of about 35-36 columns each so letters are legible, with x tick labels as human UniProt positions. Colour by chemistry (`color_scheme="chemistry"`).
4. **`fig_conservation_class_pie_or_bar.png/.svg`**: counts of positions per class (vertebrate vs. all).
5. **`fig_conservation_heatmap_per_species.png/.svg`**: heatmap of `species x column` showing whether each residue **matches the vertebrate consensus** (binary), so lineage-specific substitutions stand out; rows in fixed order, columns restricted to the variable columns only for legibility, with column labels as human UniProt positions.

### Task 4.6: Lineage-specific substitution table

Create `results/tables/lineage_specific_substitutions.csv`: for each non-invariant vertebrate column, list the residue per species and flag **diagnostic substitutions** (a residue shared by all members of a group such as Primates, Rodents, Birds, Fish and absent from all others). These feed Phase 7's evolutionary discussion. Columns: `aln_col, human_uniprot_pos, classical_pos, group, diagnostic_residue, other_residues_seen`.

### Task 4.7: Tests `tests/test_phase4.py`

```python
import pandas as pd, numpy as np
from cytc.utils import REPO_ROOT
R = REPO_ROOT / "results/conservation"

def test_per_position_shapes_and_ranges():
    df = pd.read_csv(R / "per_position_vert.csv")
    assert (df.identity_frac.between(0, 1)).all()
    assert (df.shannon_entropy_bits >= -1e-9).all()
    assert (df.norm_entropy.between(0, 1 + 1e-9)).all()

def test_invariant_means_zero_entropy():
    df = pd.read_csv(R / "per_position_vert.csv")
    inv = df[df["class"] == "invariant"]
    assert len(inv) > 0 and np.allclose(inv.shannon_entropy_bits, 0)

def test_heme_positions_conserved():
    pm = pd.read_csv(REPO_ROOT / "results/tables/position_map.csv")
    df = pd.read_csv(R / "per_position_vert.csv").merge(pm, on="aln_col")
    for classical in (14, 17, 18):          # Cys14, Cys17, His18
        row = df[df.classical_pos == classical].iloc[0]
        assert row["class"] == "invariant", classical
```

(If one of these three is not invariant in the *vertebrate* data, that is a **data problem** worth investigating, because the thioether-bond cysteines and axial His should be absolutely conserved.)

## 4. Expected outputs

| Path | Description |
|------|-------------|
| `results/tables/position_map.csv` | alignment <-> human UniProt <-> classical <-> yeast mapping |
| `results/conservation/per_position_{vert,all}.csv` | all per-column metrics |
| `results/conservation/regions.csv` | conserved and variable regions |
| `results/tables/conservation_summary.csv`, `results/conservation/summary.md` | headline numbers for RQ1/RQ2 |
| `results/tables/lineage_specific_substitutions.csv` | diagnostic substitutions |
| `results/figures/fig_conservation_profile.*`, `fig_entropy_bar.*`, `fig_sequence_logo_{vert,all}.*`, `fig_conservation_class_*.*`, `fig_conservation_heatmap_per_species.*` | figures (deliverable: "visualization of sequence conservation") |

## 5. Validation checks

- [ ] Heme Cys14/Cys17/His18 (classical) map to UniProt 15/18/19 in the human sequence, and are invariant in vertebrates.
- [ ] Entropy bounds correct (0 for invariant columns; <= log2(20)).
- [ ] Counts of classes sum to the number of non-gap-rich columns.
- [ ] Both vertebrate-only and all-17 results exist and are clearly labelled in every file name and figure title.
- [ ] Logos use the correct position labels; no off-by-one (spot-check 3 positions manually against the FASTA).
- [ ] Summary numbers (e.g. "x of y positions invariant") are recorded in `summary.md` for direct reuse in the report.
- [ ] `pytest -q tests/test_phase4.py` passes.

## 6. Common pitfalls

- **Off-by-one numbering** (UniProt vs. classical). Always test against the heme motif anchor.
- **Counting gaps as residues** in entropy. Compute over non-gap sequences and report `gap_fraction` separately.
- **Letting the yeast N-terminal extension distort the vertebrate profile.** This is why the vertebrate-only alignment is the primary basis.
- **Overclaiming from small n.** With 15 sequences, "invariant" means invariant in *this sample*. Say so in the interpretation.
- **Smoothing before classification.** Don't. Smoothing is for the plots only.
- **Sequence logo unreadable.** Split into rows; keep figure width about 12 in.

## 7. Done when

- [ ] Section 5 fully ticked.
- [ ] `make phase4` reproduces everything.
- [ ] `PROGRESS.md` updated; Phase 4 committed.
