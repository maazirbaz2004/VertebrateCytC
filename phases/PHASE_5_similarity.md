# PHASE 5: Pairwise Sequence Similarity Analysis

**Read `CLAUDE.md` first. Prerequisite: Phase 4 is done.**

## 0. Objective

Compute pairwise **identity**, **similarity**, and **evolutionary distance** matrices for all species; visualize them as heatmaps and a clustered heatmap; compare closely related and distantly related pairs quantitatively; and test whether sequence identity decreases with divergence time.

**Serves:** Objective 4; **RQ3** (how similarity varies between close and distant species) and, indirectly, RQ4.

## 1. Prerequisites

- `results/alignments/primary_all.fasta` passes Phase 3 tests.
- `data/processed/id_map.csv` and `config/species.csv` exist (with `group` and `approx_mya_from_human`).

## 2. Inputs

- `results/alignments/primary_all.fasta` (17 sequences).
- `config.yaml` -> `similarity` (`gap_open`, `gap_extend`, `substitution_matrix`).
- Biopython `Bio.Align.substitution_matrices.load("BLOSUM62")`.

## 3. Tasks

### Task 5.1: Define the metrics precisely (write these definitions into `results/similarity/metric_definitions.md`)

For each pair of sequences (i, j) in the alignment:

| Metric | Definition |
|--------|-----------|
| `n_compared` | alignment columns where **both** sequences have a residue (no gap in either) |
| `identity_pct` | `100 * identical / n_compared` |
| `identity_pct_overall` | `100 * identical / min(len_i, len_j)` (ungapped lengths; penalizes terminal extensions) |
| `similarity_pct` | `100 * (columns with BLOSUM62 score > 0) / n_compared` (a "positive" substitution, which includes identities) |
| `p_distance` | `1 - identity_pct/100` |
| `poisson_distance` | `-ln(1 - p)` (corrects for multiple hits; `p` = p_distance; handle `p=0 -> 0`; guard `p>=1`) |
| `n_substitutions` | number of non-identical compared columns |
| `blosum62_score` | sum of BLOSUM62 scores over compared columns |

State clearly which is the **headline** metric. Default: `identity_pct` over `n_compared` columns (gap-free pairwise), plus `similarity_pct`.

### Task 5.2: Compute matrices (`scripts/50_similarity.py`)

Write square, symmetric CSV matrices (rows and columns in the fixed species order) to `results/similarity/`:

- `identity_matrix.csv` (percent)
- `similarity_matrix.csv` (percent)
- `n_substitutions_matrix.csv` (integer)
- `p_distance_matrix.csv`
- `poisson_distance_matrix.csv`
- `blosum62_score_matrix.csv`
- `pairwise_long.csv` (**tidy** form: `sp1, sp2, group1, group2, relationship, identity_pct, similarity_pct, n_substitutions, poisson_distance, n_compared, mya_divergence`)

Diagonal: 100 for identity and similarity, 0 for distances.

**Cross-validation (required):** independently recompute identity for **all pairs** using a *pairwise global alignment* on the unaligned sequences (`Bio.Align.PairwiseAligner`, mode `global`, `substitution_matrix = BLOSUM62`, `open_gap_score = gap_open`, `extend_gap_score = gap_extend`). Compare to the MSA-derived values: report the correlation and the maximum absolute difference in `results/similarity/crosscheck.md`. Expect near-perfect agreement for the vertebrate pairs (differences arise only where gap placement differs; investigate any pair with > 2 percentage points difference).

Also, if Phase 3 produced `all_clustalo.distmat`, compare it to your identity matrix as a third check.

### Task 5.3: Annotate pair relationships

In `pairwise_long.csv`, label every pair's `relationship` using the species table's groups:

| Label | Rule |
|-------|------|
| `within_primates` | both Primate |
| `within_rodents` | both Rodent |
| `within_other_mammals` | both mammals but not both primate / both rodent (not within a single group) |
| `primate_vs_rodent` | one primate, one rodent |
| `mammal_vs_mammal` | any two mammals (a superset label, to be used for summaries; keep as a separate boolean column `both_mammals`) |
| `within_birds` | both Bird |
| `within_reptiles` | both Reptile |
| `within_fish` | both Fish |
| `mammal_vs_bird` | one mammal, one bird |
| `mammal_vs_reptile` | ... |
| `mammal_vs_amphibian`, `mammal_vs_fish`, `bird_vs_fish`, ... | analogous labels for every cross-class combination |
| `vertebrate_vs_drosophila` | one vertebrate vs *Drosophila* |
| `vertebrate_vs_yeast` | one vertebrate vs yeast |
| `drosophila_vs_yeast` | the two outgroups |

Implementation hint: define a `clade` column (Mammal, Bird, Reptile, Amphibian, Fish, Insect, Fungus) derived from `group`, and produce labels from sorted clade pairs.

### Task 5.4: Group-level statistics (answers RQ3)

Create `results/tables/similarity_group_summary.csv`: for each `relationship` label: number of pairs, mean, median, min, max, and SD of `identity_pct`, plus the species pair with min and max identity.

Also a **human-centric table** `results/tables/identity_to_human.csv`: each species, group, `identity_pct` to human, `similarity_pct` to human, `n_substitutions` to human, `mya_divergence`, sorted by identity (descending). This is the single most quotable table for the report.

Specific comparisons to report explicitly in `results/similarity/summary.md`:
- most similar pair overall (and whether identical, i.e. 100%);
- all **identical pairs** (e.g. 100% identity, listing species);
- least similar vertebrate pair;
- mean identity within mammals vs. mammal-vs-bird vs. mammal-vs-fish vs. vertebrate-vs-yeast;
- **closely related** (human-macaque, mouse-rat, chicken-zebra finch, zebrafish-fugu) vs. **distantly related** (human-fish, human-yeast) pairs.

### Task 5.5: Test the identity-vs-divergence-time relationship

Using `mya_divergence` (from `species.csv`, pairwise divergence times; for pairs not involving human, derive from the phylogeny: the divergence time of the pair equals the age of their most recent common ancestor. Implement with the reference tree of Phase 7 **or** a hard-coded table of node ages labelled "approximate; replace with TimeTree values when available"; record the source):

1. **Scatter:** identity (and Poisson distance) against divergence time, for all pairs and for human-vs-X only.
2. **Spearman rank correlation** (expect strongly negative for identity, positive for distance), with p-value (use `scipy.stats.spearmanr`). Also Pearson on `poisson_distance` vs time as a **molecular-clock** test; report R^2 and slope (substitutions per site per million years).
3. Identify **outliers** (species whose distance to human deviates notably from the clock expectation; for example fast-evolving lineages) and list the top 3 residuals. Discuss rate heterogeneity in Phase 7.
4. Caveat in the output: pairs are **not independent** (shared ancestry), so the p-values are descriptive only.

### Task 5.5b: Figures

1. **`fig_identity_heatmap.png/.svg`**: annotated heatmap (values shown in cells if <= 17 species, font 6-7 pt), rows/cols in the fixed biological order, `cividis` or `viridis`, with a coloured side-bar for clade membership.
2. **`fig_similarity_heatmap.png/.svg`**: the same for BLOSUM62 similarity.
3. **`fig_identity_clustermap.png/.svg`**: `seaborn.clustermap` (average linkage on `poisson_distance` or `100 - identity`) so that clusters emerge *without* using the known taxonomy. Add the clade side-colours. This is a quick preview of Phase 6.
4. **`fig_identity_to_human.png/.svg`**: ordered horizontal bar chart of identity to human, coloured by clade, annotated with divergence time.
5. **`fig_identity_vs_time.png/.svg`**: scatter with regression line, labelled outliers.
6. **`fig_identity_distribution_by_relationship.png/.svg`**: box/strip plot of `identity_pct` grouped by `relationship` (ordered from close to distant).

### Task 5.6: Tests `tests/test_phase5.py`

```python
import numpy as np, pandas as pd
from cytc.utils import REPO_ROOT
S = REPO_ROOT / "results/similarity"

def test_matrix_symmetric_and_diag():
    m = pd.read_csv(S / "identity_matrix.csv", index_col=0)
    assert (m.index == m.columns).all()
    assert np.allclose(m.values, m.values.T)
    assert np.allclose(np.diag(m.values), 100)
    assert ((m.values >= 0) & (m.values <= 100 + 1e-9)).all()

def test_similarity_geq_identity():
    i = pd.read_csv(S / "identity_matrix.csv", index_col=0)
    s = pd.read_csv(S / "similarity_matrix.csv", index_col=0)
    assert (s.values + 1e-9 >= i.values).all()

def test_long_table_complete():
    n = len(pd.read_csv(S / "identity_matrix.csv", index_col=0))
    long = pd.read_csv(S / "pairwise_long.csv")
    assert len(long) == n * (n - 1) // 2
```

## 4. Expected outputs

| Path | Description |
|------|-------------|
| `results/similarity/*_matrix.csv` (six matrices), `pairwise_long.csv` | deliverable: pairwise similarity results |
| `results/similarity/metric_definitions.md`, `crosscheck.md`, `summary.md` | documentation and headline findings |
| `results/tables/similarity_group_summary.csv`, `identity_to_human.csv` | report tables |
| `results/figures/fig_identity_heatmap.*`, `fig_similarity_heatmap.*`, `fig_identity_clustermap.*`, `fig_identity_to_human.*`, `fig_identity_vs_time.*`, `fig_identity_distribution_by_relationship.*` | figures |

## 5. Validation checks

- [ ] All matrices symmetric, diagonal correct, values within bounds, same species order as the FASTA.
- [ ] Similarity >= identity for every pair.
- [ ] MSA-derived identity matches the pairwise-alignment cross-check (documented in `crosscheck.md`).
- [ ] Every pair has a `relationship` label (no `NaN`); pair count = n(n-1)/2 = 136 for 17 sequences.
- [ ] Identical pairs are listed explicitly (not hidden).
- [ ] The identity-vs-time analysis states the divergence-time source and its caveats.
- [ ] Mean identity trend is monotone in the expected direction (close > intermediate > distant), or the exception is reported and explained.
- [ ] `pytest -q tests/test_phase5.py` passes.

## 6. Common pitfalls

- **Counting gap-gap columns** as matches. They are excluded by `n_compared`.
- **Divergence by percent identity alone is not additive.** Use Poisson distance for modelling and the clock test.
- **Heatmap order.** Never let seaborn reorder rows in the plain heatmaps. Only `clustermap` clusters, and that figure is labelled as such.
- **Misreading 100% identity** as a pipeline bug. Closely related mammals genuinely can have identical cytochrome c. Verify against the sequence and report it.
- **Treating pairs as independent samples** in significance tests. Note the caveat.

## 7. Done when

- [ ] Section 5 fully ticked.
- [ ] `make phase5` reproduces everything.
- [ ] `PROGRESS.md` updated; Phase 5 committed.
