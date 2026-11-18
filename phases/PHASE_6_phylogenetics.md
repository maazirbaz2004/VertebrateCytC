# PHASE 6: Phylogenetic Tree Construction

**Read `CLAUDE.md` first. Prerequisite: Phase 5 is done.**

## 0. Objective

Reconstruct phylogenetic trees from the cytochrome c alignment using **three complementary methods** (distance-based Neighbor-Joining, fast approximate ML with FastTree, and full maximum likelihood with IQ-TREE including model selection and ultrafast bootstrap), root them with the outgroups, assess support, run sensitivity analyses, and export publication and iTOL-ready files.

**Serves:** Objective 5; **RQ4** (can cytochrome c infer meaningful phylogeny?). The comparison to the true species tree happens in Phase 7.

## 1. Prerequisites

- `results/alignments/primary_all.fasta` (and trimmed version) pass Phase 3 tests.
- Phase 5 matrices exist (the NJ tree reuses the same distance concepts).
- Tools `iqtree*` and `FastTree` work (Phase 0).

## 2. Inputs

- `results/alignments/primary_all.fasta`, `primary_all_trimmed.fasta`, `primary_vert.fasta`.
- `config.yaml` -> `phylogenetics`, `tools`, `project.random_seed`, `project.n_threads`.
- `data/processed/id_map.csv` (for clade colours).

## 3. Background (so the agent understands the expected difficulty)

Cytochrome c is only about 104 residues and extremely conserved, so:
- many close relatives (primates, rodents, perhaps ruminants/horse/dog) will share identical or near-identical sequences, leading to **polytomies or zero-length branches**;
- bootstrap support within mammals will be **low** (there is simply not enough informative signal), while deeper splits (fish / tetrapods, vertebrate / outgroup) should be better supported;
- the tree will likely **not** reproduce every known relationship. This is an expected, reportable result and not a failure. Never tweak settings just to force the textbook topology.

## 4. Tasks

### Task 6.1: Neighbor-Joining tree with bootstrap (`scripts/60_phylogenetics.py`, part A)

Using Biopython:

```python
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.Consensus import bootstrap_trees, get_support, majority_consensus

aln = AlignIO.read("results/alignments/primary_all.fasta", "fasta")
calc = DistanceCalculator("blosum62")            # protein distance
constructor = DistanceTreeConstructor(calc, "nj")
nj_tree = constructor.build_tree(aln)
```

Requirements:
1. Build NJ trees with **two distance definitions** and compare: (a) Biopython `blosum62`; (b) a **custom Poisson-corrected distance** from `results/similarity/poisson_distance_matrix.csv`, converted to a `DistanceMatrix` (`Bio.Phylo.TreeConstruction.DistanceMatrix` takes the lower-triangular list-of-lists). Pick (b) as the primary NJ tree (it is the standard for protein NJ) and report whether (a) differs.
2. **Bootstrap**: `bootstrap_trees(aln, times=cfg["bootstrap_replicates_nj"], tree_constructor=constructor)` (set `random.seed(cfg seed)` first); compute support with `get_support(target_tree, trees)`; support values are 0-1 in Biopython. Convert to percent and store as clade confidence.
3. **Negative branch lengths** (a known NJ artefact): set to 0 and **log how many were clamped**.
4. **Root** with `tree.root_with_outgroup(primary_root_outgroup)`, and as a second view root on the *outgroup pair* (`[Saccharomyces_cerevisiae, Drosophila_melanogaster]`) only if they form a clade; otherwise log that they do not and keep the single-outgroup rooting.
5. Ladderize (`tree.ladderize()`) for display only. Do not alter topology.
6. Save: `results/trees/nj/nj_poisson.nwk` (Newick with support as internal-node labels), `nj_blosum62.nwk`, and `nj_bootstrap_consensus.nwk`.

### Task 6.2: FastTree (quick ML cross-check, part B)

```bash
FastTree -lg -gamma < results/alignments/primary_all.fasta > results/trees/fasttree/fasttree_lg_gamma.nwk
```
(Binary name from config. FastTree reports **Shimodaira-Hasegawa-like local supports** (0-1), different from bootstrap; label them accordingly.) Log the command, version, and log-likelihood from stderr.

### Task 6.3: Maximum likelihood with IQ-TREE (primary tree, part C)

```bash
iqtree2 -s results/alignments/primary_all.fasta \
        -m MFP \
        -B 1000 -alrt 1000 -bnni \
        -T AUTO --threads-max 4 \
        -seed 42 \
        -keep-ident \
        -o Saccharomyces_cerevisiae \
        --prefix results/trees/iqtree/cytc_all
```

Meaning: `-m MFP` = ModelFinder Plus (selects the best-fit substitution model by BIC), `-B 1000` = ultrafast bootstrap, `-alrt 1000` = SH-aLRT, `-bnni` = reduce UFBoot overestimation, `-keep-ident` = **keep identical sequences** in the analysis (otherwise IQ-TREE drops duplicates and the leaves vanish from the tree), `-o` = outgroup rooting.

Read parameters from config rather than hard-coding. Then:
1. Parse `cytc_all.iqtree` for: best-fit model, log-likelihood, number of parameters, tree length, and the list of **identical-sequence warnings**; save to `results/tables/iqtree_summary.csv`.
2. The `.treefile` has node labels `SH-aLRT/UFBoot`, e.g. `87.3/94`. Parse both into separate numbers per internal node.
3. Convert the tree into a table `results/tables/ml_clade_support.csv`: for each internal node: the descendant leaf set, SH-aLRT, UFBoot, and a `supported` flag (UFBoot >= `support_threshold` **and** SH-aLRT >= 80).
4. Also run IQ-TREE on **`primary_vert.fasta`** (vertebrates only; unrooted; root later on the fish lineage for display only) to check that the outgroups do not distort vertebrate relationships: prefix `results/trees/iqtree/cytc_vert`, without `-o`.
5. **IQ-TREE model sensitivity:** repeat the primary analysis with a fixed model, e.g. `-m LG+G4`, and compare to the MFP result.

### Task 6.4: Sensitivity analyses (record each result in `results/trees/sensitivity_summary.csv`)

For each variant, compute the **Robinson-Foulds distance** (use `dendropy.calculate.treecompare.symmetric_difference` with a **shared `TaxonNamespace`**, unrooted comparison via `is_bipartitions_updated` or by passing `is_rooted=False`) to the primary IQ-TREE tree and note which clades change:

| Variant | How |
|---------|-----|
| NJ (Poisson) | Task 6.1 |
| NJ (BLOSUM62 distance) | Task 6.1 |
| FastTree LG+G | Task 6.2 |
| IQ-TREE, vertebrates only | Task 6.3.4 |
| IQ-TREE, fixed LG+G4 | Task 6.3.5 |
| IQ-TREE, **MUSCLE** alignment | re-run using `all_muscle.fasta` |
| IQ-TREE, **trimmed** alignment | `primary_all_trimmed.fasta` |
| IQ-TREE, different random seeds (42, 1, 2) | stability of topology |

Report the number of variants that recover the same unrooted topology as the primary, and list the **unstable clades** (appear in some trees but not others). That stable-vs-unstable summary is a key result for RQ4.

### Task 6.5: Handling identical sequences and polytomies

- Using `results/tables/identical_sequence_groups.csv` (Phase 2), check that members of an identical group sit together with **near-zero branch lengths** in the ML tree.
- Collapse branches with length < 1e-6 (or with UFBoot < 50) into **polytomies in a separate "collapsed" version** of the tree for the honest comparison in Phase 7. Save as `results/trees/iqtree/cytc_all.collapsed.nwk`. Keep the fully-resolved version too.
- Document that resolution *among identical sequences is arbitrary*, so do not interpret it.

### Task 6.6: Tree visualisation (`results/figures/`)

Using Biopython + matplotlib (or `toytree` / `ete4` if installed and working, but **do not** make the pipeline depend on a fragile GUI library):
1. `fig_tree_ml_rectangular.png/.svg`: primary ML tree, rooted on the outgroup, tips coloured by clade (Primate, Rodent, other mammals, Bird, Reptile, Amphibian, Fish, Insect, Fungus), UFBoot support shown at nodes (hide values < 50, with a note), scale bar in substitutions/site, ladderized.
2. `fig_tree_nj_rectangular.png/.svg`: NJ tree, same styling (bootstrap %).
3. `fig_tree_comparison_NJ_vs_ML.png/.svg`: side-by-side panel of the two trees.
4. `fig_tree_unrooted_ml.png/.svg`: radial/unrooted style layout if feasible (otherwise skip, with a log note).
5. `fig_tree_with_identity_heatmap.png/.svg` (optional): tree beside a conservation or identity-to-human heat strip.

### Task 6.7: iTOL-ready exports (`results/trees/itol/`)

- `cytc_ml_rooted.nwk`: Newick with support values as internal node labels (iTOL reads them as bootstrap).
- `itol_colorstrip_clade.txt`: a **DATASET_COLORSTRIP** annotation colouring leaves by clade. Template:

```
DATASET_COLORSTRIP
SEPARATOR COMMA
DATASET_LABEL,Clade
COLOR,#ff0000
LEGEND_TITLE,Clade
LEGEND_SHAPES,1,1,1,1,1,1,1,1,1
LEGEND_COLORS,#1b9e77,#d95f02,#7570b3,#e7298a,#66a61e,#e6ab02,#a6761d,#666666,#1f78b4
LEGEND_LABELS,Primate,Rodent,OtherMammal,Bird,Reptile,Amphibian,Fish,Insect,Fungus
DATA
Homo_sapiens,#1b9e77,Primate
...
```
- `itol_labels_common_names.txt`: a **LABELS** annotation that replaces tip IDs with "Genus species (Common name)".
- `README_itol.md`: 4-line instructions (upload the Newick, drag the annotation files, set rooting, export SVG).

### Task 6.8: Tests `tests/test_phase6.py`

```python
from Bio import Phylo
import pandas as pd
from cytc.utils import REPO_ROOT

T = REPO_ROOT / "results/trees"
EXPECTED = set(pd.read_csv(REPO_ROOT / "data/processed/id_map.csv").tree_id)

def leaves(path):
    return {t.name for t in Phylo.read(path, "newick").get_terminals()}

def test_ml_tree_has_all_leaves():
    assert leaves(T / "iqtree/cytc_all.treefile") == EXPECTED

def test_nj_tree_has_all_leaves():
    assert leaves(T / "nj/nj_poisson.nwk") == EXPECTED

def test_fasttree_has_all_leaves():
    assert leaves(T / "fasttree/fasttree_lg_gamma.nwk") == EXPECTED

def test_outgroup_rooting():
    tree = Phylo.read(T / "iqtree/cytc_all.treefile", "newick")
    # root's children: one of them must be exactly the yeast outgroup
    kids = [{t.name for t in c.get_terminals()} for c in tree.root.clades]
    assert {"Saccharomyces_cerevisiae"} in kids
```

## 5. Expected outputs

| Path | Description |
|------|-------------|
| `results/trees/nj/*.nwk` | NJ trees with bootstrap |
| `results/trees/fasttree/fasttree_lg_gamma.nwk` | FastTree |
| `results/trees/iqtree/cytc_all.{treefile,iqtree,log,contree,splits.nex}` (+ `cytc_vert.*`) | primary ML outputs (deliverable: phylogenetic tree) |
| `results/trees/iqtree/cytc_all.collapsed.nwk` | polytomy-collapsed ML tree |
| `results/trees/itol/*` | iTOL files |
| `results/tables/iqtree_summary.csv`, `ml_clade_support.csv`, `results/trees/sensitivity_summary.csv` | support and robustness tables |
| `results/figures/fig_tree_*.{png,svg}` | tree figures |

## 6. Validation checks

- [ ] All trees contain exactly the 17 leaves (or the vertebrate subset), with identical-sequence species **present**.
- [ ] Root is on the yeast outgroup; Drosophila's position is reported (sister to vertebrates expected).
- [ ] ModelFinder-selected model and likelihood recorded.
- [ ] Number of NJ negative branches clamped is logged.
- [ ] RF-distance table across all sensitivity variants is produced.
- [ ] Stable and unstable clades explicitly listed.
- [ ] Support thresholds in code come from config; no ad-hoc values.
- [ ] Random seeds are logged for every stochastic step.
- [ ] `pytest -q tests/test_phase6.py` passes.

## 7. Common pitfalls

- **IQ-TREE drops identical sequences** unless you pass `-keep-ident`. Verify the leaf count.
- **Support-value semantics differ**: Biopython NJ bootstrap (0-1 -> convert to %), IQ-TREE UFBoot (%, tends to be more generous than standard bootstrap), SH-aLRT (%), FastTree local support (0-1). Never mix them in one figure legend without naming each.
- **Rooting changes interpretation, not topology.** Do not root on a vertebrate lineage in the primary tree.
- **Newick quoting.** Tip names must be `[A-Za-z0-9_]` (ensured in Phase 2).
- **Overinterpreting low support.** Collapsed/unsupported branches must be flagged in figures.
- **RF distance with different taxon namespaces** silently returns garbage. Create one `dendropy.TaxonNamespace` and pass `taxon_namespace=` when reading every tree.
- **Different IQ-TREE binary names** (`iqtree`, `iqtree2`, `iqtree3`); and option spelling (`-T`, `--threads-max`) can differ slightly between versions. Check `iqtree2 --help` if a flag is rejected and log the adjustment.

## 8. Done when

- [ ] Section 6 fully ticked.
- [ ] `make phase6` reproduces everything (allow for stochastic bootstrap differences only if seeds were not honoured, and then document).
- [ ] `PROGRESS.md` updated; Phase 6 committed.
