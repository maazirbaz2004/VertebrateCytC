# PHASE 7: Biological Interpretation and Validation Against the Known Species Tree

**Read `CLAUDE.md` first. Prerequisite: Phase 6 is done.**

## 0. Objective

Turn the numbers into biology:
1. Link conserved and variable positions to known **structural and functional features** of cytochrome c (using UniProt annotations pulled programmatically, plus the literature-based checklist below).
2. Compare the inferred trees with the **accepted vertebrate species tree** quantitatively (Robinson-Foulds, clade-by-clade monophyly checks).
3. Evaluate the rate of evolution (molecular clock) and lineage-specific changes.
4. **Answer the four research questions** with specific numbers, and write an honest limitations section.

**Serves:** Objective 6 and all research questions (RQ1-RQ4).

## 1. Prerequisites

- Phases 4, 5 and 6 outputs exist and pass tests.
- Internet access (UniProt feature lookup; optional TimeTree).

## 2. Inputs

| File | Use |
|------|-----|
| `results/conservation/per_position_{vert,all}.csv`, `regions.csv`, `results/tables/position_map.csv`, `lineage_specific_substitutions.csv` | conservation |
| `results/similarity/pairwise_long.csv`, `results/tables/identity_to_human.csv` | similarity |
| `results/trees/iqtree/cytc_all.treefile`, `cytc_all.collapsed.nwk`, NJ and FastTree trees, `sensitivity_summary.csv` | trees |
| `config/species.csv` | groups, divergence times |

## 3. Tasks

### Task 7.1: Fetch functional annotation for human cytochrome c (`scripts/70_interpretation.py`, part A)

Download the UniProt JSON for the human entry (accession from `dataset_provenance.csv`, expected P99999) and extract the `features` list. Keep these feature types: `Binding site`, `Modified residue`, `Natural variant`, `Site`, `Helix`, `Beta strand`, `Turn`, `Domain`, `Region`, `Disulfide bond` and `Cross-link` (heme attachment may appear as `Binding site` or `Modified residue`/`Cross-link`; inspect and keep everything relevant). Write `results/interpretation/uniprot_features_human.csv` with `type, description, start, end, ligand (if any)`.

Also record the entry's `cc_function`, `cc_subcellular_location` and `cc_disease` comments as text in `results/interpretation/uniprot_comments.md` (quote nothing long; paraphrase for the report).

Then **join features with conservation**: for every annotated position, attach `class`, `identity_frac`, `shannon_entropy_bits` from `per_position_vert.csv` using `position_map.csv`. Output `results/tables/functional_sites_conservation.csv`.

### Task 7.2: Literature-based functional checklist (verify, do not assume)

For each item, look up the residue in the data and record **what the alignment actually shows**. Positions are classical numbering (= UniProt - 1). The list is a **guide to what to check**, not a set of facts to paste into the report. Any claim in the report must be supported by your own alignment evidence and/or the UniProt annotation retrieved in Task 7.1.

| Feature | Residues to check (classical) | Role (as commonly described) |
|---------|------------------------------|------------------------------|
| Heme covalent attachment (thioether bonds) | Cys14, Cys17 | bind heme vinyl groups |
| Proximal axial ligand | His18 | binds heme iron |
| Distal axial ligand | Met80 | binds heme iron, key for redox potential and the peroxidase/apoptotic switch |
| Heme pocket H-bond network | Asn52, Trp59, Tyr67, Thr78 | stabilize heme and modulate reactivity |
| Close packing around heme | Gly6, Gly41, Pro30, Phe82 | glycines allow close approach; hydrophobic contacts |
| Lysine patch (partner binding) | Lys13, Lys27, Lys72, Lys73, Lys79, Lys86, Lys87 (and nearby) | electrostatic docking with complex III, cytochrome c oxidase, Apaf-1 |
| Regulatory phosphosites (human) | Tyr48, Tyr97 (UniProt numbering may differ by 1; check UniProt "Modified residue") | regulation of electron transfer / apoptosis |
| Disease variants (human, thrombocytopenia 4) | as listed under UniProt "Natural variant" | check whether the sites are conserved |

For **each item**: report vertebrate conservation class, the residue seen in each species (compact string), whether it is also conserved in Drosophila/yeast, and any exceptions. Produce `results/tables/functional_checklist_results.csv` with columns `feature, classical_pos, uniprot_pos, aln_col, human_residue, vert_class, vert_residue_counts, all17_residue_counts, supports_expected_role (yes/partly/no), comment`.

Examples of the kind of finding to make explicit:
- "All 15 vertebrates and both outgroups share Cys/Cys/His at the heme-binding positions: absolutely conserved across the deep tree."
- "Lysine patch positions X, Y, Z are conserved K in all vertebrates but vary in yeast."

If a checklist item is **not** conserved in your data, say so plainly and consider alignment error, paralog contamination, or true divergence as explanations (check in that order).

### Task 7.3: Variable-region analysis

For the most variable vertebrate positions (top ~10 from Phase 4):
- Are they at the **protein surface** vs. core? (No 3D structure is required. Use a hydropathy/charge argument, and optionally download PDB `3ZCF` (human cytochrome c) or `1HRC` (horse) and compute solvent accessibility with Biopython `Bio.PDB.SASA` / `ShrakeRupley`. **Optional.** If PDB retrieval fails, skip with a log note.)
- Are substitutions **conservative** (BLOSUM62 >= 0) or radical?
- Do they correlate with **lineages** (use `lineage_specific_substitutions.csv`; for example a residue unique to birds, primates or fish)? List a short table of diagnostic substitutions that support clades.

### Task 7.4: Reference species tree and tree comparison

**7.4.1 Build the reference (accepted) species tree** as a Newick file `config/reference_species_tree.nwk`. Use the following topology (consensus of molecular phylogenies of vertebrates; the placement of Laurasiatheria members dog/horse/cattle is contested, so it is **deliberately left as a polytomy**; the main groups reflect standard taxonomy). Branch lengths are not needed.

```
(Saccharomyces_cerevisiae,(Drosophila_melanogaster,((Danio_rerio,Takifugu_rubripes),(Xenopus_tropicalis,((Monodelphis_domestica,(((Homo_sapiens,Macaca_mulatta),(Mus_musculus,Rattus_norvegicus)),(Canis_lupus_familiaris,Equus_caballus,Bos_taurus))),(Anolis_carolinensis,(Chrysemys_picta,(Gallus_gallus,Taeniopygia_guttata))))))));
```

Before use: **parse it with Biopython and assert it has exactly the same leaf set as the dataset** (if you substituted species in Phase 1, edit the Newick accordingly and log it). Document the source for the topology (e.g. NCBI Taxonomy common tree and/or TimeTree 5) in `results/interpretation/reference_tree_notes.md`. Note that the turtle is placed as sister to birds (Archelosauria) following molecular studies, and that this differs from older morphological groupings. Only keep these claims if you can confirm them from a source consulted during the run.

**7.4.2 Quantitative comparison.** For each inferred tree (primary ML, collapsed ML, NJ, FastTree) vs. the reference tree:
- **Robinson-Foulds distance** (unrooted and rooted), **normalized RF** (RF / max possible), with a shared `TaxonNamespace`. Because the reference has a polytomy, also report RF between the reference and the **collapsed** inferred tree and explain that polytomies lower the maximum.
- **Number of reference clades recovered** / number in the reference.
- **Clade-by-clade monophyly checks** (table `results/tables/clade_recovery.csv`; columns `clade, members, recovered_in_ML (Y/N), UFBoot, recovered_in_NJ, NJ_bootstrap, recovered_in_FastTree, comment`) for at least:
  - Mammalia (incl. opossum), Eutheria, Primates (human, macaque), Rodentia (mouse, rat), Euarchontoglires, Laurasiatheria (dog, horse, cattle), Aves (chicken, zebra finch), Sauropsida (birds + reptiles), Amniota, Tetrapoda, Teleostei (zebrafish, fugu), Vertebrata (vs. Drosophila), Bilateria (vertebrates + Drosophila vs. yeast), Actinopterygii vs. tetrapods.
- **Where the trees disagree**, give a candidate explanation: too little signal, identical sequences, rate heterogeneity (long-branch attraction), lineage-specific accelerated evolution, paralog mix-up, or alignment issues. Support each explanation with evidence from the data (e.g. branch lengths, support values, number of informative sites).
- **Informative sites count**: number of parsimony-informative sites in the vertebrate alignment (a site with >= 2 residues each occurring in >= 2 sequences). Report this number: it bounds how much resolution is even possible. (IQ-TREE's `.iqtree` report also lists it.)

**7.4.3 Resolution analysis.** Compute for the ML tree: fraction of internal nodes with UFBoot >= 70/90, mean branch length in mammals vs. between classes, and the number of zero-length terminal branches. Summarize "where cytochrome c is informative": expected to be **deep** relationships (outgroup, fish vs. tetrapods, mammals vs. birds) rather than **shallow** ones (primate/rodent internals).

### Task 7.5: Rate of evolution

Using Phase 5 data:
- Estimate the substitution rate (Poisson distance per million years) from the human-vs-X regression. Compare to published expectations qualitatively only (cytochrome c is known as a **slowly evolving** protein whose rate is not constant across lineages, so do not assert numerical literature values unless you retrieve and cite them).
- Identify **branches with unusually many substitutions** from the ML tree (long terminal branches) and relate them to the lineage-specific substitution table. It is widely discussed that primates show a distinct burst of cytochrome c changes in their history (see Grossman et al. 2011 in the proposal's references); check whether your data shows human/macaque differing from other mammals at any positions and report only what your data shows.
- Optional **ancestral sequence reconstruction**: `iqtree2 -asr` (or the `--ancestral` flag in current versions) on the primary tree to list substitutions per branch (`.state` file). Summarize the number of changes on the branch leading to primates/rodents/birds, if feasible.

### Task 7.6: Answer the research questions (`results/interpretation/research_questions_answers.md`)

Write a structured answer for each, **each with specific numbers, figure references, and caveats**:

**RQ1: How conserved is cytochrome c across the selected vertebrate species?**
Provide: % identical columns in the vertebrate alignment; mean/min pairwise identity; mean entropy; whether identical sequences occur; contrast with the all-17 alignment (deep conservation).

**RQ2: Which positions or regions are highly conserved?**
Provide: list of invariant positions (both numberings); conserved regions with coordinates; overlap with heme-binding residues, Met80 loop, lysine patch, and other known functional sites (from Tasks 7.1-7.2); the most variable positions and what they have in common.

**RQ3: How does similarity vary between closely and distantly related species?**
Provide: the identity-to-human table in brief; mean identity by relationship class; correlation with divergence time (Spearman rho, p, slope); explicit examples (human-macaque, mouse-rat vs. human-zebrafish vs. human-yeast).

**RQ4: Can cytochrome c sequence information be used to infer meaningful phylogenetic relationships?**
Provide: RF distances; number of reference clades recovered; which relationships are strongly supported vs. not; conclusion framed as **"yes for deep relationships, limited for shallow ones"** *only if the data supports it*; otherwise state what the data does show.

Also write the **hypothesis check** sentence: "Closely related species show greater sequence similarity and cluster together" -> supported / partially supported / not supported, with the evidence.

### Task 7.7: Limitations (`results/interpretation/limitations.md`)

Must address, with specifics from the run:
1. Single short protein -> limited informative sites; low bootstrap for shallow nodes.
2. One sequence per species; sampling is not uniform across vertebrates.
3. Identical sequences -> unresolved relationships; arbitrary order within them.
4. Rate heterogeneity across lineages; possible compositional bias.
5. Possible paralog and annotation issues (database-derived sequences; some may be predicted).
6. Outgroup distance (yeast) is very large; rooting could be affected by long-branch effects (cross-checked with the Drosophila outgroup and the vertebrate-only tree).
7. Divergence times are approximate and from secondary sources.
8. A gene tree is not a species tree.
9. No structural/experimental validation (computational only, per the proposal's scope).

### Task 7.8: Figures

1. `fig_functional_sites_on_conservation.png/.svg`: Phase 4 profile with functional sites annotated (heme ligands, Met80, lysine patch, phosphosites, disease variants).
2. `fig_tree_vs_reference.png/.svg`: tanglegram-style or side-by-side plot (reference vs. ML) with recovered clades highlighted and conflicting clades flagged (a simple side-by-side with coloured clade labels is enough; do not rely on exotic packages).
3. `fig_clade_support_summary.png/.svg`: bar/dot plot of support for each reference clade across methods.
4. `fig_rf_distance_summary.png/.svg`: RF distances per method.
5. Optional: `fig_structure_conservation.png` (PDB colored by conservation) only if PyMOL/py3Dmol is available, otherwise skip.

### Task 7.9: Tests `tests/test_phase7.py`

```python
import pandas as pd
from Bio import Phylo
from cytc.utils import REPO_ROOT

def test_reference_tree_leaves_match_dataset():
    ref = {t.name for t in Phylo.read(REPO_ROOT / "config/reference_species_tree.nwk", "newick").get_terminals()}
    data = set(pd.read_csv(REPO_ROOT / "data/processed/id_map.csv").tree_id)
    assert ref == data

def test_required_outputs_exist():
    for rel in ["results/tables/clade_recovery.csv",
                "results/tables/functional_checklist_results.csv",
                "results/interpretation/research_questions_answers.md",
                "results/interpretation/limitations.md"]:
        assert (REPO_ROOT / rel).exists(), rel

def test_rq_file_answers_all_four():
    txt = (REPO_ROOT / "results/interpretation/research_questions_answers.md").read_text()
    for q in ("RQ1", "RQ2", "RQ3", "RQ4"):
        assert q in txt
```

## 4. Expected outputs

| Path | Description |
|------|-------------|
| `config/reference_species_tree.nwk`, `results/interpretation/reference_tree_notes.md` | reference topology + source |
| `results/interpretation/uniprot_features_human.csv`, `uniprot_comments.md` | UniProt functional annotation |
| `results/tables/functional_sites_conservation.csv`, `functional_checklist_results.csv` | function <-> conservation mapping |
| `results/tables/clade_recovery.csv`, RF summary tables | tree validation |
| `results/interpretation/research_questions_answers.md`, `limitations.md` | **the scientific conclusions** |
| `results/figures/fig_functional_sites_on_conservation.*`, `fig_tree_vs_reference.*`, `fig_clade_support_summary.*`, `fig_rf_distance_summary.*` | figures |

## 5. Validation checks

- [ ] Reference tree parses and its leaves exactly match the dataset.
- [ ] Every statement about a functional residue in the answers file is backed by a row in `functional_checklist_results.csv` or UniProt features, **not** by the checklist text alone.
- [ ] RF distances (raw and normalized) reported for every tree vs. the reference.
- [ ] Clade recovery table covers all clades listed in Task 7.4.2, with both ML and NJ.
- [ ] Disagreements between trees and reference are listed with an evidenced explanation.
- [ ] All four RQs are answered with numbers; the hypothesis is explicitly judged supported / partially / not.
- [ ] Limitations cover all nine points.
- [ ] Every number quoted in `research_questions_answers.md` can be traced to a CSV in `results/`.
- [ ] `pytest -q tests/test_phase7.py` passes.

## 6. Common pitfalls

- **Writing biology from memory.** Residue roles in the checklist are guideposts. If the data or UniProt disagrees, the data wins.
- **Reference tree errors.** A wrong reference makes RF meaningless. Validate leaf names, and document contested nodes (Laurasiatheria polytomy; turtle placement).
- **Comparing rooted and unrooted trees inconsistently.** State which one each RF number uses.
- **Polytomies inflate RF.** Compare against the collapsed tree too.
- **Over-claiming.** Low support means "not resolved", not "the tree is wrong" or "right".
- **Cherry-picking** the method that matches the reference best. Report **all** methods.
- **Numbering confusion** in the text. State numbering every time a position is given.

## 7. Done when

- [ ] Section 5 fully ticked.
- [ ] `make phase7` reproduces tables and figures (the Markdown answer files may be written by the agent from the generated numbers).
- [ ] `PROGRESS.md` updated; Phase 7 committed.
