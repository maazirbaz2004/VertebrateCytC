# Answers to Research Questions (RQ1–RQ4) & Hypothesis Evaluation

All numbers below are derived directly from `results/` CSV tables and trees. Both **UniProt-human** (Met1 = 1) and **classical horse-heart** (UniProt − 1) numbering are stated wherever positions are cited.

---

## RQ1: How conserved is cytochrome c across the selected vertebrate species?

Cytochrome c is **exceptionally conserved** across vertebrates:
- **Vertebrates only (`primary_vert.fasta`, n = 15, 105 alignment columns):**
  - **67 / 105 (63.81%)** columns are strictly invariant (zero entropy, zero gaps).
  - **83 / 105 (79.05%)** columns are >= 90% conserved.
  - Mean per-column identity is **92.99%** (median **100.0%**), mean Shannon entropy is **0.2741 bits** (out of 4.32 bits max), and mean JSD score is **0.802**.
  - Mean vertebrate pairwise identity is **88.73%** (range **71.43%** for *Homo sapiens* vs. *Takifugu rubripes* to **100.00%** for identical rodent pair *Mus musculus* & *Rattus norvegicus*).
  - Only **19 / 105** vertebrate columns are parsimony-informative.
- **Deep conservation (`primary_all.fasta`, n = 17 including *D. melanogaster* and *S. cerevisiae*, 110 columns):**
  - **52 / 110 (47.27%)** columns remain strictly invariant across ~1.1 billion years of eukaryotic divergence (**52 / 67 = 77.6%** of vertebrate-invariant sites stay invariant in fly and yeast).
  - Mean Shannon entropy increases modestly to **0.4148 bits**; **27 / 110** columns are parsimony-informative.
- **Figures:** `fig_conservation_profile.png`, `fig_entropy_bar.png`, `fig_conservation_class_bar.png`.

---

## RQ2: Which amino-acid positions or regions are highly conserved?

- **Conserved regions:** 12 contiguous conserved blocks (>= 3 aa of invariant/highly conserved columns) were detected in vertebrates (`results/conservation/regions.csv`). The longest is a **17-aa block** (`EYLENPKKYIPGTKMIF`, UniProt **67–83**, classical **66–82**, mean conservation **0.9856**), containing the longest strictly invariant stretch of **11 aa** (`NPKKYIPGTKM`, UniProt 71–81, classical 70–80, the Met80 loop).
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
- **Clade-level gradient:** Within mammals (**93.61%**) > mammal vs. bird (**90.89%**) > mammal vs. fish (**80.78%**) > vertebrate vs. yeast (**61.41%**).
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
  - **Rate heterogeneity:** *Takifugu rubripes* (`cyc-a`, residual **+0.136**, longest vertebrate terminal branch = **0.437** subs/site) evolves markedly faster than *Danio rerio*. In primates, *Homo sapiens* and *Macaca mulatta* share **5 diagnostic primate-specific substitutions** (`I12`, `M13`, `A51`, `V84`, `E90` UniProt; classical `11, 12, 50, 83, 89`) relative to all other sampled vertebrates (`results/tables/lineage_specific_substitutions.csv`), consistent with a primate-lineage substitution burst.
- **Figures:** `fig_identity_heatmap.png`, `fig_identity_to_human.png`, `fig_identity_vs_time.png`, `fig_identity_distribution_by_relationship.png`.

---

## RQ4: Can cytochrome c sequence information be used to infer meaningful phylogenetic relationships?

**Partially — cytochrome c reliably separates major outgroups and a few distinctive clades, but lacks sufficient informative sites (19 parsimony-informative sites in vertebrates) to resolve the full vertebrate species tree.**
- **Quantitative comparison vs. accepted species tree (`results/tables/rf_vs_reference.csv`, `clade_recovery.csv`):**
  - **ML (IQ-TREE primary, `Q.YEAST+G4`):** Unrooted RF = **19** (normalized **0.68**), rooted RF = **19**, recovering **5 / 14** reference clades: `Primates, Rodentia, Laurasiatheria, Vertebrata, Bilateria`.
  - **ML (collapsed polytomies, 7 weak branches collapsed):** Unrooted RF = **12** (normalized **0.43**) — collapsing unsupported branches reduces spurious conflicts.
  - **NJ (Poisson):** Unrooted RF = **13** (normalized **0.46**), recovering **8 / 14** reference clades.
  - **FastTree (LG+G):** Unrooted RF = **21** (normalized **0.75**), recovering **4 / 14** reference clades.
- **Resolution & support:** Only **2 / 15** internal ML nodes pass joint strong support (`UFBoot >= 70` and `SH-aLRT >= 80`: **Vertebrata** 89.9/70 and **Primates** 98.9/100); **4 / 15** reach `UFBoot >= 70` and **2 / 15** reach `UFBoot >= 90`. Zero-length terminal branches (6 taxa, including identical rodents *M. musculus* & *R. norvegicus*) prevent resolving shallow relationships.
- **Why disagreements occur (`Mammalia, Eutheria, Euarchontoglires, Aves, Sauropsida, Amniota, Tetrapoda, Teleostei` not recovered in ML):**
  1. **Insufficient signal:** Only 19 parsimony-informative vertebrate sites across 15 taxa leave deep mammalian and tetrapod internodes with 0–2 supporting substitutions.
  2. **Rate heterogeneity / long-branch attraction:** Fast-evolving *Takifugu rubripes* (`cyc-a`) breaks Teleostei and Sauropsida in ML (nested inside Sauropsida with weak UFBoot 43–49) and breaks Vertebrata in NJ (pulled outside *Drosophila*).
  3. **Lineage-specific substitution bursts:** The 5 primate-specific substitutions give Primates a long stem branch (`0.087` subs/site, 98.9/100 support) that makes Primates appear artificially divergent from other mammals in single-gene distance/ML trees.
- **Figures:** `fig_tree_ml_rectangular.png`, `fig_tree_comparison_NJ_vs_ML.png`, `fig_tree_vs_reference.png`, `fig_clade_support_summary.png`, `fig_rf_distance_summary.png`.

---

## Hypothesis Evaluation

> **Hypothesis:** *"Closely related species show higher sequence similarity and cluster together on the tree."*

**Verdict: PARTIALLY SUPPORTED.**
- **Supported for pairwise similarity:** Sequence identity decreases strongly and monotonically with divergence time (human-vs-X Spearman `rho = −0.958`, mammal mean `93.6%` > mammal–bird `90.9%` > mammal–fish `80.8%` > vertebrate–yeast `61.4%`), and close relatives share 99–100% identity (human–macaque `99.05%`, mouse–rat `100.0%`).
- **Only partially supported for phylogenetic clustering:** While ML recovers Vertebrata (89.9/70), Primates (98.9/100), Rodentia (0/94), and Laurasiatheria (79.9/67), it recovers only **5 / 14** reference clades (unrooted RF = **19**, dropping to **12** when 7 unsupported branches are collapsed; NJ Poisson recovers 8/14, RF = 13) and fails to recover monophyletic Mammalia, Aves, Sauropsida, Tetrapoda, or Teleostei in ML due to ultra-high conservation (only 19 parsimony-informative vertebrate sites), identical sequences, and rate heterogeneity in *Takifugu rubripes*.
