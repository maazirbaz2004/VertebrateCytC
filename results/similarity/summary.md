# Phase 5 Pairwise Similarity & Molecular-Clock Summary (RQ3)

## 1. Key Extremes & Identical Pairs
- **Most similar pair overall:** `Mus_musculus` vs `Rattus_norvegicus` = **100.00%** identity (0 substitutions).
- **Identical pairs (100.0% identity):** 1 pair(s) — `Mus_musculus` & `Rattus_norvegicus`.
- **Least similar vertebrate pair:** `Homo_sapiens` vs `Takifugu_rubripes` = **71.43%** identity (30 substitutions).

## 2. Gradient Across Evolutionary Distance
- Mean identity **within mammals**: **93.61%**
- Mean identity **mammal vs. bird**: **90.89%**
- Mean identity **mammal vs. fish**: **80.78%**
- Mean identity **vertebrate vs. yeast**: **61.41%**

## 3. Close vs. Distant Pair Highlights
| comparison | identity_pct | similarity_pct | n_substitutions | poisson_distance | mya |
| --- | --- | --- | --- | --- | --- |
| Human vs Macaque (close primate) | 99.05 | 99.05 | 1 | 0.0096 | 29 |
| Mouse vs Rat (close rodent) | 100.0 | 100.0 | 0 | 0.0 | 20 |
| Chicken vs Zebra finch (close bird) | 94.29 | 97.14 | 6 | 0.0588 | 90 |
| Zebrafish vs Fugu (teleost fish) | 75.0 | 85.58 | 26 | 0.2877 | 230 |
| Human vs Zebrafish (distant vert) | 84.62 | 91.35 | 16 | 0.1671 | 430 |
| Human vs Fugu (distant vert) | 71.43 | 82.86 | 30 | 0.3365 | 430 |
| Human vs Yeast (outgroup) | 62.5 | 75.0 | 39 | 0.47 | 1100 |

## 4. Identity to *Homo sapiens*
| species | group | clade | identity_pct | similarity_pct | n_substitutions | poisson_distance | mya_divergence |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Homo_sapiens | Primate | Mammal | 100.0 | 100.0 | 0 | 0.0 | 0.0 |
| Macaca_mulatta | Primate | Mammal | 99.05 | 99.05 | 1 | 0.0096 | 29.0 |
| Mus_musculus | Rodent | Mammal | 91.43 | 94.29 | 9 | 0.0896 | 90.0 |
| Rattus_norvegicus | Rodent | Mammal | 91.43 | 94.29 | 9 | 0.0896 | 90.0 |
| Bos_taurus | Artiodactyl | Mammal | 90.48 | 94.29 | 10 | 0.1001 | 94.0 |
| Monodelphis_domestica | Marsupial | Mammal | 90.48 | 95.24 | 10 | 0.1001 | 160.0 |
| Canis_lupus_familiaris | Carnivore | Mammal | 89.52 | 93.33 | 11 | 0.1107 | 94.0 |
| Equus_caballus | Perissodactyl | Mammal | 88.57 | 93.33 | 12 | 0.1214 | 94.0 |
| Taeniopygia_guttata | Bird | Bird | 88.57 | 92.38 | 12 | 0.1214 | 320.0 |
| Gallus_gallus | Bird | Bird | 87.62 | 92.38 | 13 | 0.1322 | 320.0 |
| Anolis_carolinensis | Reptile | Reptile | 87.62 | 94.29 | 13 | 0.1322 | 320.0 |
| Xenopus_tropicalis | Amphibian | Amphibian | 86.67 | 93.33 | 14 | 0.1431 | 352.0 |
| Chrysemys_picta | Reptile | Reptile | 84.76 | 92.38 | 16 | 0.1653 | 320.0 |
| Danio_rerio | Fish | Fish | 84.62 | 91.35 | 16 | 0.1671 | 430.0 |
| Drosophila_melanogaster | Insect | Insect | 76.92 | 83.65 | 24 | 0.2624 | 780.0 |
| Takifugu_rubripes | Fish | Fish | 71.43 | 82.86 | 30 | 0.3365 | 430.0 |
| Saccharomyces_cerevisiae | Fungus | Fungus | 62.5 | 75.0 | 39 | 0.47 | 1100.0 |

## 5. Identity vs. Divergence Time & Molecular-Clock Test
Divergence times are approximate MRCA node ages (Mya) consistent with `config/species.csv` and TimeTree consensus estimates. **Caveat:** pairwise comparisons are **not statistically independent** due to shared phylogenetic ancestry, so p-values are descriptive only.

- **All 136 pairs — Identity vs. Time:** Spearman rho = **-0.8318** (p = 4.63e-36)
- **All 136 pairs — Poisson Distance vs. Time:** Spearman rho = **0.8318** (p = 4.63e-36); Pearson R² = **0.7516**, slope = **0.000387** subs/site/Myr (**0.387** subs/site/Gyr)
- **Human-vs-X (n = 16) — Identity vs. Time:** Spearman rho = **-0.9575** (p = 6.00e-09)
- **Human-vs-X (n = 16) — Poisson Distance vs. Time:** Spearman rho = **0.9575** (p = 6.00e-09); Pearson R² = **0.8154**, slope = **0.000352** subs/site/Myr (**0.352** subs/site/Gyr)

### Top 3 Molecular-Clock Outliers (Human-vs-X Residuals)
| species | clade | mya_divergence | poisson_distance | predicted_distance | residual |
| --- | --- | --- | --- | --- | --- |
| Takifugu_rubripes | Fish | 430.0 | 0.3365 | 0.2003 | 0.1362 |
| Drosophila_melanogaster | Insect | 780.0 | 0.2624 | 0.3236 | -0.0612 |
| Macaca_mulatta | Mammal | 29.0 | 0.0096 | 0.0591 | -0.0495 |

## 6. Group-Level Summary
| relationship | n_pairs | mean_identity | median_identity | min_identity | max_identity | sd_identity | min_pair | max_pair |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| within_primates | 1 | 99.05 | 99.05 | 99.05 | 99.05 | 0.0 | Homo_sapiens vs Macaca_mulatta | Homo_sapiens vs Macaca_mulatta |
| within_rodents | 1 | 100.0 | 100.0 | 100.0 | 100.0 | 0.0 | Mus_musculus vs Rattus_norvegicus | Mus_musculus vs Rattus_norvegicus |
| primate_vs_rodent | 4 | 91.9 | 91.9 | 91.43 | 92.38 | 0.55 | Homo_sapiens vs Mus_musculus | Macaca_mulatta vs Mus_musculus |
| within_other_mammals | 22 | 93.38 | 94.29 | 88.57 | 97.14 | 2.93 | Homo_sapiens vs Equus_caballus | Mus_musculus vs Bos_taurus |
| within_birds | 1 | 94.29 | 94.29 | 94.29 | 94.29 | 0.0 | Gallus_gallus vs Taeniopygia_guttata | Gallus_gallus vs Taeniopygia_guttata |
| within_reptiles | 1 | 91.43 | 91.43 | 91.43 | 91.43 | 0.0 | Anolis_carolinensis vs Chrysemys_picta | Anolis_carolinensis vs Chrysemys_picta |
| within_fish | 1 | 75.0 | 75.0 | 75.0 | 75.0 | 0.0 | Danio_rerio vs Takifugu_rubripes | Danio_rerio vs Takifugu_rubripes |
| bird_vs_reptile | 4 | 92.14 | 91.9 | 91.43 | 93.33 | 0.91 | Taeniopygia_guttata vs Anolis_carolinensis | Gallus_gallus vs Anolis_carolinensis |
| mammal_vs_bird | 16 | 90.89 | 90.95 | 87.62 | 94.29 | 2.11 | Homo_sapiens vs Gallus_gallus | Mus_musculus vs Taeniopygia_guttata |
| mammal_vs_reptile | 16 | 88.81 | 89.05 | 84.76 | 91.43 | 2.16 | Homo_sapiens vs Chrysemys_picta | Mus_musculus vs Chrysemys_picta |
| mammal_vs_amphibian | 8 | 89.76 | 90.48 | 86.67 | 92.38 | 2.2 | Homo_sapiens vs Xenopus_tropicalis | Monodelphis_domestica vs Xenopus_tropicalis |
| bird_vs_amphibian | 2 | 88.1 | 88.1 | 87.62 | 88.57 | 0.67 | Gallus_gallus vs Xenopus_tropicalis | Taeniopygia_guttata vs Xenopus_tropicalis |
| reptile_vs_amphibian | 2 | 88.1 | 88.1 | 86.67 | 89.52 | 2.02 | Chrysemys_picta vs Xenopus_tropicalis | Anolis_carolinensis vs Xenopus_tropicalis |
| mammal_vs_fish | 16 | 80.78 | 79.45 | 71.43 | 91.35 | 8.32 | Homo_sapiens vs Takifugu_rubripes | Mus_musculus vs Danio_rerio |
| bird_vs_fish | 4 | 83.76 | 83.76 | 76.19 | 91.35 | 7.72 | Taeniopygia_guttata vs Takifugu_rubripes | Taeniopygia_guttata vs Danio_rerio |
| reptile_vs_fish | 4 | 82.09 | 81.85 | 75.24 | 89.42 | 7.41 | Anolis_carolinensis vs Takifugu_rubripes | Anolis_carolinensis vs Danio_rerio |
| amphibian_vs_fish | 2 | 82.34 | 82.34 | 73.33 | 91.35 | 12.74 | Xenopus_tropicalis vs Takifugu_rubripes | Xenopus_tropicalis vs Danio_rerio |
| vertebrate_vs_drosophila | 15 | 80.9 | 81.73 | 69.23 | 84.62 | 3.95 | Takifugu_rubripes vs Drosophila_melanogaster | Mus_musculus vs Drosophila_melanogaster |
| vertebrate_vs_yeast | 15 | 61.41 | 61.54 | 58.65 | 62.5 | 1.35 | Chrysemys_picta vs Saccharomyces_cerevisiae | Homo_sapiens vs Saccharomyces_cerevisiae |
| drosophila_vs_yeast | 1 | 61.11 | 61.11 | 61.11 | 61.11 | 0.0 | Drosophila_melanogaster vs Saccharomyces_cerevisiae | Drosophila_melanogaster vs Saccharomyces_cerevisiae |
| mammal_vs_mammal_all | 28 | 93.61 | 93.81 | 88.57 | 100.0 | 3.13 | Homo_sapiens vs Equus_caballus | Mus_musculus vs Rattus_norvegicus |
