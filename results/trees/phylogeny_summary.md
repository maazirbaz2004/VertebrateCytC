# Phase 6 Phylogenetic Reconstruction Summary (RQ4)

## 1. Methods & Random Seeds
- Random seed: `42` (plus sensitivity seeds `1` and `2` for IQ-TREE).
- NJ: Biopython Neighbor-Joining with Poisson distance (1000 bootstrap replicates; 0 negative branch(es) clamped to 0) and BLOSUM62 cross-check.
- FastTree: approximate ML with `-lg -gamma` (SH-like local supports, 0–1).
- IQ-TREE (primary): ModelFinder Plus (`MFP`), `1000` UFBoot + `1000` SH-aLRT + `-bnni` + `-keep-ident`, rooted on `Saccharomyces_cerevisiae`.
- Polytomy-collapsed ML tree (`cytc_all.collapsed.nwk`): **7** internal branch(es) collapsed (length < 1e-6 or UFBoot < 50).
- Identical-sequence check: Mus_musculus & Rattus_norvegicus patristic distance = 0.000002 (resolution arbitrary)

## 2. IQ-TREE Model & Likelihood Summary
| run | alignment | best_model | log_likelihood | n_parameters | bic | tree_length | identical_seq_warnings |
| --- | --- | --- | --- | --- | --- | --- | --- |
| cytc_all | primary_all.fasta | Q.YEAST+G4 | -999.9432 | 32 | 2150.3018 | 2.4 | none noted in .iqtree |
| cytc_vert | primary_vert.fasta | DAYHOFF+G4 | -740.1563 | 28 | 1610.6235 | 0.9431 | none noted in .iqtree |
| cytc_all_lgg4 | primary_all.fasta | LG+G4 | -1001.548 | 32 | 2153.5114 | 2.3905 | none noted in .iqtree |
| cytc_all_muscle | all_muscle.fasta | DAYHOFF+G4 | -997.4133 | 32 | 2145.242 | 2.0899 | none noted in .iqtree |
| cytc_all_trimmed | primary_all_trimmed.fasta | Q.YEAST+G4 | -967.8166 | 32 | 2084.56 | 2.2125 | none noted in .iqtree |
| cytc_all_seed1 | primary_all.fasta | Q.YEAST+G4 | -999.9431 | 32 | 2150.3015 | 2.3989 | none noted in .iqtree |
| cytc_all_seed2 | primary_all.fasta | DAYHOFF+G4 | -998.8553 | 32 | 2148.1261 | 2.0944 | none noted in .iqtree |

## 3. Primary ML Clade Support (2/15 internal nodes supported at UFBoot >= 70 & SH-aLRT >= 80)
| node_id | n_leaves | leaves | sh_alrt | ufboot | branch_length | supported |
| --- | --- | --- | --- | --- | --- | --- |
| N01 | 16 | Anolis_carolinensis;Bos_taurus;Canis_lupus_familiaris;Chrysemys_picta;Danio_rerio;Drosophila_melanogaster;Equus_caballus;Gallus_gallus;Homo_sapiens;Macaca_mulatta;Monodelphis_domestica;Mus_musculus;Rattus_norvegicus;Taeniopygia_guttata;Takifugu_rubripes;Xenopus_tropicalis | nan | nan | 1.05954 | False |
| N02 | 15 | Anolis_carolinensis;Bos_taurus;Canis_lupus_familiaris;Chrysemys_picta;Danio_rerio;Equus_caballus;Gallus_gallus;Homo_sapiens;Macaca_mulatta;Monodelphis_domestica;Mus_musculus;Rattus_norvegicus;Taeniopygia_guttata;Takifugu_rubripes;Xenopus_tropicalis | 89.9 | 70.0 | 0.138949 | True |
| N03 | 5 | Anolis_carolinensis;Chrysemys_picta;Gallus_gallus;Taeniopygia_guttata;Takifugu_rubripes | 87.2 | 49.0 | 0.048656 | False |
| N04 | 4 | Anolis_carolinensis;Gallus_gallus;Taeniopygia_guttata;Takifugu_rubripes | 28.2 | 25.0 | 0.022317 | False |
| N05 | 3 | Anolis_carolinensis;Gallus_gallus;Takifugu_rubripes | 0.0 | 17.0 | 2e-06 | False |
| N06 | 2 | Gallus_gallus;Takifugu_rubripes | 87.9 | 43.0 | 0.023059 | False |
| N07 | 10 | Bos_taurus;Canis_lupus_familiaris;Danio_rerio;Equus_caballus;Homo_sapiens;Macaca_mulatta;Monodelphis_domestica;Mus_musculus;Rattus_norvegicus;Xenopus_tropicalis | 0.0 | 32.0 | 2e-06 | False |
| N08 | 2 | Mus_musculus;Rattus_norvegicus | 0.0 | 94.0 | 1e-06 | False |
| N09 | 8 | Bos_taurus;Canis_lupus_familiaris;Danio_rerio;Equus_caballus;Homo_sapiens;Macaca_mulatta;Monodelphis_domestica;Xenopus_tropicalis | 86.4 | 33.0 | 0.011087 | False |
| N10 | 3 | Bos_taurus;Canis_lupus_familiaris;Equus_caballus | 79.9 | 67.0 | 0.010912 | False |
| N11 | 2 | Bos_taurus;Equus_caballus | 78.7 | 70.0 | 0.011645 | False |
| N12 | 5 | Danio_rerio;Homo_sapiens;Macaca_mulatta;Monodelphis_domestica;Xenopus_tropicalis | 0.0 | 9.0 | 2e-06 | False |
| N13 | 2 | Homo_sapiens;Macaca_mulatta | 98.9 | 100.0 | 0.086692 | True |
| N14 | 3 | Danio_rerio;Monodelphis_domestica;Xenopus_tropicalis | 84.6 | 56.0 | 0.022679 | False |
| N15 | 2 | Danio_rerio;Xenopus_tropicalis | 84.6 | 66.0 | 0.027772 | False |

## 4. Sensitivity Analysis (Unrooted Robinson-Foulds Distance vs. Primary ML Tree; 0/9 identical topology)
| variant | n_taxa | rf_distance_to_primary | max_rf | normalized_rf | identical_topology | notes |
| --- | --- | --- | --- | --- | --- | --- |
| nj_poisson | 17 | 18 | 28 | 0.6429 | False | 9 differing split(s) |
| nj_blosum62 | 17 | 18 | 28 | 0.6429 | False | 9 differing split(s) |
| fasttree_lg_gamma | 17 | 12 | 28 | 0.4286 | False | 6 differing split(s) |
| iqtree_vert_only | 15 | 4 | 24 | 0.1667 | False | compared on 15-vertebrate pruned namespace |
| iqtree_fixed_lgg4 | 17 | 6 | 28 | 0.2143 | False | 3 differing split(s) |
| iqtree_muscle | 17 | 12 | 28 | 0.4286 | False | 6 differing split(s) |
| iqtree_trimmed | 17 | 2 | 28 | 0.0714 | False | 1 differing split(s) |
| iqtree_seed1 | 17 | 2 | 28 | 0.0714 | False | 1 differing split(s) |
| iqtree_seed2 | 17 | 2 | 28 | 0.0714 | False | 1 differing split(s) |

## 5. Stable vs. Unstable Clades Across All 17-Taxon Variants
### Stable clades (recovered in 100% of 17-taxon variants)
- `{Bos_taurus, Equus_caballus}`
- `{Mus_musculus, Rattus_norvegicus}`
- `{Danio_rerio, Xenopus_tropicalis}`
- `{Bos_taurus, Canis_lupus_familiaris, Equus_caballus}`

### Unstable primary-tree clades (absent in >= 1 sensitivity variant)
- `(8/9) {Homo_sapiens, Macaca_mulatta}`
- `(7/9) {Danio_rerio, Monodelphis_domestica, Xenopus_tropicalis}`
- `(6/9) {Gallus_gallus, Takifugu_rubripes}`
- `(6/9) {Bos_taurus, Canis_lupus_familiaris, Danio_rerio, Equus_caballus, Homo_sapiens, Macaca_mulatta, Monodelphis_domestica, Xenopus_tropicalis}`
- `(5/9) {Drosophila_melanogaster, Saccharomyces_cerevisiae}`
- `(5/9) {Anolis_carolinensis, Gallus_gallus, Takifugu_rubripes}`
- `(5/9) {Anolis_carolinensis, Gallus_gallus, Taeniopygia_guttata, Takifugu_rubripes}`
- `(5/9) {Anolis_carolinensis, Chrysemys_picta, Gallus_gallus, Taeniopygia_guttata, Takifugu_rubripes}`
- `(5/9) {Anolis_carolinensis, Chrysemys_picta, Drosophila_melanogaster, Gallus_gallus, Saccharomyces_cerevisiae, Taeniopygia_guttata, Takifugu_rubripes}`
- `(2/9) {Danio_rerio, Homo_sapiens, Macaca_mulatta, Monodelphis_domestica, Xenopus_tropicalis}`
