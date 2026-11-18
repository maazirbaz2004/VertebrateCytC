# Comparative Sequence and Phylogenetic Analysis of Cytochrome c Across Vertebrate Species

[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-3110/)
[![Bioinformatics](https://img.shields.io/badge/domain-Bioinformatics-success.svg)](#)
[![Build Status](https://img.shields.io/badge/status-in--progress-orange.svg)](#)

## 1. Project Overview

Cytochrome *c* (encoded by the *CYCS* gene in vertebrates) is a small (~104–105 amino acids), nuclear-encoded, mitochondrial electron-transfer protein that is essential for cellular respiration and apoptosis. Due to strong functional and structural constraints, cytochrome *c* is remarkably conserved across eukaryotes, making it a classical model system for molecular evolution and comparative sequence analysis.

This computational biology project collects, curates, aligns, and analyzes somatic cytochrome *c* protein sequences across **15 vertebrate species** spanning major taxonomic lineages, alongside **2 non-vertebrate outgroups** (17 species total). We quantify position-specific evolutionary conservation, compute pairwise sequence similarity, infer phylogenetic trees (Neighbor-Joining and Maximum Likelihood), and test whether sequence-derived relationships match accepted evolutionary history.

### Academic Context
- **Institution:** National University of Computer and Emerging Sciences (FAST-NUCES), Lahore, Pakistan
- **Department:** Department of Computer Science
- **Course:** Bioinformatics
- **Instructor:** Ms. Rushda Muneer
- **Team Members:**
  - Maaz Mahmood (23L-0728)
  - Qasim Shakil (23I-0687)

---

## 2. Research Questions

| ID | Research Question |
|:---|:---|
| **RQ1** | How conserved is cytochrome *c* across the selected vertebrate species? |
| **RQ2** | Which specific amino-acid positions, motifs, or structural regions are invariant or highly conserved? |
| **RQ3** | How does sequence similarity correlate between closely related vs. distantly related vertebrate species? |
| **RQ4** | Can cytochrome *c* primary sequence data alone infer meaningful and robust phylogenetic relationships? |

---

## 3. Dataset & Species Selection

The dataset consists of **15 vertebrate species** across diverse evolutionary lineages plus **2 non-vertebrate outgroups** for tree rooting:

| Taxonomic Group | Species | Common Name | Taxon ID | Accession | Database / Status |
|:---|:---|:---|:---|:---|:---|
| **Primate** | *Homo sapiens* | Human | 9606 | `P99999` | UniProt (Swiss-Prot) |
| **Primate** | *Macaca mulatta* | Rhesus macaque | 9544 | `P00002` | UniProt (Swiss-Prot) |
| **Rodent** | *Mus musculus* | House mouse | 10090 | `P62897` | UniProt (Swiss-Prot) |
| **Rodent** | *Rattus norvegicus* | Brown rat | 10116 | `P62898` | UniProt (Swiss-Prot) |
| **Carnivore** | *Canis lupus familiaris* | Domestic dog | 9615 | `P00011` | UniProt (Swiss-Prot) |
| **Artiodactyl** | *Bos taurus* | Cattle | 9913 | `P62894` | UniProt (Swiss-Prot) |
| **Perissodactyl** | *Equus caballus* | Horse | 9796 | `P00004` | UniProt (Swiss-Prot) |
| **Marsupial** | *Monodelphis domestica* | Gray short-tailed opossum | 13616 | `A0A5F8HEK2` | UniProt (TrEMBL) |
| **Bird** | *Gallus gallus* | Chicken | 9031 | `P67881` | UniProt (Swiss-Prot) |
| **Bird** | *Taeniopygia guttata* | Zebra finch | 59729 | `B5FZX0` | UniProt (TrEMBL) |
| **Reptile** | *Anolis carolinensis* | Green anole | 28377 | `A0A803TKE1` | UniProt (TrEMBL) |
| **Reptile** | *Chrysemys picta* | Painted turtle | 8478 | `A0A8C3PCA0` | UniProt (TrEMBL) |
| **Amphibian** | *Xenopus tropicalis* | Western clawed frog | 8364 | `Q640U4` | UniProt (Swiss-Prot) |
| **Fish** | *Danio rerio* | Zebrafish | 7955 | `Q6IQM2` | UniProt (Swiss-Prot) |
| **Fish** | *Takifugu rubripes* | Torafugu | 31033 | `Q1KL06` | UniProt (Swiss-Prot) |
| **Insect (Outgroup)** | *Drosophila melanogaster* | Fruit fly | 7227 | `P84029` | UniProt (Swiss-Prot) |
| **Fungus (Outgroup)** | *Saccharomyces cerevisiae* | Baker's yeast | 559292 | `P00044` | UniProt (Swiss-Prot) |

> **Quality Control & Paralog Filtering:** All sequences are verified to possess the signature `C-x-x-C-H` heme-c covalent attachment motif and standard 20 amino acids. Non-somatic paralogs (e.g., testis-specific *Cyct* in rodents/amphibians, *Cyt-c-d* in fruit fly, and *cyc-B* in fugu) were strictly filtered out.

---

## 4. Repository Structure

```text
cytc-phylo/
├── CLAUDE.md / AGENTS.md        # Master protocol and single source of truth
├── README.md                    # Project documentation
├── PROGRESS.md                  # Milestone and phase execution tracker
├── environment.yml              # Conda environment definition
├── Makefile                     # Stage execution entry points
├── conftest.py                  # Pytest path configuration
├── config/
│   ├── config.yaml              # Centralized configuration (parameters, thresholds)
│   ├── species.csv              # Curated species table with accessions
│   └── species.csv.bak          # Pre-retrieval backup
├── phases/                      # Phase-by-phase protocol instructions
├── data/
│   ├── raw/                     # Untouched raw API downloads (JSONs, fasta, audit log)
│   │   ├── uniprot/             # Raw JSON responses per accession
│   │   ├── cytc_raw.fasta       # Verified sequence FASTA
│   │   └── candidates_all.csv   # Comprehensive audit log of all screened candidates
│   ├── interim/                 # Cleaned and standardized intermediate files (Phase 2)
│   └── processed/               # Final curated datasets
├── results/
│   ├── alignments/              # Clustal Omega and MUSCLE v5 alignments (Phase 3)
│   ├── conservation/            # Per-position scores, logos, invariant regions (Phase 4)
│   ├── similarity/              # Pairwise similarity & identity matrices (Phase 5)
│   ├── trees/                   # NJ and ML phylogenetic trees (Phase 6)
│   ├── interpretation/          # Biological concordance and species tree tests (Phase 7)
│   ├── figures/                 # Publication-ready plots (300 DPI PNG & SVG)
│   └── tables/                  # Data tables and provenance reports
│       └── dataset_provenance.csv
├── src/cytc/                    # Reusable Python helper package
│   ├── __init__.py
│   └── utils.py                 # Config loader, logger, and path resolution
├── scripts/                     # Standalone numbered phase scripts
│   ├── 00_check_env.py          # Phase 0 environment and toolchain smoke test
│   └── 10_select_and_fetch.py   # Phase 1 sequence retrieval and validation
├── tests/                       # Pytest verification suites
│   ├── test_phase0.py           # Phase 0 scaffolding tests
│   └── test_phase1.py           # Phase 1 sequence integrity and provenance tests
└── logs/                        # Execution run logs with timestamps
```

---

## 5. Environment & Toolchain Setup

The project relies on Python 3.11 and native command-line bioinformatics packages. On Windows systems, execution is run in **WSL2 (Ubuntu)** to ensure compatibility with bioconda CLI tools.

### Prerequisites
- Conda or Mamba (Miniforge recommended)
- Git

### Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/maazirbaz2004/VertebrateCytC.git
   cd VertebrateCytC
   ```

2. **Create and activate the environment:**
   ```bash
   mamba env create -f environment.yml
   conda activate cytc
   ```

3. **Verify the environment:**
   ```bash
   python scripts/00_check_env.py
   ```
   *Verified tools:*
   - `python` 3.11.18
   - `biopython` 1.88
   - `clustalo` 1.2.4 (Clustal Omega)
   - `muscle` 5.3 (MUSCLE v5)
   - `iqtree3` 3.1.4 (IQ-TREE)
   - `FastTree` 2.2.0
   - `trimal` 1.5.rev1

---

## 6. Reproducing the Workflow

The project is structured with a `Makefile` supporting step-by-step reproduction:

```bash
# Check environment and CLI tools
make check

# Phase 1: Retrieve and validate raw sequences
make phase1

# Run automated tests
make test
```

---

## 7. Project Roadmap & Status

- [x] **Phase 0: Environment Setup & Scaffolding** (Completed)
- [x] **Phase 1: Species Selection & Sequence Retrieval** (Completed — 17 species collected)
- [ ] **Phase 2: Sequence Preparation & Quality Control** (Standardized headers, length QC)
- [ ] **Phase 3: Multiple Sequence Alignment** (Clustal Omega vs. MUSCLE v5 comparison)
- [ ] **Phase 4: Conservation Analysis** (Information entropy, sequence logos, invariant residues)
- [ ] **Phase 5: Pairwise Sequence Similarity** (BLOSUM62 matrices, divergence heatmaps)
- [ ] **Phase 6: Phylogenetic Tree Construction** (NJ, Maximum Likelihood with IQ-TREE & FastTree)
- [ ] **Phase 7: Evolutionary Interpretation** (Concordance with NCBI taxonomy / TimeTree)
- [ ] **Phase 8: Final Report & Figure Generation** (Publication assets and synthesis)

---

## 8. License

This project is developed for academic and computational research purposes at FAST-NUCES.
