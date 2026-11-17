# AGENTS.md: Master Context for the Cytochrome c Project

> **Read this file completely before doing anything.** It is the single source of truth for project goals, conventions, and workflow rules. Each phase has its own file in `phases/`. Never start a phase without reading its file in full.
> (`CLAUDE.md` is an identical copy of this file. If you edit one, edit both.)

---

## 1. Project summary

**Title:** Comparative Sequence and Phylogenetic Analysis of Cytochrome c Across Vertebrate Species

**Team:** Maaz Mahmood (23L-0728), Qasim Shakil (23I-0687)
**Instructor:** Ms. Rushda Muneer
**Institution:** National University of Computer and Emerging Sciences (FAST-NUCES), Department of Computer Science, Lahore, Pakistan

**What the project does.** Cytochrome c (gene *CYCS* in vertebrates) is a small (~104-105 aa), heme-containing, mitochondrial electron-transfer protein. It is extremely conserved across eukaryotes. We collect its protein sequence from ~15 vertebrates plus 2 non-vertebrate outgroups, align the sequences, quantify per-position conservation, compute pairwise similarity, build phylogenetic trees, and test whether the sequence-derived relationships match the accepted species tree.

**This is a computational-only project.** No wet-lab work. No whole-genome analysis.

### Research questions (every phase must serve at least one)

| ID | Question |
|----|----------|
| RQ1 | How conserved is cytochrome c across the selected vertebrate species? |
| RQ2 | Which amino-acid positions or regions are highly conserved? |
| RQ3 | How does sequence similarity vary between closely and distantly related species? |
| RQ4 | Can cytochrome c sequence information be used to infer meaningful phylogenetic relationships? |

### Objectives (from the proposal)

1. Collect cytochrome c protein sequences from selected vertebrates.
2. Perform multiple sequence alignment.
3. Identify conserved and variable positions/regions.
4. Calculate and compare pairwise sequence similarity.
5. Construct a phylogenetic tree.
6. Analyze whether sequence-based relationships match known evolutionary relationships.

### Hypothesis to test (not assume)

Closely related species show higher sequence similarity and cluster together on the tree. **Test this honestly.** Cytochrome c is so conserved that many mammal sequences will be identical or nearly so, which limits tree resolution. If the data does not support the hypothesis in places, report that. Do not tune parameters to force a "nice" result.

---

## 2. Phase roadmap

Phases must be executed **in order**. Phase N may only start when Phase N-1's "Done when" checklist is fully satisfied.

| Phase | File | Goal |
|-------|------|------|
| 0 | `phases/PHASE_0_setup.md` | Environment, repo scaffolding, config, logging, tool installation |
| 1 | `phases/PHASE_1_data_collection.md` | Species selection and sequence retrieval (UniProt/NCBI) |
| 2 | `phases/PHASE_2_sequence_preparation.md` | QC, cleaning, standardization of the dataset |
| 3 | `phases/PHASE_3_alignment.md` | Multiple sequence alignment (Clustal Omega + MUSCLE) and QC |
| 4 | `phases/PHASE_4_conservation.md` | Per-position conservation, conserved regions, sequence logo |
| 5 | `phases/PHASE_5_similarity.md` | Pairwise identity/similarity matrices and heatmaps |
| 6 | `phases/PHASE_6_phylogenetics.md` | NJ and ML trees, bootstrap, rooting, visualization exports |
| 7 | `phases/PHASE_7_interpretation.md` | Biological interpretation, comparison to the known species tree |
| 8 | `phases/PHASE_8_reporting.md` | Final figures, report, README, reproducibility check |

Track progress in `PROGRESS.md` (created in Phase 0). After finishing a phase, tick its boxes and write a 3-6 line summary of what was done, key numbers, and any deviations.

---

## 3. Repository layout (created in Phase 0)

```
cytc-phylo/
├── CLAUDE.md / AGENTS.md        # this file
├── README.md                    # written in Phase 8
├── PROGRESS.md                  # phase tracker, updated after each phase
├── environment.yml              # conda environment definition
├── Makefile                     # one target per phase + `all`
├── config/
│   ├── config.yaml              # ALL tunable parameters live here
│   └── species.csv              # curated species list (Phase 1)
├── phases/                      # the phase instruction files
├── data/
│   ├── raw/                     # untouched downloads (never edit by hand)
│   ├── interim/                 # intermediate cleaned files
│   └── processed/               # final curated dataset used downstream
├── results/
│   ├── alignments/
│   ├── conservation/
│   ├── similarity/
│   ├── trees/
│   ├── interpretation/
│   ├── figures/                 # PNG (300 dpi) + SVG/PDF
│   └── tables/                  # CSV tables for the report
├── src/cytc/                    # importable helper modules
│   ├── __init__.py
│   └── utils.py
├── scripts/                     # one numbered script per phase step
├── tests/                       # pytest sanity tests
├── logs/                        # run logs (one per script run)
└── report/                      # final report (Phase 8)
```

---

## 4. Conventions (apply to every phase)

### Code
- **Python 3.11**, Biopython, pandas, numpy, scipy, matplotlib, seaborn, dendropy.
- Scripts live in `scripts/` and are named `NN_short_name.py` (e.g. `10_fetch_uniprot.py`). The number prefix matches the phase (Phase 1 -> 10-19, Phase 2 -> 20-29, and so on).
- Every script: has a module docstring, uses `argparse` or reads `config/config.yaml`, uses `logging` through `cytc.utils.get_logger` (no bare `print` except in tiny CLI summaries), and is runnable from the repo root as `python scripts/NN_name.py`.
- **No hard-coded absolute paths.** All paths come from `config.yaml` and are resolved relative to the repo root.
- **No magic numbers.** Thresholds (e.g. conservation cutoffs, bootstrap replicates) live in `config.yaml`.
- Set `random_seed` from config wherever randomness is involved (bootstraps, tie-breaking).
- Type hints on functions. Short docstrings. Functions under about 50 lines where reasonable.
- Idempotent: re-running a script overwrites its own outputs and does not corrupt others.

### Data integrity
- `data/raw/` is **read-only after download**. Never modify raw files. Every transformation writes to `data/interim/` or `data/processed/`.
- **Never fabricate or "recall from memory" sequences.** All sequences must be downloaded from UniProt/NCBI and carry an accession. Any accession written in these phase files is a *candidate* to be verified programmatically, not a trusted fact.
- Record provenance: accession, database, download date, release/version, URL used.
- **Do not deduplicate across species.** Identical sequences from different species are real biological signal (e.g. human and chimp cytochrome c are identical). Only remove duplicate *records for the same species*.

### Outputs
- Figures: PNG at 300 dpi **and** SVG or PDF. Readable fonts (>= 9 pt), labelled axes, titles, and colour-blind-safe palettes (e.g. `viridis`, `cividis`, seaborn `colorblind`).
- Tables: CSV with a header row, UTF-8, no index column unless meaningful.
- FASTA headers: `>Genus_species|ACCESSION` (see Phase 2). Species names use underscores, no spaces.
- Alignment files: FASTA format, one sequence per record, unwrapped or wrapped at 60 columns (be consistent).

### Numbering convention (critical, easy to get wrong)
Cytochrome c residue numbers appear in two systems in the literature:
- **UniProt numbering** counts the initiator Met as residue 1 (human CYCS: 105 aa).
- **Classical (horse-heart) numbering** drops the initiator Met, so classical position = UniProt position - 1 for vertebrates (Cys14/Cys17/His18/Met80 are classical numbers).
- Yeast iso-1 cytochrome c has an N-terminal extension, so its numbering is offset differently.

Whenever you report a position, **state which numbering you use**. Default: report **both** UniProt-human numbering and classical numbering, and also the alignment column index (1-based).

### Logging and errors
- Log to `logs/<script_name>_<YYYYmmdd_HHMMSS>.log` and to stdout.
- Fail loudly: raise on missing inputs, unexpected record counts, invalid residues. Do not silently skip.
- Each script ends with a one-line `SUMMARY:` log message stating counts and output paths.

### Testing
- Each phase adds at least 2-3 `pytest` tests in `tests/` (file exists, record count correct, no invalid characters, matrix symmetric, tree has the expected leaves, and so on). Run `pytest -q` before ticking a phase as done.

### Git (if the repo is a git repo)
- Commit at the end of each phase with message `Phase N: <summary>`.
- Do not commit `data/raw/` if large (these files are tiny, so committing them is fine and helps reproducibility). Do not commit `logs/` or conda envs. Add a `.gitignore`.

---

## 5. How to work (workflow rules for the agent)

1. **Read the phase file fully** before writing code. Restate the goal and inputs in your own words in the first log message.
2. **Check prerequisites.** Verify the previous phase's outputs exist and pass their validation checks. If not, stop and fix that phase first.
3. **Implement step by step**, running each script as you go. Inspect outputs (print heads of tables, check file sizes), not just exit codes.
4. **Validate** using the phase's "Validation checks" section. Every check must pass or be explicitly justified in `PROGRESS.md`.
5. **Handle deviations transparently.** If a candidate accession is wrong, a species has no good sequence, or a tool version behaves differently, choose a documented alternative and record it in `PROGRESS.md` under "Deviations".
6. **Do not skip ahead** and do not rewrite earlier phases' outputs unless fixing a bug (then re-run all downstream phases).
7. **Ask the user** only if blocked on something that cannot be resolved from these files (e.g. network access denied, license issue). Otherwise make a reasonable documented decision and continue.
8. **Be honest in interpretation.** Phase 7 must report where the tree disagrees with the known species tree and why.

---

## 6. Biological background cheat-sheet (for sanity checks)

- Vertebrate somatic cytochrome c: ~104-105 aa (UniProt length usually 105 including Met1).
- **Heme c attachment motif:** `C-x-x-C-H` (vertebrates: `CxQCH`). The two Cys (classical Cys14 and Cys17) form thioether bonds with the heme vinyl groups; His18 is the proximal (axial) heme ligand.
- **Met80** (classical) is the distal axial heme ligand. Its loop (about residues 70-85) is a key conserved region.
- Conserved surface **lysine patch** around the heme crevice mediates binding to cytochrome *c* oxidase, cytochrome *bc1*, and Apaf-1.
- Other classically conserved residues include Trp59, Tyr67, Asn52, Thr78, Gly41 and Phe82; **verify against the alignment, do not assume**.
- Some species have more than one cytochrome c paralog (e.g. a testis-specific one in mouse/rat, and several in Drosophila). **Use only the somatic/main cytochrome c ortholog** (see Phase 1).
- Yeast (*S. cerevisiae*) iso-1 cytochrome c (CYC1) has an N-terminal extension and some internal differences. It is used here as a non-vertebrate **outgroup**.
- Expected sanity trend: identity to human decreases with evolutionary distance (primates > rodents > other mammals > birds/reptiles > amphibians > fish > insects > yeast).

---

## 7. Definition of "done" for the entire project

- [ ] All nine phases ticked in `PROGRESS.md`.
- [ ] `make all` reproduces every table and figure from `data/raw/` in a fresh environment.
- [ ] `pytest -q` passes.
- [ ] All 8 proposal deliverables exist (see Phase 8 checklist).
- [ ] All 4 research questions are answered with specific numbers and figures.
- [ ] Limitations are stated honestly in the report.

