# PHASE 8: Final Figures, Report, README and Reproducibility Check

**Read `CLAUDE.md` first. Prerequisite: Phase 7 is done.**

## 0. Objective

Finalize the project for submission: consolidate publication-quality figures and tables, assemble the final research report (Markdown source, with optional DOCX/PDF export), write the repository README, run a **from-scratch reproducibility check**, and verify every proposal deliverable against a checklist.

**Serves:** all objectives; delivers the "final research report documenting the methodology and results".

## 1. Prerequisites

- Phases 0-7 are ticked in `PROGRESS.md` and their tests pass.
- `results/interpretation/research_questions_answers.md` and `limitations.md` exist.

## 2. Inputs

Everything under `results/`, `data/processed/`, `config/`, `PROGRESS.md`, and the original proposal's structure (sections below).

## 3. Tasks

### Task 8.1: Figure and table consolidation (`scripts/80_make_report_assets.py`)

1. Create `report/figures/` and copy (or regenerate at the config DPI) the **final selected figures** with stable report names:

| Report figure | Source | Caption (draft) |
|---------------|--------|-----------------|
| Fig 1: Workflow | new: flowchart (matplotlib or Mermaid-rendered) of Data Collection -> Preparation -> MSA -> Conservation/Similarity -> Phylogenetics -> Interpretation | Overview of the analysis pipeline |
| Fig 2: Dataset overview | `fig_length_by_species` | Sequence lengths by species and clade |
| Fig 3: Alignment overview | `fig_alignment_overview` | Multiple sequence alignment coloured by residue class |
| Fig 4: Conservation profile | `fig_conservation_profile` | Per-position conservation; heme motif and Met80 marked |
| Fig 5: Sequence logo | `fig_sequence_logo_vert` | Logo of vertebrate cytochrome c |
| Fig 6: Identity heatmap | `fig_identity_heatmap` | Pairwise identity |
| Fig 7: Identity to human | `fig_identity_to_human` | Identity vs. divergence |
| Fig 8: Identity vs. time | `fig_identity_vs_time` | Molecular-clock behaviour |
| Fig 9: ML tree | `fig_tree_ml_rectangular` | ML tree with support |
| Fig 10: NJ vs ML | `fig_tree_comparison_NJ_vs_ML` | Method comparison |
| Fig 11: Tree vs. reference | `fig_tree_vs_reference` | Agreement with the accepted species tree |
| Fig 12: Functional sites | `fig_functional_sites_on_conservation` | Conservation at known functional sites |

2. Create `report/tables/` with report-ready CSV and **Markdown** versions of:
   - T1 dataset provenance (`dataset_provenance`: species, accession, database, length);
   - T2 sequence QC summary;
   - T3 alignment statistics;
   - T4 conservation summary and conserved regions;
   - T5 identity to human; T6 group-level similarity summary;
   - T7 model selection and sensitivity summary; T8 clade recovery; T9 functional checklist results.
3. Figure quality rules: 300 dpi PNG plus vector (SVG or PDF); fonts >= 9 pt at final size; colour-blind-safe; every axis labelled with units; legends not covering data; titles optional if the caption carries it; consistent clade colours **across all figures** (define the palette once in `src/cytc/plotting.py` and import it everywhere. Refactor earlier scripts if they used different colours).
4. A manifest `report/asset_manifest.csv` listing each report asset, its source file, and creation timestamp.

### Task 8.2: Write the final report (`report/final_report.md`)

Structure (mirrors the proposal so reviewers can match it; length about 8-15 pages when rendered):

1. **Title page**: title, group members with roll numbers, instructor, institution, date.
2. **Abstract** (150-250 words): problem, data, methods, headline numbers, main conclusions, key limitation.
3. **Introduction**: cytochrome c's role (cellular respiration, electron transfer), why it is a good molecular-evolution marker, aim and research questions.
4. **Background / Related Work**: the four studies cited in the proposal (Fitch 1976; Syvanen et al. 1989; Grossman et al. 2011; the 2006 *Gene* paper on metazoan gene families) **after the citation check in 8.3**, and how this project builds on them.
5. **Problem Statement and Objectives**.
6. **Methodology**: data sources and selection criteria (with Table T1), sequence QC, alignment tools and versions, conservation metrics (definitions), similarity metrics (definitions), tree-building methods and parameters (including model, bootstrap, rooting), validation approach (reference tree, RF), software versions and seeds (pull from logs).
7. **Results**:
   7.1 Dataset; 7.2 Alignment; 7.3 Conservation (RQ1, RQ2); 7.4 Pairwise similarity (RQ3); 7.5 Phylogenetic trees; 7.6 Agreement with known evolutionary relationships (RQ4); 7.7 Functional interpretation.
8. **Discussion**: the four RQs answered; hypothesis check; relation to prior studies; what cytochrome c can and cannot resolve; the biological meaning of the conserved regions.
9. **Limitations** (from `limitations.md`) and **Future work** (additional genes, e.g. concatenation with COX1/CYTB; more taxa; Bayesian inference (MrBayes/BEAST); structure-based analysis; testing for positive selection with PAML/HyPhy on the nucleotide coding sequences).
10. **Conclusion**.
11. **References** (IEEE style, as in the proposal), plus the software tools cited (Clustal Omega, MUSCLE v5, IQ-TREE 2, ModelFinder, UFBoot, FastTree, Biopython, DendroPy, UniProt, NCBI) **with their actual publications. Look up and verify each citation; do not invent bibliographic details.** If a reference cannot be verified, mark it `[VERIFY]` rather than guessing.
12. **Appendix**: reproducibility instructions, full species table, full pairwise matrices, command log.

Writing rules:
- **Every number in the text must come from a file in `results/`.** Generate the text's numbers by reading the CSVs (a script may fill a template, or the agent writes it while looking at the tables). No invented values.
- Use cautious language for inference ("consistent with", "suggests"), explicit about support levels.
- Explicitly report negative or unexpected findings (unresolved clades, identical sequences, disagreements with the reference tree).
- Refer to every figure and table in the text, in order of appearance.
- State residue numbering convention at first use.

Also render the report if tools are available: `pandoc report/final_report.md -o report/final_report.docx --citeproc` and, optionally, PDF (`--pdf-engine=xelatex`). If pandoc is not installed, tell the user the Markdown is complete and can be converted.

### Task 8.3: Citation audit (important)

The proposal's reference list contains at least one entry that appears inconsistent with how it is described:
- **Reference [2]** (Syvanen, Hartman, Stevens, *J. Mol. Evol.* 28:536-544, 1989) is titled "Classical plant taxonomic ambiguities extend to the molecular level", which suggests a **plant** study, whereas the proposal describes it as an analysis of vertebrate cytochrome c. 

Task for the agent: attempt to verify **all four** references via web search or Crossref (`https://api.crossref.org/works?query.bibliographic=...`), compare title, journal, year, volume and pages, and the claim each is used for. Write `report/citation_audit.md` with a table `ref, claim in proposal, verified title/venue, matches claim?, action`. For mismatched ones, **propose** a replacement (searching for a genuine vertebrate cytochrome c phylogeny paper) but mark proposals as suggestions for the **user to approve**. **Do not silently change the group's references.** If there is no network access, list the unresolved references in the audit file and tell the user.

### Task 8.4: README (`README.md`)

Sections: project title and one-paragraph summary; team and course; **key results at a glance** (5-8 bullet numbers pulled from results files, with the main figure embedded); repository layout; quick start (`mamba env create -f environment.yml && conda activate cytc && make all`); detailed phase table linking to `phases/*.md`; outputs map (where each proposal deliverable lives); software versions; data sources and licensing notes (UniProt CC BY 4.0, NCBI public domain; cite them); limitations summary; citation / acknowledgement; contact.

### Task 8.5: Reproducibility check

1. **Clean-room run:** in a fresh clone/copy of the repo (or after `make clean-results` and deleting `data/interim`, `data/processed`, `logs`), create the environment from `environment.yml` and run `make all`. All phases must complete without manual intervention.
2. **Compare** outputs with the originals: identical FASTA checksums (`dataset_manifest.json`), identical alignments (if tools deterministic), identical conservation/similarity tables (byte-for-byte or within 1e-9), and trees with the **same topology** (RF = 0 vs. originals; branch lengths may differ trivially; bootstrap supports may vary if seeds were not honoured. Document).
3. Write `results/tables/reproducibility_check.csv`: `artifact, original_sha256, rerun_sha256, identical (Y/N), note`.
4. Pin versions: export `conda env export --no-builds > environment.lock.yml` and record tool versions in `report/software_versions.md`.
5. **Network-dependent steps** (Phase 1 downloads): the raw data are committed under `data/raw/`; provide a `make phase1-offline` target that skips downloading and verifies checksums, so the analysis is reproducible even if UniProt changes.
6. Run `pytest -q` (all phases) and `ruff`/`flake8` if installed (non-blocking).

### Task 8.6: Final deliverables checklist (map every proposal deliverable)

Write `report/deliverables_checklist.md`:

| # | Proposal deliverable | Location | Verified |
|---|---------------------|----------|----------|
| 1 | Curated dataset of cytochrome c sequences | `data/processed/cytc_curated.fasta`, `id_map.csv`, `report/tables/T1` | [ ] |
| 2 | Multiple sequence alignment | `results/alignments/primary_all.fasta`, `primary_vert.fasta` | [ ] |
| 3 | Conserved and variable region analysis | `results/conservation/regions.csv`, `summary.md` | [ ] |
| 4 | Pairwise sequence similarity results | `results/similarity/*_matrix.csv`, `pairwise_long.csv` | [ ] |
| 5 | Visualization of sequence conservation | `fig_conservation_profile`, `fig_sequence_logo_*` | [ ] |
| 6 | Phylogenetic tree | `results/trees/iqtree/cytc_all.treefile`, `fig_tree_ml_rectangular`, iTOL files | [ ] |
| 7 | Biological interpretation of the findings | `results/interpretation/*`, report sections 7-8 | [ ] |
| 8 | Final research report | `report/final_report.md` (+ `.docx` / `.pdf`) | [ ] |

Also map **objectives** (1-6) and **research questions** (RQ1-RQ4) to the report sections and figures that address them.

### Task 8.7: Package for submission

- `make package` -> `dist/cytc_project_submission.zip` containing: report (md/docx/pdf), `README.md`, `environment.yml`, `config/`, `data/processed/`, `results/` (minus huge intermediates), `scripts/`, `src/`, `tests/`, `phases/`, `PROGRESS.md`.
- Add the `package` target to the Makefile.

### Task 8.8: Tests `tests/test_phase8.py`

```python
from cytc.utils import REPO_ROOT

def test_report_exists_and_has_sections():
    txt = (REPO_ROOT / "report/final_report.md").read_text()
    for h in ["Abstract", "Introduction", "Methodology", "Results", "Discussion",
              "Limitations", "Conclusion", "References"]:
        assert h in txt, h

def test_all_figures_in_manifest_exist():
    import pandas as pd
    man = pd.read_csv(REPO_ROOT / "report/asset_manifest.csv")
    for p in man["report_path"]:
        assert (REPO_ROOT / p).exists(), p

def test_deliverables_checklist_complete():
    txt = (REPO_ROOT / "report/deliverables_checklist.md").read_text()
    assert txt.count("[x]") >= 8

def test_no_placeholder_text_left():
    txt = (REPO_ROOT / "report/final_report.md").read_text()
    for bad in ["TODO", "TBD", "XXX", "CHANGE_ME", "lorem"]:
        assert bad not in txt, bad
```

(`[VERIFY]` markers on unresolved references are allowed, but the final message to the user must list them.)

## 4. Expected outputs

| Path | Description |
|------|-------------|
| `report/final_report.md` (+ `.docx`, `.pdf` if possible) | **final research report** |
| `report/figures/`, `report/tables/`, `report/asset_manifest.csv` | report assets |
| `report/citation_audit.md`, `report/software_versions.md`, `report/deliverables_checklist.md` | audits |
| `README.md`, `environment.lock.yml` | repo documentation and pins |
| `results/tables/reproducibility_check.csv` | proof of reproducibility |
| `dist/cytc_project_submission.zip` | submission package |

## 5. Validation checks

- [ ] Clean-room `make all` succeeds; reproducibility table shows identical checksums or documented, explained differences.
- [ ] `pytest -q` (all phases) passes.
- [ ] All eight proposal deliverables exist and are ticked in the checklist.
- [ ] Every figure referenced in the report exists, with caption and in-text citation; same clade colours throughout.
- [ ] No number in the report lacks a source file in `results/`.
- [ ] Citation audit completed; unresolved or mismatched references are flagged to the user (especially reference [2]).
- [ ] Report states limitations and negative findings honestly.
- [ ] README quick start works verbatim in a fresh environment.
- [ ] No placeholder text remains.

## 6. Common pitfalls

- **Fabricated references.** Never generate bibliographic details from memory. Verify or mark `[VERIFY]`.
- **Numbers drifting** between report and results after re-running an earlier phase. Regenerate the report numbers after any upstream change.
- **Inconsistent colours or species order** between figures.
- **Reproducibility theatre.** Don't just re-run in the same dirty directory. Do it clean.
- **Over-long report.** Keep the main text focused. Move bulky tables to the Appendix.
- **Leaving agent notes** ("as the agent I...") in user-facing documents.

## 7. Done when

- [ ] Section 5 fully ticked.
- [ ] `PROGRESS.md` shows all phases complete, with deviations documented.
- [ ] Final commit `Phase 8: final report and package`, plus a final message to the user summarizing: headline findings, any deviations, unresolved `[VERIFY]` items, and where to find each deliverable.
