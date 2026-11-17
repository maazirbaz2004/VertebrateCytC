# PHASE 2: Sequence Preparation, QC and Curation

**Read `CLAUDE.md` first. Prerequisite: Phase 1 is done.**

## 0. Objective

Turn the raw download into a clean, standardized, fully validated **curated dataset** that every later phase consumes. Produce QC reports (lengths, amino-acid composition, motif presence, identical-sequence groups) and a single canonical FASTA.

**Serves:** Objective 1 (curated dataset deliverable); quality foundation for RQ1-RQ4.

## 1. Prerequisites

- `data/raw/cytc_raw.fasta`, `results/tables/dataset_provenance.csv` and `config/species.csv` exist and passed Phase 1 tests.

## 2. Inputs

| File | Used for |
|------|----------|
| `data/raw/cytc_raw.fasta` | sequences (headers `>Genus_species|ACC|source`) |
| `results/tables/dataset_provenance.csv` | metadata |
| `config/species.csv` | group, role, divergence time |
| `config/config.yaml` -> `sequence_prep` | thresholds |

## 3. Tasks

### Task 2.1: Parse and standardize (`scripts/20_prepare_sequences.py`)

1. Parse with `Bio.SeqIO`.
2. Upper-case sequences; strip whitespace, `*` stop symbols and gap characters.
3. **Standardize headers** to `>Genus_species|ACCESSION` (drop the `source` field from the FASTA header, as it stays in the provenance CSV). Species name must exactly match `species.csv` -> `species`. IQ-TREE and Biopython both choke on `|` in some contexts, so also produce a **tree-safe** FASTA where the header is just `Genus_species` (the ID used in all alignments and trees). Keep an `id_map.csv` mapping `tree_id` <-> `accession`.
   - Decision: **the canonical ID used downstream is the species name** (e.g. `Homo_sapiens`). It is unique by construction (one sequence per species).
4. Validate each sequence:
   - only characters in `allowed_alphabet`; fail if any ambiguous residues remain (`max_ambiguous_fraction` = 0);
   - contains the `C..CH` motif; record its start position (1-based, UniProt-style) and the exact 5-mer;
   - length in the plausible range.
5. **Duplicate handling (careful):**
   - *Within a species*: there must be exactly one record. If duplicates exist, stop and report.
   - *Across species*: **do NOT remove** identical sequences. Instead, detect them and write `results/tables/identical_sequence_groups.csv` (group id, member species). This is biologically informative (e.g. two primates with identical cytochrome c), and matters for Phase 6 (zero-length branches).
6. Sort records in a **fixed, biologically meaningful order** (by `group` order: Primate, Rodent, Carnivore, Artiodactyl, Perissodactyl, Marsupial, Bird, Reptile, Amphibian, Fish, then outgroups). This ordering is reused in every table and heatmap.

### Task 2.2: QC report

Compute and write these tables to `results/tables/`:

**`seq_qc.csv`**: one row per species:
`species, accession, group, role, length, n_term_residue, c_term_residue, heme_motif, heme_motif_start, mw_kda, isoelectric_point, gravy, n_lys, n_arg, n_his, n_cys, n_met, fraction_charged`
Use `Bio.SeqUtils.ProtParam.ProteinAnalysis` for `molecular_weight()`, `isoelectric_point()`, `gravy()`.

**`aa_composition.csv`**: species x 20 amino acids (percent composition), long and wide forms.

**`length_summary.csv`**: min / max / mean / median / mode of length; note how many vertebrates have the modal length and how the outgroups differ.

**Biological expectations to check and log (warn, not fail):**
- Vertebrate cytochrome c is strongly **basic** (pI around 9-10; many Lys). Flag any vertebrate with pI < 8.5.
- Vertebrate molecular weight about 11.5-12 kDa.
- Cysteine count is small (typically 2, the two heme-binding Cys) in most vertebrates. Flag outliers.
- Yeast iso-1 is longer than vertebrates due to an N-terminal extension. Log its length difference.

### Task 2.3: QC figures (`results/figures/`)

1. `fig_length_by_species.png/.svg`: horizontal bar chart of sequence length, bars coloured by `group`, outgroups hatched.
2. `fig_aa_composition_heatmap.png/.svg`: heatmap species x amino acid (%), rows in the fixed order. Annotate.
3. `fig_pI_vs_length.png/.svg`: scatter of pI vs. length, labelled by species.

### Task 2.4: Write the final curated dataset (`data/processed/`)

| File | Content |
|------|---------|
| `cytc_curated.fasta` | headers `>Genus_species|ACCESSION`, 60-col wrap, fixed order |
| `cytc_curated_treeid.fasta` | headers `>Genus_species` only (**used for alignment and phylogenetics**) |
| `cytc_vertebrates.fasta` | same, ingroup only (for vertebrate-only conservation in Phase 4) |
| `id_map.csv` | `tree_id, accession, group, role, taxon_id, order_index` |
| `dataset_manifest.json` | counts, SHA-256 of each FASTA, creation date, tool versions |

### Task 2.5: Tests `tests/test_phase2.py`

```python
import re, pandas as pd
from Bio import SeqIO
from cytc.utils import REPO_ROOT

P = REPO_ROOT / "data/processed"

def test_counts_and_unique_ids():
    recs = list(SeqIO.parse(P / "cytc_curated_treeid.fasta", "fasta"))
    assert 12 <= len(recs) <= 20
    assert len({r.id for r in recs}) == len(recs)

def test_alphabet_and_motif():
    for r in SeqIO.parse(P / "cytc_curated_treeid.fasta", "fasta"):
        s = str(r.seq)
        assert set(s) <= set("ACDEFGHIKLMNPQRSTVWY"), r.id
        assert re.search("C..CH", s), r.id

def test_order_matches_id_map():
    ids = [r.id for r in SeqIO.parse(P / "cytc_curated_treeid.fasta", "fasta")]
    m = pd.read_csv(P / "id_map.csv").sort_values("order_index")
    assert ids == list(m.tree_id)

def test_outgroups_present():
    m = pd.read_csv(P / "id_map.csv")
    assert {"Saccharomyces_cerevisiae", "Drosophila_melanogaster"} <= set(m.tree_id)
```

## 4. Expected outputs

| Path | Description |
|------|-------------|
| `data/processed/cytc_curated.fasta` | curated dataset (deliverable #1) |
| `data/processed/cytc_curated_treeid.fasta` | species-ID FASTA for downstream tools |
| `data/processed/cytc_vertebrates.fasta` | ingroup only |
| `data/processed/id_map.csv`, `dataset_manifest.json` | metadata and checksums |
| `results/tables/seq_qc.csv`, `aa_composition.csv`, `length_summary.csv`, `identical_sequence_groups.csv` | QC tables |
| `results/figures/fig_length_by_species.*`, `fig_aa_composition_heatmap.*`, `fig_pI_vs_length.*` | QC figures |

## 5. Validation checks

- [ ] One and only one sequence per species; IDs unique.
- [ ] 0 ambiguous or non-standard residues anywhere in `data/processed/`.
- [ ] Every sequence has the heme motif; start positions logged (expected about residue 14-17 in vertebrates, in UniProt numbering about 15-18).
- [ ] Vertebrate pI/MW/Cys expectations checked, and warnings explained.
- [ ] Identical-sequence groups table written (may be empty or may list several species; both are valid).
- [ ] No cross-species deduplication was applied (record count in processed equals record count in raw).
- [ ] Fixed species order is identical across FASTA, `id_map.csv`, and tables.
- [ ] `pytest -q tests/test_phase2.py` passes.

## 6. Common pitfalls

- **Silently dropping identical sequences.** This would delete human/primate (or mouse/rat) twins and corrupt every downstream statistic. Don't.
- **Header characters.** `|`, `:`, `(`, `)`, spaces and `,` break Newick parsing and PHYLIP. Use only `[A-Za-z0-9_]` in `tree_id`.
- **Initiator methionine.** Do not strip Met1. Numbering conversions are handled later, explicitly.
- **Sorting alphabetically by accident.** Use the fixed biological order everywhere.
- **Checksums.** Record SHA-256 so Phase 8 can prove the dataset did not change.

## 7. Done when

- [ ] Section 5 fully ticked.
- [ ] `make phase2` runs from a clean `data/processed/`.
- [ ] `PROGRESS.md` updated; Phase 2 committed.
