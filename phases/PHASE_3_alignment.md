# PHASE 3: Multiple Sequence Alignment and Alignment QC

**Read `CLAUDE.md` first. Prerequisite: Phase 2 is done.**

## 0. Objective

Produce a reliable multiple sequence alignment (MSA) of the curated cytochrome c sequences using **two independent aligners** (Clustal Omega and MUSCLE v5), quantify how well they agree, analyse gaps, decide on the **primary alignment** used downstream, and build an optional **trimmed** alignment for a sensitivity analysis.

**Serves:** Objective 2; foundation for RQ1, RQ2, RQ4.

## 1. Prerequisites

- `data/processed/cytc_curated_treeid.fasta` and `cytc_vertebrates.fasta` exist and pass Phase 2 tests.
- Tools `clustalo` and `muscle` (v5) work (Phase 0).

## 2. Inputs

- `data/processed/cytc_curated_treeid.fasta` (17 sequences: vertebrates + outgroups).
- `data/processed/cytc_vertebrates.fasta` (vertebrates only).
- `config.yaml` -> `alignment`, `tools`.

## 3. Tasks

### Task 3.1: Run the aligners (`scripts/30_align.py`)

Use `subprocess.run([...], check=True, capture_output=True, text=True)` (never `shell=True`), log the exact command line and tool version, and write alignments to `results/alignments/`.

Run **both aligners on two datasets**: `all` (vertebrates + outgroups) and `vert` (vertebrates only):

```bash
# Clustal Omega
clustalo -i data/processed/cytc_curated_treeid.fasta \
         -o results/alignments/all_clustalo.fasta \
         --seqtype Protein --outfmt fasta --force --threads 4 -v

# MUSCLE v5
muscle -align data/processed/cytc_curated_treeid.fasta \
       -output results/alignments/all_muscle.fasta -threads 4
```

and the same for `cytc_vertebrates.fasta` -> `vert_clustalo.fasta`, `vert_muscle.fasta`.

Notes:
- Clustal Omega also gives a percent identity matrix: `--percent-id --distmat-out=results/alignments/all_clustalo.distmat --full` (useful as a cross-check against Phase 5, not as the primary source).
- MUSCLE v5 does **not** preserve input order by default. **Re-order** all alignment records to the fixed species order from `id_map.csv` after reading them (do this for both tools so alignments are directly comparable).
- Rewrite alignment FASTA with fixed 60-column wrapping and fixed order via Biopython (`SeqIO.write`).
- Also write **Clustal-format (`.aln`) and PHYLIP-relaxed copies** of the primary alignment for the convenience of MEGA/iTOL users.

### Task 3.2: Basic alignment statistics

For every alignment, compute and save to `results/tables/alignment_stats.csv`:

| Metric | Definition |
|--------|------------|
| `n_seq`, `aln_length` | dimensions |
| `n_gap_columns` | columns containing >= 1 gap |
| `n_gapfree_columns` | columns with no gaps |
| `n_invariant_columns` | gap-free columns with a single residue |
| `pct_identity_columns` | invariant / aln_length |
| `mean_pairwise_identity` | mean over all pairs (identity over gap-free pairwise-aligned positions) |
| `total_gap_chars`, `n_gap_blocks` | total "-" characters and contiguous gap runs |
| `sum_of_pairs_blosum62` | SP score with BLOSUM62 and a simple gap penalty (document formula) |

### Task 3.3: Aligner-agreement analysis

Compare Clustal Omega and MUSCLE for each dataset:

1. **Column-level agreement:** map each alignment column to a tuple of "residue index per sequence" (or `None` for gaps). A column is *shared* if the identical tuple exists in both alignments. Report **% shared columns** (the standard "column score" / TC score).
2. **Pairwise-residue agreement (SP-style):** fraction of aligned residue *pairs* found in both alignments.
3. **Disagreement map:** list the alignment regions (in primary-alignment coordinates) where the two disagree; expected to be mainly around the N-terminus (yeast extension) and a few loop positions.
4. If agreement is >= 95% for the vertebrate-only alignment, the alignment is considered robust. Otherwise, inspect and describe the disputed regions in `results/alignments/alignment_notes.md`.

Implementation sketch for column tuples:

```python
def column_signatures(aln) -> list[tuple]:
    """Each column -> tuple of per-sequence residue indices (None for gaps)."""
    counters = [0] * len(aln)
    sigs = []
    for col in range(aln.get_alignment_length()):
        sig = []
        for i, rec in enumerate(aln):
            if rec.seq[col] == "-":
                sig.append(None)
            else:
                sig.append(counters[i]); counters[i] += 1
        sigs.append(tuple(sig))
    return sigs
```

A column counts as agreeing only if it contains **at least 2 residues** (all-gap and single-residue columns are uninformative). Document this.

### Task 3.4: Gap analysis

- List every gap block: alignment coordinates, length, which sequences carry the gap.
- Check the expected pattern: the **yeast N-terminal extension** (yeast has residues the vertebrates lack) and perhaps one-residue indels in Drosophila or fish. Vertebrate-only alignment should have **few or zero internal gaps**, since cytochrome c is nearly indel-free in vertebrates.
- Output `results/tables/gap_blocks.csv`: `alignment, start, end, length, sequences_with_gap, sequences_without_gap`.
- Figure: `fig_alignment_overview.png/.svg`, a colour-coded overview (rows = species in fixed order, columns = alignment positions; colour by amino-acid class, white for gaps). Build it with matplotlib `imshow` on an integer-encoded matrix, with a legend by physico-chemical class (hydrophobic, polar, positive, negative, aromatic, special [G/P/C]). Label x-axis ticks every 10 columns.

### Task 3.5: Choose the primary alignment (decision rule, record it!)

Primary = `alignment.primary` from config (default **Clustal Omega on the `all` dataset**). Switch to MUSCLE only if it has a clearly higher SP score *and* fewer implausible gaps (e.g. gaps inside the conserved heme motif or Met80 loop). Whatever you choose, write the decision and justification into `results/alignments/alignment_notes.md` and copy the chosen files to canonical names:

```
results/alignments/primary_all.fasta       # vertebrates + outgroups
results/alignments/primary_vert.fasta      # vertebrates only
```

**Sanity check that must pass:** in the primary alignment the heme motif `C..CH` occupies the **same alignment columns for all 17 sequences** (with no gaps), and the Met80-containing region is aligned without gaps in vertebrates.

### Task 3.6: Optional trimmed alignment (sensitivity analysis only)

Run trimAl: `trimal -in primary_all.fasta -out primary_all_trimmed.fasta -gt 0.8 -colnumbering > trimmed_columns.txt`. This removes poorly supported columns (mainly the yeast extension). It is used **only** in Phase 6 to test whether the tree topology changes, never as the primary alignment, because the conserved protein is short and trimming discards information.

### Task 3.7: Tests `tests/test_phase3.py`

```python
from Bio import AlignIO
from cytc.utils import REPO_ROOT
import re

A = REPO_ROOT / "results/alignments"

def test_equal_length_rows():
    aln = AlignIO.read(A / "primary_all.fasta", "fasta")
    assert len({len(r.seq) for r in aln}) == 1 and len(aln) >= 12

def test_ungapped_sequences_match_curated():
    from Bio import SeqIO
    cur = {r.id: str(r.seq) for r in SeqIO.parse(REPO_ROOT / "data/processed/cytc_curated_treeid.fasta", "fasta")}
    for r in AlignIO.read(A / "primary_all.fasta", "fasta"):
        assert str(r.seq).replace("-", "") == cur[r.id], r.id   # alignment must not alter residues

def test_heme_motif_columns_consistent():
    aln = AlignIO.read(A / "primary_all.fasta", "fasta")
    starts = set()
    for r in aln:
        m = re.search("C-*.-*.-*C-*H", str(r.seq))  # tolerate gaps just to detect problems
        assert m, r.id
        assert "-" not in m.group(0), f"gap inside heme motif in {r.id}"
        starts.add(m.start())
    assert len(starts) == 1
```

## 4. Expected outputs

| Path | Description |
|------|-------------|
| `results/alignments/{all,vert}_{clustalo,muscle}.fasta` | four raw alignments (fixed order) |
| `results/alignments/primary_all.fasta`, `primary_vert.fasta` | chosen primary alignments (deliverable #2) |
| `results/alignments/primary_all.aln`, `.phy` | alternate formats |
| `results/alignments/primary_all_trimmed.fasta`, `trimmed_columns.txt` | sensitivity-analysis alignment |
| `results/alignments/alignment_notes.md` | decision log, tool versions, commands, disagreements |
| `results/tables/alignment_stats.csv`, `gap_blocks.csv`, `aligner_agreement.csv` | tables |
| `results/figures/fig_alignment_overview.*` | alignment overview figure |

## 5. Validation checks

- [ ] Both aligners ran on both datasets; commands and versions logged.
- [ ] Residues in each alignment row, after removing gaps, are **identical** to the curated sequences (no tool altered residues).
- [ ] Column agreement between aligners reported; >= 95% on the vertebrate alignment, or disputed regions are explained.
- [ ] Heme motif aligned in the same columns, gap-free, for all sequences.
- [ ] Vertebrate-only alignment has no gap inside the motif or Met80 region.
- [ ] Primary alignment chosen with a written justification.
- [ ] `pytest -q tests/test_phase3.py` passes.

## 6. Common pitfalls

- **MUSCLE v3 syntax** (`-in/-out`) fails on v5, which uses `-align/-output`.
- **Record order changes.** MUSCLE v5 and Clustal Omega may reorder; always reorder to the fixed order (otherwise heatmaps and tables misalign silently).
- **Lower-case residues** from some tools. Normalize to upper-case.
- **Over-trimming.** Do not use the trimmed alignment for the main results.
- **Gap characters.** Standardize to `-`. MUSCLE and Clustal both use `-`, but check for `.`.
- **Mistaking an alignment artefact for biology.** A gap placed inside the CxxCH motif is almost certainly an alignment error. Investigate.

## 7. Done when

- [ ] Section 5 fully ticked.
- [ ] `make phase3` reproduces everything.
- [ ] `PROGRESS.md` updated; Phase 3 committed.
