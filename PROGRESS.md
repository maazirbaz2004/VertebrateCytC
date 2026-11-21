# PROGRESS

| Phase | Status | Date | Notes |
|-------|--------|------|-------|
| 0 Setup | [x] | 2026-11-18 | python 3.11.18, biopython 1.88, clustalo 1.2.4, muscle 5.3, iqtree 3.1.4, fasttree 2.2.0, trimal 1.5.rev1. WSL moved to D:. |
| 1 Data collection | [x] | 2026-11-18 | 17/17 species fetched from UniProt (13 Swiss-Prot reviewed, 4 TrEMBL); 45 candidates audited; lengths 104-105 aa (vertebrates), 108 aa (D. melanogaster), 109 aa (S. cerevisiae); pytest 6 passed. |
| 2 Sequence preparation | [x] | 2026-11-21 | 17/17 sequences curated (15 vertebrates, 2 outgroups) in fixed biological order with 0 ambiguous residues; 14/15 vertebrates 105 aa (`Danio_rerio` 104 aa, MW 11.456 kDa slightly below 11.5 kDa note), outgroups 108/109 aa; all vertebrates basic (pI 9.46-9.79) with 2 Cys and `CxxCH` at UniProt pos 15 (classical pos 14); 1 identical-sequence pair (`Mus_musculus` & `Rattus_norvegicus`); 3 QC figures + 5 QC tables + manifest SHA-256 written; pytest 11 passed. |
| 3 Alignment | [ ] | | |
| 4 Conservation | [ ] | | |
| 5 Similarity | [ ] | | |
| 6 Phylogenetics | [ ] | | |
| 7 Interpretation | [ ] | | |
| 8 Reporting | [ ] | | |

## Deviations (anything that differs from the phase files, with justification)
- Phase 1 (`config/species.csv`): Quoted `"teleost, second fish"` in the `Takifugu_rubripes` `notes` field so CSV retains 9 columns cleanly. Four vertebrates lacking Swiss-Prot entries (*Monodelphis domestica* `A0A5F8HEK2`, *Taeniopygia guttata* `B5FZX0`, *Anolis carolinensis* `A0A803TKE1`, *Chrysemys picta* `A0A8C3PCA0`) were retrieved from TrEMBL after passing all validation checks (105 aa, `CxxCH`, standard residues, highest somatic identity).

## Open issues

