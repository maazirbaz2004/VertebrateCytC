# Phase 4 Conservation Analysis Summary (RQ1 & RQ2)

## 1. Headline Conservation Statistics (RQ1)
- Sequence weights: **unweighted** (n = 15 vertebrates, n = 17 all; no pseudo-count smoothing applied).
- Henikoff sequence-weighting sensitivity check (vertebrates): **13** / 105 columns changed conservation class under Henikoff weights.
- JSD background: published BLOSUM62 amino-acid frequencies (Capra & Singh 2007, window = 0, multiplied by `1 - gap_fraction`).

| dataset | total_columns | gapfree_columns | n_invariant | pct_invariant | n_highly_conserved | n_ge_90_conserved | pct_ge_90_conserved | n_moderately_conserved | pct_moderately_conserved | n_variable | pct_variable | n_gap_rich | mean_identity_frac | median_identity_frac | mean_entropy_bits | mean_conservation_score | mean_jsd_score | n_invariant_stretches_ge2 | longest_invariant_stretch_len | longest_invariant_stretch_uniprot |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| vertebrates_only | 105 | 104 | 67 | 63.81 | 16 | 83 | 79.05 | 11 | 10.48 | 11 | 10.48 | 0 | 0.9299 | 1.0 | 0.2741 | 0.9366 | 0.802 | 16 | 11 | 71-81 |
| all_17_species | 110 | 104 | 52 | 47.27 | 18 | 70 | 63.64 | 24 | 21.82 | 11 | 10.0 | 5 | 0.8925 | 0.9412 | 0.4148 | 0.904 | 0.748 | 10 | 11 | 71-81 |

## 2. Deep Conservation (Vertebrates -> All 17 Species)
- **52 / 67 (77.6%)** of vertebrate-invariant positions remain strictly invariant when *D. melanogaster* and *S. cerevisiae* are included.
- Total strictly invariant positions across all 17 species: **52 / 110 (47.3%)**.
- All-17 invariant residues (UniProt-human / classical): G2 (cl.1), K6 (cl.5), G7 (cl.6), F11 (cl.10), C15 (cl.14), Q17 (cl.16), C18 (cl.17), H19 (cl.18), T20 (cl.19), V21 (cl.20), G25 (cl.24), H27 (cl.26), K28 (cl.27), G30 (cl.29), P31 (cl.30), N32 (cl.31), L33 (cl.32), G35 (cl.34), G38 (cl.37), R39 (cl.38), G42 (cl.41), Q43 (cl.42), A44 (cl.43), G46 (cl.45), Y49 (cl.48), T50 (cl.49), A52 (cl.51), N53 (cl.52), K56 (cl.55), W60 (cl.59), E62 (cl.61), Y68 (cl.67), L69 (cl.68), N71 (cl.70), P72 (cl.71), K73 (cl.72), K74 (cl.73), Y75 (cl.74), I76 (cl.75), P77 (cl.76), G78 (cl.77), T79 (cl.78), K80 (cl.79), M81 (cl.80), F83 (cl.82), G85 (cl.84), K88 (cl.87), R92 (cl.91), D94 (cl.93), I96 (cl.95), Y98 (cl.97), K100 (cl.99).

## 3. Detected Conserved and Variable Regions (Vertebrates, RQ2)
Conserved regions = maximal runs of >= 3 consecutive `invariant` or `highly_conserved` columns (strict, unmerged). Variable regions = maximal runs of >= 2 consecutive `variable` columns.

| region_type | aln_start | aln_end | human_uniprot_start | human_uniprot_end | classical_start | classical_end | length | mean_conservation | sequence_in_human |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| conserved | 6 | 8 | 1 | 3 | nan | 2 | 3 | 1.0 | MGD |
| conserved | 10 | 16 | 5 | 11 | 4.0 | 10 | 7 | 0.9766 | EKGKKIF |
| conserved | 22 | 28 | 17 | 23 | 16.0 | 22 | 7 | 0.9766 | QCHTVEK |
| conserved | 30 | 33 | 25 | 28 | 24.0 | 27 | 4 | 0.9796 | GKHK |
| conserved | 35 | 38 | 30 | 33 | 29.0 | 32 | 4 | 1.0 | GPNL |
| conserved | 40 | 49 | 35 | 44 | 34.0 | 43 | 10 | 0.9918 | GLFGRKTGQA |
| conserved | 53 | 55 | 48 | 50 | 47.0 | 49 | 3 | 0.9727 | SYT |
| conserved | 57 | 59 | 52 | 54 | 51.0 | 53 | 3 | 0.9727 | ANK |
| conserved | 61 | 63 | 56 | 58 | 55.0 | 57 | 3 | 1.0 | KGI |
| conserved | 72 | 88 | 67 | 83 | 66.0 | 82 | 17 | 0.9856 | EYLENPKKYIPGTKMIF |
| conserved | 90 | 94 | 85 | 89 | 84.0 | 88 | 5 | 0.9673 | GIKKK |
| conserved | 99 | 105 | 94 | 100 | 93.0 | 99 | 7 | 0.9766 | DLIAYLK |
| variable | 109 | 110 | 104 | 105 | 103.0 | 104 | 2 | 0.7104 | NE |

## 4. Top 10 Most Variable Vertebrate Positions
| aln_col | human_uniprot_pos | classical_pos | entropy_bits | identity_frac | residues |
| --- | --- | --- | --- | --- | --- |
| 95 | 90 | 89 | 2.04 | 0.533 | G:8, E:2, T:2, N:1, R:1, S:1 |
| 50 | 45 | 44 | 1.746 | 0.467 | P:7, E:4, A:3, D:1 |
| 106 | 101 | 100 | 1.688 | 0.6 | K:9, D:3, E:1, S:1, V:1 |
| 109 | 104 | 103 | 1.563 | 0.467 | N:7, S:6, A:1, K:1 |
| 56 | 51 | 50 | 1.426 | 0.667 | D:10, A:2, E:2, Q:1 |
| 39 | 34 | 33 | 1.371 | 0.6 | H:9, N:3, W:3 |
| 64 | 59 | 58 | 1.242 | 0.667 | T:10, I:3, V:2 |
| 98 | 93 | 92 | 1.237 | 0.733 | A:11, E:2, L:1, V:1 |
| 60 | 55 | 54 | 1.159 | 0.667 | N:10, S:4, E:1 |
| 66 | 61 | 60 | 1.039 | 0.8 | G:12, D:1, K:1, S:1 |

## 5. Diagnostic Lineage-Specific Substitutions (Multi-Species Groups)
| aln_col | human_uniprot_pos | classical_pos | group | diagnostic_residue | other_residues_seen |
| --- | --- | --- | --- | --- | --- |
| 17 | 12 | 11 | Primate | I | V |
| 18 | 13 | 12 | Primate | M | Q |
| 34 | 29 | 28 | Fish | V | T |
| 56 | 51 | 50 | Primate | A | D,E,Q |
| 56 | 51 | 50 | Reptile | E | A,D,Q |
| 89 | 84 | 83 | Primate | V | A |
| 95 | 90 | 89 | Primate | E | G,N,R,S,T |

## 6. Differences vs. *Homo sapiens* (Vertebrates Quick Check)
| species | aa_differences_from_human |
| --- | --- |
| Homo_sapiens | 0 |
| Macaca_mulatta | 1 |
| Mus_musculus | 9 |
| Rattus_norvegicus | 9 |
| Canis_lupus_familiaris | 11 |
| Bos_taurus | 10 |
| Equus_caballus | 12 |
| Monodelphis_domestica | 10 |
| Gallus_gallus | 13 |
| Taeniopygia_guttata | 12 |
| Anolis_carolinensis | 13 |
| Chrysemys_picta | 16 |
| Xenopus_tropicalis | 14 |
| Danio_rerio | 17 |
| Takifugu_rubripes | 30 |
