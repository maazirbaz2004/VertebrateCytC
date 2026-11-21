# Phase 5 Pairwise Similarity Metric Definitions

**Headline metrics:** `identity_pct` (gap-free pairwise aligned columns) and `similarity_pct` (BLOSUM62 positive substitutions over gap-free pairwise aligned columns).

| Metric | Definition |
| --- | --- |
| `n_compared` | Number of alignment columns where **both** sequences contain a standard amino acid (gap-gap and residue-gap columns excluded). |
| `identity_pct` | `100 * identical / n_compared` (headline percent identity). |
| `identity_pct_overall` | `100 * identical / min(len_i, len_j)` using ungapped sequence lengths (penalizes terminal extensions/indels). |
| `similarity_pct` | `100 * (compared columns with BLOSUM62 score > 0) / n_compared` (includes identical residues and biochemically conservative substitutions). |
| `p_distance` | `1 - identity_pct / 100` (observed proportion of differing sites). |
| `poisson_distance` | `-ln(1 - p_distance)` (Poisson correction for multiple substitutions; `0` when `p = 0`; guarded if `p >= 1`). |
| `n_substitutions` | Number of non-identical compared columns (`n_compared - identical`). |
| `blosum62_score` | Sum of `BLOSUM62` substitution scores across all `n_compared` columns. |
