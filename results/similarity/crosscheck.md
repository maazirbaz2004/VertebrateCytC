# Phase 5 Cross-Validation of Pairwise Sequence Identities

## 1. MSA-Derived Identity vs. Independent Global Pairwise Alignment (`Bio.Align.PairwiseAligner`)
- Settings: `mode=global`, `substitution_matrix=BLOSUM62`, `open_gap_score=-10.0`, `extend_gap_score=-0.5`.
- Pairs compared: **136**
- Pearson correlation: **r = 0.999032**
- Spearman correlation: **rho = 0.999831**
- Mean absolute difference: **0.3179** percentage points
- Maximum absolute difference: **1.9231** percentage points

### Pairs differing by > 2.0 percentage points
_None (all 136 pairs agree within <= 2.0 percentage points)._

## 2. Third Cross-Check: Clustal Omega Distance Matrix
- **Clustal Omega `--percent-id` distmat (`all_clustalo.distmat`) vs. MSA identity matrix:** Pearson r = **0.999779**, max absolute difference = **0.8059** percentage points.
