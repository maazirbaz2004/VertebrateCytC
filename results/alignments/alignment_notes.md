# Phase 3 Alignment Notes and Primary Alignment Decision

## 1. Tool Versions
- Clustal Omega: `1.2.4`
- MUSCLE: `muscle 5.3.linux64 []`
- trimAl: `trimAl v1.5.rev1 build[2025-11-25]`

## 2. Exact Commands Executed
```bash
clustalo -i /mnt/c/Users/Mr Laptop point/OneDrive/Desktop/BioInfoProject/data/processed/cytc_curated_treeid.fasta -o /mnt/c/Users/Mr Laptop point/OneDrive/Desktop/BioInfoProject/results/alignments/all_clustalo.fasta --seqtype Protein --outfmt fasta --force --threads 4 -v --percent-id --distmat-out=/mnt/c/Users/Mr Laptop point/OneDrive/Desktop/BioInfoProject/results/alignments/all_clustalo.distmat --full
muscle -align /mnt/c/Users/Mr Laptop point/OneDrive/Desktop/BioInfoProject/data/processed/cytc_curated_treeid.fasta -output /mnt/c/Users/Mr Laptop point/OneDrive/Desktop/BioInfoProject/results/alignments/all_muscle.fasta -threads 4
clustalo -i /mnt/c/Users/Mr Laptop point/OneDrive/Desktop/BioInfoProject/data/processed/cytc_vertebrates.fasta -o /mnt/c/Users/Mr Laptop point/OneDrive/Desktop/BioInfoProject/results/alignments/vert_clustalo.fasta --seqtype Protein --outfmt fasta --force --threads 4 -v --percent-id --distmat-out=/mnt/c/Users/Mr Laptop point/OneDrive/Desktop/BioInfoProject/results/alignments/vert_clustalo.distmat --full
muscle -align /mnt/c/Users/Mr Laptop point/OneDrive/Desktop/BioInfoProject/data/processed/cytc_vertebrates.fasta -output /mnt/c/Users/Mr Laptop point/OneDrive/Desktop/BioInfoProject/results/alignments/vert_muscle.fasta -threads 4
trimal -in /mnt/c/Users/Mr Laptop point/OneDrive/Desktop/BioInfoProject/results/alignments/primary_all.fasta -out /mnt/c/Users/Mr Laptop point/OneDrive/Desktop/BioInfoProject/results/alignments/primary_all_trimmed.fasta -gt 0.8 -colnumbering > /mnt/c/Users/Mr Laptop point/OneDrive/Desktop/BioInfoProject/results/alignments/trimmed_columns.txt
```

## 3. Alignment Statistics Summary
Sum-of-pairs (SP) score uses BLOSUM62 with pairwise affine gap penalties (gap_open=-10.0, gap_extend=-0.5; gap-gap columns ignored for each pair).

| alignment | dataset | aligner | n_seq | aln_length | n_gap_columns | n_gapfree_columns | n_invariant_columns | pct_identity_columns | mean_pairwise_identity | total_gap_chars | n_gap_blocks | sum_of_pairs_blosum62 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all_clustalo | all | clustalo | 17 | 110 | 6 | 104 | 52 | 47.27 | 84.65 | 79 | 19 | 65030.5 |
| all_muscle | all | muscle | 17 | 110 | 6 | 104 | 52 | 47.27 | 84.71 | 79 | 19 | 65088.5 |
| vert_clustalo | vert | clustalo | 15 | 105 | 1 | 104 | 67 | 63.81 | 88.73 | 1 | 1 | 52941.0 |
| vert_muscle | vert | muscle | 15 | 105 | 1 | 104 | 67 | 63.81 | 88.68 | 1 | 1 | 52913.0 |
| primary_all_trimmed | all_trimmed | clustalo+trimal | 17 | 105 | 1 | 104 | 52 | 49.52 | 84.67 | 3 | 3 | 65398.0 |

## 4. Aligner Agreement (Clustal Omega vs. MUSCLE v5)
Column-level agreement (TC score) counts only informative columns containing at least 2 non-gap residues (all-gap and single-residue columns are excluded as uninformative).

| dataset | primary_aligner | secondary_aligner | n_informative_primary | n_informative_secondary | n_shared_columns | pct_shared_columns_vs_primary | pct_shared_columns_jaccard | n_pairs_primary | n_shared_pairs | pct_shared_pairs_vs_primary | pct_shared_pairs_jaccard | disagreed_primary_cols |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all | clustalo | muscle | 109 | 109 | 105 | 96.33 | 92.92 | 14239 | 14195 | 99.69 | 99.38 | 2-3, 109-110 |
| vert | clustalo | muscle | 105 | 105 | 103 | 98.1 | 96.26 | 11011 | 10997 | 99.87 | 99.75 | 104-105 |

## 5. Biological Anchor Verification
- **Primary all (`primary_all.fasta`)**: Heme motif `CxxCH` aligned gap-free across all 17 species at alignment columns **20–24** (UniProt-human 15–19; classical 14–18). Distal axial ligand **Met80** at alignment column **86** (UniProt-human 81; classical 80).
- **Primary vertebrates (`primary_vert.fasta`)**: Heme motif `CxxCH` aligned gap-free across all 15 vertebrates at alignment columns **15–19**. Met80 loop (classical 70–85; UniProt-human 71–86) spans alignment columns **71–86** with zero gaps.

## 6. Primary Alignment Decision
- **Selected primary aligner:** `clustalo` (per `config.yaml` -> `alignment.primary`).
- Both Clustal Omega and MUSCLE v5 produce identical or near-identical vertebrate alignments, place the heme motif `CxxCH` and Met80 loop in gap-free homologous columns, and confine gaps primarily to the outgroup N-terminal extensions (*S. cerevisiae* and *D. melanogaster*) and the 1-aa C-terminal truncation in *Danio rerio*.

## 7. Gap Blocks in Primary Alignments
| alignment | start | end | length | sequences_with_gap | sequences_without_gap |
| --- | --- | --- | --- | --- | --- |
| all_clustalo | 1 | 1 | 1 | Homo_sapiens;Macaca_mulatta;Mus_musculus;Rattus_norvegicus;Canis_lupus_familiaris;Bos_taurus;Equus_caballus;Monodelphis_domestica;Gallus_gallus;Taeniopygia_guttata;Anolis_carolinensis;Chrysemys_picta;Xenopus_tropicalis;Danio_rerio;Takifugu_rubripes;Drosophila_melanogaster | Saccharomyces_cerevisiae |
| all_clustalo | 2 | 5 | 4 | Homo_sapiens;Macaca_mulatta;Mus_musculus;Rattus_norvegicus;Canis_lupus_familiaris;Bos_taurus;Equus_caballus;Monodelphis_domestica;Gallus_gallus;Taeniopygia_guttata;Anolis_carolinensis;Chrysemys_picta;Xenopus_tropicalis;Danio_rerio;Takifugu_rubripes | Drosophila_melanogaster;Saccharomyces_cerevisiae |
| all_clustalo | 110 | 110 | 1 | Danio_rerio;Drosophila_melanogaster;Saccharomyces_cerevisiae | Homo_sapiens;Macaca_mulatta;Mus_musculus;Rattus_norvegicus;Canis_lupus_familiaris;Bos_taurus;Equus_caballus;Monodelphis_domestica;Gallus_gallus;Taeniopygia_guttata;Anolis_carolinensis;Chrysemys_picta;Xenopus_tropicalis;Takifugu_rubripes |
| vert_clustalo | 105 | 105 | 1 | Danio_rerio | Homo_sapiens;Macaca_mulatta;Mus_musculus;Rattus_norvegicus;Canis_lupus_familiaris;Bos_taurus;Equus_caballus;Monodelphis_domestica;Gallus_gallus;Taeniopygia_guttata;Anolis_carolinensis;Chrysemys_picta;Xenopus_tropicalis;Takifugu_rubripes |

## 8. Trimmed Alignment (Sensitivity Analysis Only)
- Generated with trimAl (`-gt 0.8 -colnumbering`) as `primary_all_trimmed.fasta`.
- Used **only** in Phase 6 phylogenetic sensitivity testing, never as the primary alignment.
