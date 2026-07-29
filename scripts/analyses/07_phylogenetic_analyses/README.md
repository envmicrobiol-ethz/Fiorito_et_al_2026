## Scripts

- `01_build_all_MAGs_phylogeny_gtotree.sh`: constructs a genome phylogeny for the 1,081 species-level dereplicated MAGs using GToTree v1.8.6 and the combined bacterial and archaeal marker set.

- `02_build_EET_MAGs_and_reference_phylogeny_gtotree.sh`: constructs a phylogeny containing the 120 EET-bearing MAGs and 60 electroactive reference genomes using GToTree v1.8.6.

- `03_build_MHC_or_porin_protein_phylogeny.sh`: aligns selected MHC or porin homologs with MAFFT L-INS-i, trims the alignment with trimAl and infers a maximum-likelihood phylogeny with IQ-TREE 3 using ModelFinder and 1,000 ultrafast bootstrap replicates.

- `04_build_McrA_phylogeny.py`: identifies McrA candidates in proteins predicted from 1,081 dereplicated MAGs and assembled contigs from 32 metagenomes, applies reference-calibrated HMM score and coverage thresholds, collapses exact amino-acid duplicates while retaining source provenance, adds the resulting sequences to a published McrA reference alignment and infers a maximum-likelihood phylogeny using MAFFT, trimAl and IQ-TREE 3.
