## Scripts


- `01_identify_putative_MHCs.py`: scans predicted proteins from MAGs, reference genomes or metagenomic assemblies for CXXCH, CXXCX₀–₇₀H and CX₁–₇₀CH motifs and retains proteins containing at least three motifs of the same type as candidate multiheme cytochromes. This script was used for screening peat metagenomes and electroactive reference microbes mentioned in the study.
- `02_predict_MHC_subcellular_localization_psortb.sh`: predicts the subcellular localization of candidate MHCs using PSORTb v3.0. The Gram-negative model was used for the main metagenomic dataset, while Gram-negative, Gram-positive and archaeal models were applied as appropriate to electroactive reference microorganisms.
- `03_cluster_MHC_proteins_mmseqs2.sh`: clusters candidate MHCs from the metagenomic dataset together with 73 reference MHCs from *Geobacter sulfurreducens* PCA and 23 reference MHCs from *Shewanella oneidensis* MR-1 using MMseqs2 v14-7e284 (`--min-seq-id 0.50`, `-c 0.9`, `--cluster-reassign`), and exports cluster representatives and protein membership.
- `04_score_MHC_EET_credibility.py`: assigns cluster-level credibility scores to candidate MHC/EET protein families using multiheme motif content, PSORTb localization, PFAM support, heme density and normalized metatranscriptomic expression. Clusters are classified as Tier 1 credible, Tier 2 putative or Tier 3 uncertain.
