## Scripts


- `01_identify_putative_MHCs.py`: scans predicted proteins from MAGs, reference genomes or metagenomic assemblies for CXXCH, CXXCX₀–₇₀H and CX₁–₇₀CH motifs and retains proteins containing at least three motifs of the same type as candidate multiheme cytochromes. This script was used for screening peat metagenomes and electroactive reference microbes mentioned in the study.
- `02_predict_MHC_subcellular_localization_psortb.sh`: predicts the subcellular localization of candidate MHCs using PSORTb v3.0. The Gram-negative model was used for the main metagenomic dataset, while Gram-negative, Gram-positive and archaeal models were applied as appropriate to electroactive reference microorganisms.
