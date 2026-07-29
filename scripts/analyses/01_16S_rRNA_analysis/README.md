## Scripts

- `01_process_16S_amplicon_reads.sh`: processes raw paired-end 16S rRNA amplicon reads, including PhiX removal, primer trimming, read merging, quality filtering, ASV inference, chimera removal, target verification, OTU clustering, count-table generation and taxonomic classification against SILVA v138.
- `02_alpha_diversity_analysis.R`: filters the OTU table, performs 100 rarefactions to the minimum sequencing depth, calculates mean observed richness and Shannon diversity, derives Pielou's evenness, applies bestNormalize transformations and performs two-way ANOVA.
