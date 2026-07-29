## Scripts

- `01_process_16S_amplicon_reads.sh`: processes raw paired-end 16S rRNA amplicon reads, including PhiX removal, primer trimming, read merging, quality filtering, ASV inference, chimera removal, target verification, OTU clustering, count-table generation and taxonomic classification against SILVA v138.
- `02_alpha_diversity_analysis.R`: filters the OTU table, performs 100 rarefactions to the minimum sequencing depth, calculates mean observed richness and Shannon diversity, derives Pielou's evenness, applies bestNormalize transformations and performs two-way ANOVA.
- `03_beta_diversity_analysis.R`: calculates Bray–Curtis and binary Jaccard dissimilarities between consecutive-depth samples across 100 rarefied OTU tables and summarizes the resulting values by peatland, sampling profile and depth.
- `04_prepare_16S_community_table.R`: combines the 16S rRNA OTU count and taxonomy tables, removes non-bacterial, non-archaeal, chloroplast and mitochondrial OTUs, and calculates sample-level relative abundances for the taxonomic composition figures.
