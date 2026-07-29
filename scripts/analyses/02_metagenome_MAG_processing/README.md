## Scripts

- `01_singlem_raw_reads.sh`: estimates taxonomic profiles and microbial read fractions from raw paired-end metagenomic reads using SingleM.
- `02_fastqc_raw_reads.sh`: performs FastQC quality assessment of raw paired-end metagenomic reads before trimming.
- `03_bbduk_trimming_and_fastqc.sh`: performs adapter removal and quality trimming using BBDuk, followed by FastQC assessment of the cleaned reads.
- `04_metaspades_assembly.sh`: assembles cleaned paired-end metagenomic reads for each sample using metaSPAdes with k-mer sizes 21, 33, 55, 77, 99 and 127.
- `05a_bbmap_assembly_stats.sh`: calculates per-sample metaSPAdes assembly statistics using BBMap `stats.sh`.
- `05b_seqkit_assembly_stats.sh`: calculates per-sample metaSPAdes assembly statistics using SeqKit.
- `06_filter_scaffolds_min1000.sh`: removes scaffolds shorter than 1,000 nucleotides from each metaSPAdes assembly using BBMap `reformat.sh`.
- `07_predict_genes_with_prodigal.sh`: predicts protein-coding genes from scaffolds ≥1,000 nucleotides using Prodigal v2.6.3 in metagenomic mode.
- `08_differential_coverage_mapping_design.tsv`: defines the cross-sample read-to-assembly combinations used to generate differential scaffold-coverage profiles for genome binning.
- `08_generate_differential_coverage_profiles.sh`: maps selected cleaned metagenomic read sets to each filtered assembly using BBMap and produces sorted, indexed BAM files and scaffold-coverage profiles for differential-coverage binning.
- `09_metabat2_binning.sh`: generates differential scaffold-depth profiles from the selected cross-sample BAM files and performs per-assembly genome binning using MetaBAT2 v2.17.
- `10_maxbin2_binning.sh`: extracts all sample-specific scaffold-coverage profiles from the MetaBAT2 depth tables and performs per-assembly genome binning using MaxBin2 v2.2.7.
- `11_concoct_binning.sh`: generates differential scaffold-coverage tables from the selected cross-sample BAM files and performs per-assembly genome binning using CONCOCT v1.1.0 with 10-kb subcontigs.
- `12_prepare_bins_for_das_tool.sh`: standardizes MetaBAT2, MaxBin2 and CONCOCT bin names and generates the corresponding contig-to-bin tables required by DAS Tool.
- `13_das_tool_bin_selection.sh`: integrates the MetaBAT2, MaxBin2 and CONCOCT results and selects non-redundant genome bins for each assembly using DAS Tool v1.1.7.
- `14_drep_species_level_dereplication.sh`: collects DAS Tool-selected bins, assigns unique sample-prefixed filenames and performs species-level dereplication using dRep v3.4.5 with 95% primary and secondary ANI thresholds, ≥50% completeness and ≤10% contamination.
- `15_checkm2_quality_assessment.sh`: estimates completeness and contamination of the 1,081 species-level dereplicated MAGs using CheckM2 v1.0.2.
