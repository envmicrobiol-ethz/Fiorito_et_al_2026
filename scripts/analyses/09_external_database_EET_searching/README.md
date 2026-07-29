## External-database EET searching

## BLASTn screening of external peatland MAGs

- `01_BLASTn_external_MAGs/01_run_BLASTn_external_MAGs.sh`: builds a nucleotide BLAST database from an external MAG dataset and searches it with nucleotide sequences from the conserved EET gene cassettes using BLASTn-plus v2.14.1.

- `01_BLASTn_external_MAGs/02_build_contig_to_MAG_mapping.py`: maps every target contig identifier to its source external MAG.

- `01_BLASTn_external_MAGs/03_call_credible_EET_cassettes_from_BLASTn.py`: applies the final nucleotide criteria (E-value ≤1 × 10⁻¹⁰, identity ≥95%, query coverage ≥95% and zero gap openings) and calls a cassette match credible when an essential EET gene and at least one additional gene from the same query cassette occur on the same target contig.









## mOTUs classification

The MAG classification can be documented as a command rather than a separate script:

motus classify -i "$INPUT_GENOMES_LIST" -o "$OUTPUT_TABLE" -t "$THREADS"
