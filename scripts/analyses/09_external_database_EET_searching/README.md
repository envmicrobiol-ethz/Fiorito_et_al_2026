## External-database EET searching

## BLASTn screening of external peatland MAGs (SPRUCE and Stordalen MIRE) - nucleotide level

- `01_BLASTn_external_MAGs/01_run_BLASTn_external_MAGs.sh`: builds a nucleotide BLAST database from an external MAG dataset and searches it with nucleotide sequences from the conserved EET gene cassettes using BLASTn-plus v2.14.1.

- `01_BLASTn_external_MAGs/02_build_contig_to_MAG_mapping.py`: maps every target contig identifier to its source external MAG.

- `01_BLASTn_external_MAGs/03_call_credible_EET_cassettes_from_BLASTn.py`: applies the final nucleotide criteria (E-value ≤1 × 10⁻¹⁰, identity ≥95%, query coverage ≥95% and zero gap openings) and calls a cassette match credible when an essential EET gene and at least one additional gene from the same query cassette occur on the same target contig.


## tBLASTn screening of external peatland MAGs (SPRUCE and Stordalen MIRE) - protein level

- `02_TBLASTN_external_MAGs/01_translate_EET_CDS_to_proteins.py`: translates nucleotide CDS sequences from the conserved EET cassettes using bacterial genetic code 11 while retaining the original query identifiers.

- `02_TBLASTN_external_MAGs/02_run_TBLASTN_external_MAGs.sh`: searches the translated EET query proteins against an external MAG nucleotide dataset using TBLASTN from BLAST+ v2.14.1.

- `02_TBLASTN_external_MAGs/03_build_explicit_EET_query_cassette_map.py`: reconstructs explicit EET query cassettes by splitting query genes on each MAG contig into consecutive runs of Prodigal gene numbers.

- `02_TBLASTN_external_MAGs/04_call_credible_EET_cassettes_from_TBLASTN.py`: applies the final protein-level criteria (E-value ≤1 × 10⁻¹⁰, amino-acid identity ≥50% and query coverage ≥90%), clusters overlapping HSPs into target loci and calls a cassette credible when an essential EET gene and a different gene from the same query cassette map to distinct loci on the same target contig.

The contig-to-MAG mapping script is shared with the nucleotide workflow:

01_BLASTn_external_MAGs/02_build_contig_to_MAG_mapping.py




## mOTUs classification

The MAG classification can be documented as a command rather than a separate script:

motus classify -i "$INPUT_GENOMES_LIST" -o "$OUTPUT_TABLE" -t "$THREADS"
