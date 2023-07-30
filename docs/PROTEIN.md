# Protein table format

## Use cases

The Protein table is a [Parquet](https://github.com/apache/parquet-format) file that contains the details of the proteins identified and quantified .

- Store proteins identified and quantified from mzTab file, with the corresponding abundance and search engine scores.
- Enable easy visualization and scanning on protein level.

## Format

- `accession`: Protein accession or comma-separated list of accessions for indistinguishable groups -> `string`
- `description`: The protein’s name and or description line -> `string`
- `best_search_engine_score`: The best search engine score for the given protein across all replicates reported -> `double`
- `modifications`: Modification names using a comma delimited list of modifications found in the given protein. -> `list[string]`
- `protein_coverage`: A value between 0 and 1 defining the protein coverage. -> `float`
- `protein_abundance_{sample identifier}/{file_name_channel}`: The protein's abundance as measured in the given samples in LFQ experiments or in defined spectrum file name combined with channel in TMT tag experiments, which usually consistent with `study variable` column from mzTab -> `double`
- `optional columns`: Additional required columns can be added to the end of the protein table for specific use-cases. These column headers MUST start with the prefix `opt_`
  - `opt_global_nr_found_peptides`: The number of peptides in the given protein -> `int`
  - `opt_global_cv_PRIDE:0000303_decoy_hit`: Indicates whether the protein is decoy hit -> `boolean (0/1)`

Example:

- LFQ

| accession | ...... | best_search_engine_score | protein_abundance_PXD002854-Sample-Erythrocytes-1 | protein_abundance_PXD002854-Sample-Erythrocytes-2 |
| -- | -- | ------ | ---------- | ---------- |
| sp\|P11277\|SPTB1_HUMAN | ......   | 0.9 | 14837130 | 0 |

- TMT/ITRAQ

| accession | ...... | best_search_engine_score | protein_abundance_20150820_Haura-Pilot-TMT1-bRPLC01-1_TMT126 | protein_abundance_20150820_Haura-Pilot-TMT1\-bRPLC01-1_TMT127 |
| -- | -- | ------ | ---------- | ---------- |
| sp\|Q9P287\|BCCIP_HUMAN | ...... | 0.8 | 137571.703125  | NA |