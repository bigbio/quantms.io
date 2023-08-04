# Protein table format

## Use cases

The Protein table is a [Parquet](https://github.com/apache/parquet-format) file that contains the details of the proteins identified and quantified .

- Store proteins identified and quantified from mzTab file, with the corresponding abundance and search engine scores.
- Enable easy visualization and scanning on protein level.

## Format

- `protein_accession`: A list protein's accessions -> `list[string] (e.g. [P02768, P02769])`
- `best_id_score_{}`: The best search engine score for the given protein -> `double`
- `abundance_{label}`: The protein's abundance as measured in the given samples in LFQ experiments or in defined spectrum file name combined with channel in TMT tag experiments. `{}` can be omitted in LFQ `(e.g. abundance_TMT126)` -> `double`
- `global_qvalue`: Global q-value from quantms -> `double`
- `is_decoy`: Indicates whether the protein is decoy -> `boolean (0/1)`

Optional fields:

- `gene_accessions`: A list of gene accessions -> `list[string] (e.g. [ENSG00000139618, ENSG00000139618])`
- `gene_names`: A list of gene names -> `list[string] (e.g. [APOA1, APOA1])`