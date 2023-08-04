# Peptide table format

## Use cases

The Peptide table aims to cover detail on peptide level including peptide intensity. The most of content are from peptide part of mzTab.

- Store peptide intensity to perform down-stream analysis and integration.
- Enable easy visualization and scanning on peptide level.

## Format

For large-scale datasets, peptide section would be very large. Therefore, Parquet format is adopted and its data section mainly consists of the following column:

- `sequence`: Peptide sequence -> `string`
- `protein_accessions`: A list protein's accessions -> `list[string] (e.g. [P02768, P02769])`
- `unique`: Indicates whether the peptide is unique for this protein in respect to the searched database -> `boolean (0/1)`
- `best_search_engine_score_{}`: The best search engine score for the given peptide -> `double`
- `reference_file_name`: The reference file name that contains the spectrum. -> `string` `e.g. Adult_Frontalcortex_bRP_Elite_85_f09.mzML`
- `modifications`: A list of modifications for a give peptide -> `[modification1, modification2, ...]`. A modification should be recorded as string similarly to mztab like:
  - `{position}({Probabilistic Score:0.9})|{position2}|..-{modification accession or name}` -> e.g `1(Probabilistic Score:0.9)|2|3-UNIMOD:35`
- `retention_time`: Retention time (seconds) -> `double`
- `charge`: Precursor charge -> `int`
- `exp_mass_to_charge`: The precursor’s experimental mass to charge (m/z) -> `double`
- `scan_number`: The scan number of the spectrum for the quantification. The scan number or index of the spectrum in the file -> `string` `e.g. 1`
- `abundance_{label}`: The peptide’s abundance in the given sample in LFQ experiments or in defined spectrum file name combined with channel in TMT experiments. `{}` can be omitted in LFQ `(e.g. abundance_TMT126)` -> `double`
- `peptidoform`: Peptidoform of the peptide `PEPTIDE[+80.0]FORM` -> `string`
- `is_decoy`: Indicates whether the peptide sequence is decoy -> `boolean (0/1)`

Optional fields:

- `gene_accessions`: A list of gene accessions -> `list[string] (e.g. [ENSG00000139618, ENSG00000139618])`
- `gene_names`: A list of gene names -> `list[string] (e.g. [APOA1, APOA1])`
- `consensus_support`: Global consensus support scores for multiple search engines -> `double`
- `posterior_error_probability`: Posterior Error Probability scores -> `double`
- `search_engine_score_{}`: The search engine scores are the only fields that are not defined in the schema. Each search engine score will be an optional column with the prefix `search_engine_score_` and the value of the column is the corresponding value of the `(e.g. MS-GF:RawScore -> search_engine_score_MS-GF:RawScore)` -> `double`