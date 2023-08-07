# Peptide table format

## Use cases

The Peptide table aims to cover detail on peptide level including peptide intensity. The most of content are from peptide part of mzTab.

- Store peptide intensity to perform down-stream analysis and integration.
- Enable easy visualization and scanning on peptide level.

## Format

For large-scale datasets, a peptide section would be very large. Therefore, a Parquet format is adopted and its data section mainly consists of the following column:

- `sequence`: Peptide sequence -> `string`
- `protein_accessions`: A list protein's accessions -> `list[string] (e.g. [P02768, P02769])`
- `unique`: Indicates whether the peptide is unique for this protein in respect to the searched database -> `boolean (0/1)`
- `best_id_score`: A key value pair of the best search engine score selected by the algorithm `(e.g. "MS-GF:RawScore": 234.0)` -> `string`
- `posterior_error_probability`: Posterior Error Probability scores -> `double`
- `modifications`: A list of modifications for a give peptide -> `[modification1, modification2, ...]`. A modification should be recorded as string similarly to mztab like:
  - `{position}({Probabilistic Score:0.9})|{position2}|..-{modification accession or name}` -> e.g `1(Probabilistic Score:0.9)|2|3-UNIMOD:35`
- `charge`: Precursor charge -> `int`
- `exp_mass_to_charge`: The precursor’s experimental mass to charge (m/z) -> `double`
- `peptidoform`: Peptidoform of the peptide `PEPTIDE[+80.0]FORM` -> `string`
- `sample_accession`: A unique sample accession corresponding to the source name in the SDRF-> `string`
- `abundance`: The peptide’s abundance in the given sample -> `float`
- `is_decoy`: Indicates whether the peptide sequence is decoy -> `boolean (0/1)`

Optional fields:

- `number_of_psms`: Number of PSMs for the peptide in the given sample `sample_accession` -> `int`
- `retention_time`: Retention time (seconds), it can be the median across all retention times in the [Feature File](FEATURE.md) -> `float`
- `gene_accessions`: A list of gene accessions -> `list[string] (e.g. [ENSG00000139618, ENSG00000139618])`
- `gene_names`: A list of gene names -> `list[string] (e.g. [APOA1, APOA1])`
- `consensus_support`: Global consensus support scores for multiple search engines -> `float`
- `id_scores`: A list of identification scores, search engine, percolator etc. Each search engine score will be a key/value pair `(e.g. "MS-GF:RawScore": 78.9)` -> `list[string]`
- `reference_file_name`: The reference file name that contains the spectrum. -> `string` 
- `scan_number`: The scan number of the spectrum. The scan number or index of the spectrum in the file -> `string` 
