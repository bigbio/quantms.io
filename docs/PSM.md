# PSM table format

## Use cases

- The PSM table aims to cover detail on PSM level for AI/ML training and other use-cases. 
- Most of the content is similar to mzTab, a PSM would be a peptide identification in an specific msrun file.
- The representation should be in a parquet file. *
- Store details on PSM level including spectrum mz/intensity for specific use-cases such as AI/ML training.
- Fast and easy visualization and scanning on PSM level.


_* [Parquet](https://github.com/apache/parquet-format)_ file is a columnar storage format that supports nested data. For these large-scale analyses, Parquet has helped its users reduce storage requirements by at least one-third on large datasets, in addition, it greatly improved scan and deserialization time (web use-cases), hence the overall costs. The following table compares the savings as well as the speedup obtained by converting data into Parquet from CSV.

| Dataset    | Size on Amazon S3 | Query Run Time | Data Scanned |
| ---------  |--------------|-----------|--------|
|Data stored as CSV files | 1 TB     | 236 seconds     | 1.15 TB |
|Data stored in Apache Parquet Format | 130 GB     | 6.78 seconds     | 2.51 GB |

## Format

The [Avro PSM schema](psm.avsc) is used to define the PSM table format. The following table describes the fields of the PSM table.

- `sequence`: The peptide's sequence corresponding to the PSM -> `string`
- `psm_id`: A unique identifier for a PSM within the file -> `integer`
- `protein_accessions`: A list protein's accessions -> `list[string] (e.g. [P02768, P02769]`)` 
- `protein_start_positions`: A list of protein's start positions -> `list[int] (e.g. [1, 2])`
- `protein_end_positions`: A list of protein's end positions -> `list[int] (e.g. [10, 20])`
- `unique`: Indicates whether the peptide sequence (coming from the PSM) is unique for this protein in respect to the searched database -> `boolean (0/1)`
- `modifications`: A list of modifications for a give peptide -> `[modification1, modification2, ...]`. A modification should be recorded as string similarly to mztab like:
  - `{position}({Probabilistic Score:0.9})|{position2}|..-{modification accession or name}` -> e.g `1(Probabilistic Score:0.9)|2|3-UNIMOD:35`
- `retention_time`: The retention time of the spectrum -> `double`
- `charge`: The charge assigned by the search engine/software -> `integer`
- `exp_mass_to_charge`: The PSM’s experimental mass to charge (m/z) -> `double`
- `calc_mass_to_charge`: The PSM’s calculated (theoretical) mass to charge (m/z) -> `double`
- `reference_file_name`: The reference file name that contains the spectrum. -> `string` `e.g. Adult_Frontalcortex_bRP_Elite_85_f09.mzML`
- `scan_number`: The scan number of the spectrum. The scan number or index of the spectrum in the file -> `string` `e.g. 1`
- `peptidoform`: Peptidoform of the PSM `PEPTIDE[+80.0]FORM` -> `string`
- `posterior_error_probability`: Posterior Error Probability score from quantms -> `double`
- `global_qvalue`: Global q-value from quantms -> `double`
- `is_decoy`: Indicates whether the peptide sequence (coming from the PSM) is decoy -> `boolean (0/1)`

Optional fields:

- `gene_accessions`: A list of gene accessions -> `list[string] (e.g. [ENSG00000139618, ENSG00000139618])`
- `gene_names`: A list of gene names -> `list[string] (e.g. [APOA1, APOA1])`
- `consensus_support`: Global consensus support scores for multiple search engines -> `double`
- `mz`: A list of mz values for the spectrum -> `list[double]`
- `intesity`: A list of intensity values for the spectrum ->  `list[double]`
- `num_peaks`: The number of peaks in the spectrum, this is the size of previous lists intensity and mz -> `integer`
- `search_engine_score_{}`: The search engine scores are the only fields that are not defined in the schema. Each search engine score will be an optional column with the prefix `search_engine_score_` and the value of the column is the corresponding value of the `(e.g. MS-GF:RawScore -> search_engine_score_MS-GF:RawScore)` -> `double`