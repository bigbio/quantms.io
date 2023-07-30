# PSM table format

## Use cases

The PSM table aims to cover detail on PSM level for AI/ML training and other use-cases. The most of content are similar to mzTab. [Parquet](https://github.com/apache/parquet-format) file is a columnar storage format that supports nested data. For these large-scale analysis, Parquet has helped its users reduce storage requirements by at least one-third on large datasets, in addition, it greatly improved scan and deserialization time (web use-cases), hence the overall costs. The following table compares the savings as well as the speedup obtained by converting data into Parquet from CSV.

- Store details on PSM level including spectrum mz/intensity for specific use-cases such as AI/ML training.
- Fast and easy visualization and scanning on PSM level.

| Dataset    | Size on Amazon S3 | Query Run Time | Data Scanned |
| ---------  |--------------|-----------|--------|
|Data stored as CSV files | 1 TB     | 236 seconds     | 1.15 TB |
|Data stored in Apache Parquet Format | 130 GB     | 6.78 seconds     | 2.51 GB |

## Format

- `sequence`: The peptide's sequence corresponding to the PSM -> `string`
- `PSM_ID`: A unique identifier for a PSM within the file -> `integer`
- `accession`: The protein's accession the corresponding peptide sequence -> `string`
- `unique`: Indicates whether the peptide sequence (coming from the PSM) is unique for this protein in respect to the searched database -> `boolean (0/1)`
- `database`: The protein database used for the search and the peptide sequence comes from -> `string`
- `search_engines`: A `,` delimited list of search engine(s) that identified the PSM -> `list[string]`
- `search_engine_score`: The search engine score for the given PSM -> `double`
- `modifications`: The peptide's (coming from the PSM) modifications. -> `string`
- `retention_time`: The retention time of the spectrum -> `double`
- `charge`: The charge assigned by the search engine/software -> `integer`
- `exp_mass_to_charge`: The PSM’s experimental mass to charge (m/z) -> `double`
- `calc_mass_to_charge`: The PSM’s calculated (theoretical) mass to charge (m/z) -> `double`
- `spectra_ref`: USI of the PSM `mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555` -> `string`
- `pre`: Amino acid preceding the peptide (coming from the PSM) in the protein sequence -> `string`
- `post`: Amino acid following the peptide (coming from the PSM) in the protein sequence -> `string`
- `start`: The start position of the peptide (coming from the PSM) within the protein, counting 1 as the N-terminus of the protein -> `string`
- `end`: The end position of the peptide (coming from the PSM) within the protein, counting 1 as the N-terminus of the protein -> `string`
- `optional columns`: Additional required columns can be added to the end of the PSM table for specific use-cases. These column headers MUST start with the prefix `opt_`
  - `opt_global_Posterior_Error_Probability_score`: Posterior Error Probability score from quantms -> `double`
  - `opt_global_q-value`: Global q-value from quantms -> `double`
  - `opt_global_consensus_support`: Global consensus support scores for multiple search engines -> `float`
  - `opt_global_cv_MS:1002217_decoy_peptide`: Indicates whether the peptide sequence (coming from the PSM) is decoy -> `boolean (0/1)`
  - `opt_global_cv_MS:1000889_peptidoform_sequence`: Peptidoform sequence-> `string`
  - `opt_global_cv_MS:1000514_mz_array`: -> A data array of m/z values -> `array`
  - `opt_global_cv_MS:1000515_intensity_array`: A data array of intensity values -> `array`