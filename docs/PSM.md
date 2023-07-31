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
- `accessions`: A list protein's accessions with the start and end position for peptide: `[{accession:string, start:int, end:int}, ...]`}] 
- `unique`: Indicates whether the peptide sequence (coming from the PSM) is unique for this protein in respect to the searched database -> `boolean (0/1)`
- `search_engine_scores`: A list search engine score for the given PSM -> `[{name:string, score:double}, ...]`
- `modifications`: A list of modifications for a give peptide -> `[{name:string, position:int, localization_score:double}, ...]`
- `retention_time`: The retention time of the spectrum -> `double`
- `charge`: The charge assigned by the search engine/software -> `integer`
- `exp_mass_to_charge`: The PSM’s experimental mass to charge (m/z) -> `double`
- `calc_mass_to_charge`: The PSM’s calculated (theoretical) mass to charge (m/z) -> `double`
- `spectra_usi`: USI of the PSM `mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555` -> `string`
- `peptifoform`: Peptiform of the PSM `PEPTIDE[+80.0]FORM` -> `string`
- `posterior_error_probability`: Posterior Error Probability score from quantms -> `double`
- `global_qvalue`: Global q-value from quantms -> `double`
- `consensus_support`: Global consensus support scores for multiple search engines -> `double`
- `isdecoy`: Indicates whether the peptide sequence (coming from the PSM) is decoy -> `boolean (0/1)`
- `spectrum`: A record with two array list mz and intensity ->  `{[double], [double]}`