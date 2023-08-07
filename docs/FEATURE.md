# Peptide Features table format

## Use cases

The features table (peptide features) aims to cover detail on peptide level, including peptide intensity in relation to the sample metadata. The `feature parquet file` is the combination of between the MSstats, mzTab and Triqler peptide tables:

- Store peptide intensities in relation to the sample metadata to perform down-stream analysis and integration. This file can be used as input of MSstats and ibaqpy for protein quantification. 
- Enable peptide level statistics and algorithms to move from peptide level to protein level.

*NOTE*: quantms also release the peptide table for MSstats. The objective of the feature table is to provide a more general peptide table and improve the annotations of the peptides with more columns. 

## Format

Peptide properties and columns: 

- `sequence`: The peptide's sequence corresponding to the PSM -> `string`
- `protein_accessions`: A list protein's accessions -> `list[string]` 
- `protein_start_positions`: A list of protein's start positions -> `list[int]`
- `protein_end_positions`: A list of protein's end positions -> `list[int]`
- `unique`: Indicates whether the peptide sequence (coming from the PSM) is unique for this protein in respect to the searched database -> `boolean (0/1)`
- `modifications`: A list of modifications for a give peptide `[modification1, modification2, ...]`. A modification should be recorded as string like [modification definition](README.md#modifications)-> `list[string]`
- `retention_time`: The retention time of the spectrum -> `float`
- `charge`: The charge assigned by the search engine/software -> `integer`
- `exp_mass_to_charge`: The PSM’s experimental mass to charge (m/z) -> `double`
- `calc_mass_to_charge`: The PSM’s calculated (theoretical) mass to charge (m/z) -> `double`
- `peptidoform`: Peptidoform of the PSM. See more [documentation here](README.md#peptidoform). -> `string`
- `posterior_error_probability`: Posterior Error Probability score from quantms -> `double`
- `global_qvalue`: Global q-value from quantms -> `double`
- `is_decoy`: Indicates whether the peptide sequence (coming from the PSM) is decoy -> `boolean (0/1)`
- `best_id_score_{}`: The best search engine score is the only field that is not defined in the schema. Each search engine score will be an optional column with the prefix `best_id_score_` and the value of the column is the corresponding value of the `(e.g. MS-GF:RawScore -> best_id_score_MS-GF:RawScore)` -> `double`
- `intensity`: The abundance of the peptide in the sample -> `double`
- `spectral_count`: The number of spectra that match the peptide -> `integer`

Properties and columns from sample: 

- `sample_accession`: The sample accession in the sdrf which column is called `source name` -> `string`
- `condition`: The value for the factor value column in the sdrf, for example, the tissue name for the given sample in the column `factor value[organism part]` -> `string`
- `fraction`: The index value in the SDRF for the fraction column -> `string`
- `biological_replicate`: The value of the biological replicate column in the SDRF in relation with the condition -> `string`
- `fragment_ion`: The column defines a spectral feature: fragment ions `e.g. y7`. If information for the column is not available or not applicable, it should be set to a constant value `NA` -> `string`
- `isotope_label_type`: The column indicates whether the measurement is based on an endogenous peptide (indicated by value `L` or `light`) or reference peptide (indicated by value `H` or `heavy`) -> `string`
- `run`: The column stores IDs of mass spectrometry runs for LFQ experiments `e.g. 1`. For TMT/iTRAQ experiments, it is a identifier of mixture combined with technical replicate and fractions `{mixture}_{technical_replicate}_{fraction}` `e.g. 1_2_3` -> `string`
- `channel`: The channel used to label the sample (e.g. TMT115)-> `string`

Optional fields:

- `gene_accessions`: A list of gene accessions -> `list[string]`
- `gene_names`: A list of gene names -> `list[string]`
- `consensus_support`: Global consensus support scores for multiple search engines -> `double`
- `id_score_{}`: The search engine scores are the only fields that are not defined in the schema. Each search engine score will be an optional column with the prefix `id_score_` and the value of the column is the corresponding value of the `(e.g. MS-GF:RawScore -> id_score_MS-GF:RawScore)` -> `double`
- `reference_file_name`: The reference file name that contains the spectrum. -> `string` 
- `scan_number`: The scan number of the spectrum. The scan number or index of the spectrum in the file -> `string` 
