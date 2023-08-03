# Peptide evidence table format

## Use cases

The Peptide evidence table aims to cover detail on peptide level, including peptide intensity in relation to the sample metadata. The peptide evidence is the combination of between the MSstats, mzTab and Triqler peptide tables:

- Store peptide intensities in relation to the sample metadata to perform down-stream analysis and integration. This file can be used as input of MSstats and ibaqpy for protein quantification. 
- Enable peptide level statistics and algorithms to move from peptide level to protein level.

*NOTE*: quantms also release the peptide table for MSstats. The objective of the PEPTIDE_EVIDENCE table is to provide a more general peptide table and improve the annotations of the peptides with more columns. 

## Format

For large-scale datasets, peptide section would be very large. Therefore, Parquet format is adopted and its data section mainly consists of the following column:

- `sequence`: Peptide sequence -> `string`
- `accession`: Protein accession -> `string`
- `unique`: Indicates whether the peptide is unique for this protein in respect to the searched database -> `boolean (0/1)`
- `best_search_engine_score`: The best search engine score for the given peptide -> `double`
- `search_engine_score_{file_name}`: The search engine score for the given peptide in the defined ms run -> `double`
- `modifications`: The peptide's modifications or substitutions -> `string`
- `retention_time`: Retention time (seconds) -> `double`
- `charge`: Precursor charge -> `int`
- `mass_to_charge`: The precursor’s experimental mass to charge (m/z) -> `double`
- `spectra_ref`: USI of the best PSM -> `string`
- `peptide_abundance_{sample identifier}/{file_name_channel}`: The peptide’s abundance in the given sample in LFQ experiments or in defined spectrum file name combined with channel in TMT experiments, which usually consistent with `study variable` column from mzTab -> `double`
- `optional columns`: Additional required columns can be added to the end of the peptide table. These column headers MUST start with the prefix `opt_`
  - `opt_global_consensus_support`: Global consensus support scores for multiple search engines -> `float`
  - `opt_global_Posterior_Error_Probability_score`: Posterior Error Probability scores -> `double`
  - `opt_global_cv_MS:1002217_decoy_peptide`: Indicates whether the peptide is decoy -> `boolean (0/1)`
  - `opt_global_cv_MS:1000889_peptidoform_sequence`: Peptidoform sequence -> `string`

Example:

- LFQ

| sequence | accession | ...... | spectra_ref | peptide_abundance_PXD002854-Sample-Erythrocytes-1 | peptide_abundance_PXD002854-Sample-Erythrocytes-2 |
| --------------- | -- | -- | ------ | ---------- | ---------- |
|VHNYCVDCEETSK | sp\|P11277\|SPTB1_HUMAN | ......   | mzspec:PXD002854:20141201_qEp2_LC11_PhGe_SA_4_49_10_1:scan:3541 | 14837130 | 0 |

- TMT/ITRAQ

| sequence | accession | ...... | spectra_ref | peptide_abundance_20150820_Haura-Pilot-TMT1-bRPLC01-1_TMT126 | peptide_abundance_20150820_Haura-Pilot-TMT1\-bRPLC01-1_TMT127 |
| --------------- | -- | -- | ------ | ---------- | ---------- |
| TFVEAGKNNSK | sp\|Q9P287\|BCCIP_HUMAN | ...... | mzspec:PXD004683:20150820_Haura-Pilot-TMT1-bRPLC01-1:scan:31 | 137571.703125  | NA |