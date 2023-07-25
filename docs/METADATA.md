# Metadata Files

The `.qms` folder will contain multiple metadata files that will be used to describe the project, the samples, the data acquisition and the data processing. 

## Project file

The project file is a json file that contains the metadata of the project. The project file is used to link the different files of the project and to store the metadata of the project. 
The project file is a json file that contains the following fields:

- `project_accession` -> ProteomeXchange Identifier -> `string`
- `project_title` -> Project title -> `string` 
- `project_description` -> Project description -> `string`
- `project_sample_description` -> Sample description of the project -> `string` 
- `project_data_description` -> Data description of the project -> `string`
- `project_pubmed_id` -> PubMed identifier -> `string`
- `organism` -> Organism name -> `list[string]`
- `organism_part` -> Organism part -> `list[string]`
- `disease` -> Disease -> `list[string]`
- `cell line` -> Cell line (if available) -> `list[string]`
- `instrument` -> Instrument name -> `list[string]`
- `data_acquisition` -> Data acquisition type (DDA, DIA, TMT) -> `list[string]`
- `fragmentation_type` -> The fragmentation type, mode (HCD, ...) -> `list[string]`
- `enzyme` -> Protease type for digest -> `list[string]`
- `experiment type` -> Keywords in ProteomeXchange or PRIDE around the dataset. -> `list[string]`

## Sample file

We only provide here the SDRF format used to analyze the data with quantms. The SDRF file is a tab-delimited file that contains the metadata of the samples. 
The SDRF file is used to link the different files of the project and to store the metadata of the samples.

Read [here](https://github.com/bigbio/proteomics-sample-metadata/tree/master/sdrf-proteomics) more about SDRF. 
