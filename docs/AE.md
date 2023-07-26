# Absolute expression format

The absolute expression format by quantms is a tab-delimited file format that contains the following fields:

- `Protein` -> Protein Accession
- `Sample_ID` -> Label for the Conditions combination.
- `Condition` -> Condition name
- `iBAQ` -> iBAQ value
- `riBAQ` -> Relative iBAQ value

Example: 

| Protein    | Sample_ID     | Condition | iBAQ   | riBAQ  |
| ---------  |---------------|-----------|--------| -------|
|LV861_HUMAN | Sample-1      | heart     | 1234.1 | 12.34  |

## AE Header 

By default, the MSstats format does not have any header of metadata. We suggest adding a header to the output for better understanding of the file. By default, MSstats allows comments in the file if the line starts with `#`. The quantms output will start with some key value pairs that describe the project, the workflow and also the columns in the file. For example: 

`#project_accession=PXD000000`

In addition, for each `Default` column of the matrix the following information should be added: 

```
#INFO=<ID=Protein, Number=1, Type=String, Description="Protein Accession">
#INFO=<ID=Sample_ID, Number=1, Type=String, Description="Sample Accession in the SDRF">
#INFO=<ID=Condition, Number=1, Type=String, Description="Value of the factor value">
#INFO=<ID=iBAQ, Number=1, Type=Float, Description="Intensity based absolute quantification">
#INFO=<ID=riBAQ, Number=1, Type=Float, Description="relative iBAQ">
```

The `ID` is the column name in the matrix, the `Number` is the number of values in the column (separated by `;`), the `Type` is the type of the values in the column and the `Description` is a description of the column.

We suggest including the following properties in the header: 

- project_accession: The project accession in PRIDE Archive
- project_title: The project title in PRIDE Archive
- project_description: The project description in PRIDE Archive
- quanmts_version: The version of the quantms workflow used to generate the file
- factor_value: The factor values used in the analysis (e.g. `tissue`)


Please check also the differential expression example for more information [DE](DE.md)



