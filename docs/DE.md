# Differential expression format

The differential expression format by quantms is based on the [MSstats](https://msstats.org/wp-content/uploads/2017/01/MSstats_v3.7.3_manual.pdf) output. The MSstats format is a tab-delimited file that contains the following fields - see example [file](include/PXD004683.csv):

- `Protein` -> Protein Accession
- `Label` -> Label for the Conditions combination.	
- `log2FC` -> Log2 Fold Change	
- `SE` -> Standard error of the log2 fold change 	
- `DF` -> Degree of freedom of the Student test	
- `pvalue`	-> Raw p-values
- `adj.pvalue`	->  P-values adjusted among all the proteins in the specific comparison using the approach by Benjamini and Hochberg
- `issue` -> Issue column shows if there is any issue for inference in corresponding protein and comparison,  for example, OneConditionMissing or CompleteMissing. 

Example: 

| Protein    | Label                          | log2FC | SE | DF | pvalue | adj.pvalue | issue |
| ---------  |--------------------------------| ------ | -- | -- | ------ | ---------- |-------|
|LV861_HUMAN | normal-squamous cell carcinoma | 0.60   | 0.87 | 8  | 0.51   | 0.62       | NA  |

## MSstats Header 

By default, the MSstats format does not have any header of metadata. We suggest adding a header to the output for better understanding of the file. By default, MSstats allows comments in the file if the line starts with `#`. The quantms output will start with some key value pairs that describe the project, the workflow and also the columns in the file. For example: 

`#project_accession=PXD000000`

In addition, for each `Default` column of the matrix the following information should be added: 

```
#INFO=<ID=Protein, Number=1, Type=String, Description="Protein Accession">
#INFO=<ID=Label, Number=1, Type=String, Description="Label for the Conditions combination">
#INFO=<ID=log2FC, Number=1, Type=Float, Description="Log2 Fold Change">
#INFO=<ID=SE, Number=1, Type=Float, Description="Standard error of the log2 fold change">
#INFO=<ID=DF, Number=1, Type=Integer, Description="Degree of freedom of the Student test">
#INFO=<ID=pvalue, Number=1, Type=Float, Description="Raw p-values">
#INFO=<ID=adj.pvalue, Number=1, Type=Float, Description="P-values adjusted among all the proteins in the specific comparison using the approach by Benjamini and Hochberg">
#INFO=<ID=issue, Number=1, Type=String, Description="Issue column shows if there is any issue for inference in corresponding protein and comparison">
```

The `ID` is the column name in the matrix, the `Number` is the number of values in the column (separated by `;`), the `Type` is the type of the values in the column and the `Description` is a description of the column.

We suggest including the following properties in the header: 

- project_accession: The project accession in PRIDE Archive
- project_title: The project title in PRIDE Archive
- project_description: The project description in PRIDE Archive
- quanmts_version: The version of the quantms workflow used to generate the file
- factor_value: The factor values used in the analysis (e.g. `phenotype`)
- fdr_threshold: The FDR threshold used to filter the protein lists (e.g. `adj.pvalue < 0.05`)


A complete example of a quantms output file can be seen [here](include/PXD004683-quantms.csv).



