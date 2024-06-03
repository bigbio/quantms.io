# quantms.io
[![Python application](https://github.com/bigbio/quantms.io/actions/workflows/python-app.yml/badge.svg?branch=dev)](https://github.com/bigbio/quantms.io/actions/workflows/python-app.yml)
[![Upload Python Package](https://github.com/bigbio/quantms.io/actions/workflows/python-publish.yml/badge.svg)](https://github.com/bigbio/quantms.io/actions/workflows/python-publish.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e71a662e8d4f483094576c1d8f8888c3)](https://app.codacy.com/gh/bigbio/quantms.io/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/e71a662e8d4f483094576c1d8f8888c3)](https://app.codacy.com/gh/bigbio/quantms.io/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_Coverage)
[![Documentation Status](https://readthedocs.org/projects/quantmsio/badge/?version=latest)](https://quantmsio.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/quantmsio.svg)](https://badge.fury.io/py/quantmsio)

[quantms](https://quantms.org) is a nextflow pipeline for the analysis of quantitative proteomics data. The pipeline is based on the [OpenMS](https://www.openms.de/) framework and [DIA-NN](https://github.com/vdemichev/DiaNN); and it is designed to analyze large scale experiments. The main outputs of quantms workflow are the following: 

- [mzTab](https://github.com/HUPO-PSI/mzTab) files with the identification and quantification information.
- [MSstats](https://msstats.org/wp-content/uploads/2017/01/MSstats_v3.7.3_manual.pdf) input file with the peptide quantification values needed for the MSstats analysis.
- [MSstats](https://msstats.org/wp-content/uploads/2017/01/MSstats_v3.7.3_manual.pdf) output file with the differential expression values for each protein. 
- The input [SDRF](https://github.com/bigbio/proteomics-sample-metadata) of the pipeline if available. 

While all the previous formats are well-known standards and popular formats in the proteomics community; they are difficult to use in big data analysis projects. In addition, these file formats are difficult to extend and provide multiple views of the underlying data. For example, in mzTab it is extremely hard for big datasets to retrieve the identified peptides and features and the corresponding intensities. At the same time it is difficult to get the protein quantification values for a given sample.  

Here, we aim to formalize and develop a more standardized format that enables better representation of the identification and quantification results but also enables new and novel use cases for proteomics data analysis. The main use cases for the format are:  

- Fast and easy visualization of the identification and quantification results.
- Easy integration with other omics data.
- Easy integration with sample metadata.
- AI/ML model development based on identification and quantification results.
- Easy data retrieval for big datasets and large-scale collections of proteomics data.

>**Note**: We are not trying to replace the mzTab format, but to provide a new format that enables AI-related use cases. Most of the features of the mzTab format will be included in the new format.  

## Data model

quantms.io could be seen as a **multiple view** representation of a proteomics data analysis results. Each view of the format can be serialized in different formats depending on the use case. the **data model** of quantms.io defines two main things, the **view** and how the view is **serialized**. 

- The **data model view** defines the structure, the fields and properties that will be included in a view for each peptide, psms, feature or protein, for example.    
- The **data serialization** defines the format in which the view will be serialized and what features of serialization will be supported, for example compression, indexing or slicing.

| view         | file class        | serialization format | definition                                                      | example                                                                                                       |
|:-------------|:------------------|:---------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------|
| psm          | psm_file          | _parquet_            | [psm](docs/psm.rst)                                             | [psm example](docs/include/PXD002854-80934754-57c1-47e2-9951-787ef703a484.psm.parquet)                        |
| feature      | feature_file      | _parquet_            | [feature](docs/feature.rst)                                     | [feature example](docs/include/PXD004683-219a8a0a-d6a8-44c9-9e51-1851876d2f69.feature.parquet)                |
| absolute     | absolute_file     | _tsv_                | [absolute](docs/ae.rst)                                         | [absolute example](docs/include/PXD004683-quantms.tsv)                                                        |
| differential | differential_file | _tsv_                | [differential](docs/differential.rst)                           | [differential example](docs/include/PXD004683-219a8a0a-d6a8-44c9-9e51-1851876d2f69.differential.tsv)          |
| sdrf         | sdrf_file         | _tsv_                | [metadata](https://github.com/bigbio/proteomics-sample-metadata)| [sdrf example](https://github.com/bigbio/proteomics-sample-metadata/tree/master/annotated-projects/PXD000612) |
|project | - | _json_ | [project](docs/project.rst) | -- |

> **Note**: Views can be extended and new views can be added to the format.

### Introduction to quantms.io

A quantms.io file is a collection of views, and they are aggregated into a folder `.qms` and inside that folder a file collect `project.json` MUST be present. Please read about the [project view](docs/project.rst) for more information. 

The introduction to the format, concepts and more details topics about serialization can be read in the introduction to the format [here](docs/introduction.rst).

## How to contribute

External contributors, researchers and the proteomics community are more than welcome to contribute to this project.

Contribute with the specification: you can contribute to the specification with ideas or refinements by adding an issue into the [issue tracker](https://github.com/bigbio/proteomics-quant-formats/issues) or performing a PR.

## Core contributors and collaborators

The project is run by different groups:

- Yasset Perez-Riverol (PRIDE Team, European Bioinformatics Institute - EMBL-EBI, U.K.)

IMPORTANT: If you contribute with the following specification, please make sure to add your name to the list of contributors.

## Code of Conduct

As part of our efforts toward delivering open and inclusive science, we follow the [Contributor Covenant Code of Conduct for Open Source Projects](https://www.contributor-covenant.org/version/2/0/code_of_conduct/).

## How to cite

## Copyright notice


    This information is free; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This information is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this work; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

