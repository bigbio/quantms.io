# quantms.io
[![Python application](https://github.com/bigbio/quantms.io/actions/workflows/python-app.yml/badge.svg?branch=dev)](https://github.com/bigbio/quantms.io/actions/workflows/python-app.yml)

[quantms](https://docs.quantms.org) is a nf-core pipeline for the analysis of quantitative proteomics data. The pipeline is based on the [OpenMS](https://www.openms.de/) framework and [DIA-NN](https://github.com/vdemichev/DiaNN); and it is designed to analyze large scale experiments. the main outputs of quantms tools are the following: 

- [mzTab](https://github.com/HUPO-PSI/mzTab) files with the identification and quantification information.
- [MSstats](https://msstats.org/wp-content/uploads/2017/01/MSstats_v3.7.3_manual.pdf) input file with the peptide quantification values needed for the MSstats analysis.
- [MSstats](https://msstats.org/wp-content/uploads/2017/01/MSstats_v3.7.3_manual.pdf) output file with the differential expression values for each protein. 
- The input [SDRF](https://github.com/bigbio/proteomics-sample-metadata) of the pipeline if available. 

Here, we aim to formalize and develop a more standardized format that enables better representation of the identification and quantification results but also enables new and novel use cases for proteomics data analysis: 

- Fast and easy visualization of the identification and quantification results.
- Easy integration with other omics data.
- Easy integration with sample metadata.
- AI/ML model development based on identification and quantification results.

**Note**: We are not trying to replace the mzTab format, but to provide a new format that enables AI-related use cases. Most of the features of the mzTab format will be included in the new format.  

## Data model

The GitHub repository aims to provide multiple formats for serialization of the data model, including:

- `Tab-delimited` format similar to mzTab. 
- `JSON` format to enable integration with other bioinformatics resources. 
- `Parquet` format to enable integration with big data frameworks and large-scale data integration. 

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

