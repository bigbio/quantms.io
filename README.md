# quantms.io

The proteomics quantification formats, is a Github repository that aims to formalize some existing proteomics quantification formats for large scale quantification experiments. Different from previous efforts like (mzTab, quantML, etc), the present representation aims to reuse existing popular formats and work in the following use cases: 

**Note**: Before starting, for a more generic/extended MS format for quantitative proteomics please check the [mzTab](https://github.com/HUPO-PSI/mzTab) format. mzTab is a generic format for MS data, including quantification data for proteomics and metabolomics experiments. Aim to capture not only the quantitative information but also, the identification information, including the peptide spectrum matches (psms), post-translational modifications, etc.   

## Why a new format?

Why other efforts like mzTab, quantML, have been developed to represent quantitative proteomics data, we believe those formats are not enough to represent the following information, and also fails to handle the following cases: 

- QuantML and mzTab are design for DDA experiments, with lower number of ms_runs.  
   - In mzTab, when the number of ms_runs increases, the number of columns in the peptide table with null values increases making difficult the reliability of the format.
   - In mzTab, when the number of ms_runs and samples increases, the metadata section increases making difficult the reliability of the format. for each ms_run, at least 5 sections are needed, in an experiment with 1000 ms_runs, 5000 sections are needed. 
   - QuantML was never designed to handle large scale experiments, and the format is not flexible enough to handle large scale experiments.
- mzTab is a large format making difficult to handle and visualize quantitative information as simple as: 
   - Different expression results tables. 
   - Raw intensities tables at peptide/protein information.
- Sample metadata integration with [SDRF-Proteomics](https://github.com/bigbio/proteomics-sample-metadata) format is not possible.

More important, both formats and previous efforts do not provide enough tooling framework to enable bioinformatic software packages, main reason why the formats are yet popular within the bioinformatics community. 

## Gols and Use cases

The main goals of this repository are:

- Provide a data model to represent quantitative proteomics data for absolute quantification and differential expression experiments.
- Provide a data model to represent quantitative proteomics data for DIA and DDA experiments.
- Provide a data model to represent quantitative proteomics data for large scale experiments.
- Provide a data model that enable integration with the Sample to Data Relationship Format (SDRF-Proteomics) for proteomics experiments.
- Provide a data model to represent protein and peptide quantification data.

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

