Differential expression format
==============================

Use cases
---------

-  Store the differential express proteins between two contrasts, with
   the corresponding fold changes and p-values.
-  Enable easy visualization using tools like `Volcano
   Plot <https://en.wikipedia.org/wiki/Volcano_plot_(statistics)>`__.
-  Enable easy integration with other omics data resources.
-  Store metadata information about the project, the workflow and the
   columns in the file.

Format
------

The differential expression format by quantms is based on the
`MSstats <https://msstats.org/wp-content/uploads/2017/01/MSstats_v3.7.3_manual.pdf>`__
output. The MSstats format is a tab-delimited file that contains the
following fields - see example `file <../include/PXD004683.csv>`__:

-  ``protein`` -> Protein Accession
-  ``label`` -> Label for the contrast on which the fold changes and
   p-values are based on
-  ``log2fc`` -> Log2 Fold Change
-  ``se`` -> Standard error of the log2 fold change
-  ``df`` -> Degree of freedom of the Student test
-  ``pvalue`` -> Raw p-values
-  ``adj.pvalue`` -> P-values adjusted among all the proteins in the
   specific comparison using the approach by Benjamini and Hochberg
-  ``issue`` -> Issue column shows if there is any issue for inference
   in corresponding protein and comparison, for example,
   OneConditionMissing or CompleteMissing.

Example:

+---------+-------------------------+-----+----+---+----+-------+----+
| protein | label                   | log | se | d | pv | adj.p | i  |
|         |                         | 2fc |    | f | al | value | ss |
|         |                         |     |    |   | ue |       | ue |
+=========+=========================+=====+====+===+====+=======+====+
| LV86    | normal-squamous cell    | 0   | 0. | 8 | 0. | 0.62  | NA |
| 1_HUMAN | carcinoma               | .60 | 87 |   | 51 |       |    |
+---------+-------------------------+-----+----+---+----+-------+----+

DE Header
---------

By default, the MSstats format does not have any header of metadata. We
suggest adding a header to the output for better understanding of the
file. By default, MSstats allows comments in the file if the line starts
with ``#``. The quantms output will start with some key value pairs that
describe the project, the workflow and also the columns in the file. For
example:

``#project_accession=PXD000000``

In addition, for each ``Default`` column of the matrix the following
information should be added:

::

   #INFO=<ID=protein, Number=inf, Type=String, Description="Protein Accession">
   #INFO=<ID=label, Number=1, Type=String, Description="Label for the Conditions combination">
   #INFO=<ID=log2fc, Number=1, Type=Double, Description="Log2 Fold Change">
   #INFO=<ID=se, Number=1, Type=Double, Description="Standard error of the log2 fold change">
   #INFO=<ID=df, Number=1, Type=Integer, Description="Degree of freedom of the Student test">
   #INFO=<ID=pvalue, Number=1, Type=Double, Description="Raw p-values">
   #INFO=<ID=adj.pvalue, Number=1, Type=Double, Description="P-values adjusted among all the proteins in the specific comparison using the approach by Benjamini and Hochberg">
   #INFO=<ID=issue, Number=1, Type=String, Description="Issue column shows if there is any issue for inference in corresponding protein and comparison">

-  The ``ID`` is the column name in the matrix, the ``Number`` is the
   number of values in the column (separated by ``;``), the ``Type`` is
   the type of the values in the column and the ``Description`` is a
   description of the column. The number of values in the column can go
   from 1 to ``inf`` (infinity).
-  Protein groups are written as a list of protein accessions separated
   by ``;`` (e.g. ``P12345;P12346``)

We suggest including the following properties in the header:

-  project_accession: The project accession in PRIDE Archive
-  project_title: The project title in PRIDE Archive
-  project_description: The project description in PRIDE Archive
-  quanmts_version: The version of the quantms workflow used to generate
   the file
-  factor_value: The factor values used in the analysis
   (e.g. ``phenotype``)
-  fdr_threshold: The FDR threshold used to filter the protein lists
   (e.g. ``adj.pvalue < 0.05``)

A complete example of a quantms output file can be seen
`here <../include/PXD004683-quantms.csv>`__.
