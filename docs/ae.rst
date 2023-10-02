Absolute expression format
==========================

The absolute expression format aims to cover the following use cases:

-  Fast and easy visualization absolute expression (AE) results using
   iBAQ values.
-  Store the AE results of each protein on each sample.
-  Provide information about the condition (factor value) of each sample
   for easy integration.
-  Store metadata information about the project, the workflow and the
   columns in the file.

Format
------

The absolute expression format by quantms is a tab-delimited file format
that contains the following fields:

-  ``protein`` -> Protein accession or semicolon-separated list of
   accessions for indistinguishable groups
-  ``sample_accession`` -> Sample accession in the SDRF.
-  ``condition`` -> Condition name
-  ``ibaq`` -> iBAQ value
-  ``ribaq`` -> Relative iBAQ value

Example:

=========== ================ ========= ====== =====
protein     sample_accession condition ibaq   ribaq
=========== ================ ========= ====== =====
LV861_HUMAN Sample-1         heart     1234.1 12.34
=========== ================ ========= ====== =====

AE Header
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
   #INFO=<ID=sample_accession, Number=1, Type=String, Description="Sample Accession in the SDRF">
   #INFO=<ID=condition, Number=1, Type=String, Description="Value of the factor value">
   #INFO=<ID=ibaq, Number=1, Type=Float, Description="Intensity based absolute quantification">
   #INFO=<ID=ribaq, Number=1, Type=Float, Description="relative iBAQ">

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
   (e.g. ``tissue``)

Please check also the differential expression example for more
information :doc:`de`
