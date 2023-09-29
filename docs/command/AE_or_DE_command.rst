Generate AE or DE file
======================

DE
--

use case
~~~~~~~~

-  PRIDE projet(make sure you have run the ``project_command.py``)

::

   python differential_expression_command.py
   --msstats_file PXD014414.sdrf_openms_design_msstats_in_comparisons.csv
   --project_file result/PXD014414.json
   --sdrf_file PXD014414.sdrf.tsv
   --output_folder result

-  Non-PRIDE project(Don’t not need to run the ``project_command.py``)

::

   python differential_expression_command.py
   --msstats_file PXD014414.sdrf_openms_design_msstats_in_comparisons.csv
   --generate_project False
   --sdrf_file PXD014414.sdrf.tsv
   --output_folder result

Optional parameter
~~~~~~~~~~~~~~~~~~

-  –fdr_threshold FDR threshold to use to filter the results
-  –output_prefix_file Prefix of the df expression file
-  –delete_existing Delete existing files in the output folder

AE
--

.. _use-case-1:

use case
~~~~~~~~

-  PRIDE projet(make sure you have run the ``project_command.py``)

::

   python absolute_expression_command.py
   --ibaq_file PXD004452-ibaq.csv
   --project_file result/PXD004452.json
   --output_folder result

-  Non-PRIDE project(Don’t not need to run the ``project_command.py``)

::

   python absolute_expression_command.py
   --ibaq_file PXD004452-ibaq.csv
   --generate_project False
   --output_folder result

.. _optional-parameter-1:

Optional parameter
~~~~~~~~~~~~~~~~~~

-  –output_prefix_file Prefix of the df expression file
-  –delete_existing Delete existing files in the output folder
