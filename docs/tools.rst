quantms.io tools
=================================

project converter tool
-------------------------

.. code:: python

   python project_command.py
      --project_accession PXD014414
      --sdrf PXD014414.sdrf.tsv
      --quantms_version 1.12
      --output_folder result

DE converter tool
-------------------------

PRIDE project (make sure you have run the ``project_command.py``)

.. code:: python

   python differential_expression_command.py
      --msstats_file PXD014414.sdrf_openms_design_msstats_in_comparisons.csv
      --project_file result/PXD014414.json
      --sdrf_file PXD014414.sdrf.tsv
      --output_folder result


Non-PRIDE project(Don’t not need to run the ``project_command.py``)

.. code:: python

   python differential_expression_command.py
      --msstats_file PXD014414.sdrf_openms_design_msstats_in_comparisons.csv
      --generate_project False
      --sdrf_file PXD014414.sdrf.tsv
      --output_folder result

Optional parameter
~~~~~~~~~~~~~~~~~~

.. code:: python
   
   --fdr_threshold FDR threshold to use to filter the results
   --output_prefix_file Prefix of the df expression file
   --delete_existing Delete existing files in the output folder

AE converter tool
-------------------------

PRIDE project (make sure you have run the ``project_command.py``)

.. code:: python

   python absolute_expression_command.py
      --ibaq_file PXD004452-ibaq.csv
      --project_file result/PXD004452.json
      --output_folder result

Non-PRIDE project(Don’t not need to run the ``project_command.py``)

.. code:: python

   python absolute_expression_command.py
     --ibaq_file PXD004452-ibaq.csv
     --generate_project False
     --output_folder result


Optional parameter
~~~~~~~~~~~~~~~~~~

.. code:: python

   --output_prefix_file Prefix of the df expression file
   --delete_existing Delete existing files in the output folder


feature converter tool
-------------------------

-  pride projet (make sure you have run the ``project_command.py``)

.. code:: python

   python feature_command.py
      --sdrf_file PXD014414.sdrf.tsv
      --msstats_file PXD014414.sdrf_openms_design_msstats_in.csv
      --mztab_file PXD014414.sdrf_openms_design_openms.mzTab
      --output_folder result

-  Non-PRIDE project(Don’t not need to run the ``project_command.py``)

.. code:: python

   python feature_command.py
     --sdrf_file PXD014414.sdrf.tsv
     --msstats_file PXD014414.sdrf_openms_design_msstats_in.csv
     --mztab_file PXD014414.sdrf_openms_design_openms.mzTab
     --generate_project False
     --output_folder result

Optional parameter
~~~~~~~~~~~~~~~~~~

.. code:: python

   --use_cache Whether to use disk instead of memory.
   --output_prefix_file The prefix of the result file.
   --consensusxml_file The consensusXML file used to retrieve the mz/rt


psm converter tool
-------------------------

Note:: Make sure before generating the psm feature file that you generate the project.json
    
.. code:: python

   python psm_command.py convert-psm-file
      --mztab_file PXD014414.sdrf_openms_design_openms.mzTab
      --output_folder result

-  Non-PRIDE project(Don’t not need to run the ``project_command.py``)

.. code:: python

   python feature_command.py convert-psm-file
      --mztab_file PXD014414.sdrf_openms_design_openms.mzTab
      --generate_project False
      --output_folder result

Optional parameter
~~~~~~~~~~~~~~~~~~

.. code:: python

   --output_prefix_file The prefix of the result file.
   --verbose Output debug information.

compare psm.parquet
-------------------


.. code:: python

   python feature_command.py compare-set-of-psms
      --parquets PXD014414-comet.parquet PXD014414-sage.parquet PXD014414-msgf.parquet
      --tags comet sage msgf

generate spectra message
------------------------

generate_spectra_message support psm and parquet. Since the result file
is too large, you can specify –partition to split the result file.

.. code:: python

   python generate_spectra_message_command.py 
      --parquet_path PXD014414-f4fb88f6-0a45-451d-a8a6-b6d58fb83670.psm.parquet
      --mzml_directory mzmls
      --output_path psm/PXD014414.parquet
      --label psm
      --partition charge

map proteins accessions
------------------


get_unanimous_name support parquet and tsv. For parquet, map_parameter
have two option (map_protein_name or map_protein_accession), and the
label controls whether it is PSM or Feature.

-  parquet

.. code:: python

   python get_unanimous_command.py map-unanimous-for-parquet
      --parquet_path PXD014414-f4fb88f6-0a45-451d-a8a6-b6d58fb83670.psm.parquet
      --fasta Reference fasta database
      --output_path psm/PXD014414.psm.parquet
      --map_parameter map_protein_name
      --label psm


.. code:: python

   python get_unanimous_command.py get-unanimous-for-tsv
      --path PXD014414-c2a52d63-ea64-4a64-b241-f819a3157b77.differential.tsv
      --fasta Reference fasta database
      --output_path psm/PXD014414.de.tsv
      --map_parameter map_protein_name

compare two parquet files
--------------------------


.. code:: python

   python parquet_command.py
      --parquet_path_one res_lfq2_discache.parquet
      --parquet_path_two res_lfq2_no_cache.parquet
      --report_path report.txt

