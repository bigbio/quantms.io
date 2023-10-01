quantms.io tools
=================================
quantms.io tools provides a standardized 
set of commands to generate different files for your project.
You can generate separate files or complete project files depending on your needs.
A complete project file should contain 
``absolute_expression.tsv`` or ``differential_expression.tsv``, ``feature.parquet``, ``psm.parquet``, ``sdrf.tsv``.

Project converter tool
-------------------------
If your project comes from the PRIDE database, 
you can generate a ``project.json`` that contains 
descriptive information about the entire project.
For the project of the PRIDE database, our subsequent operations will be based on ``project.json``. 
So, if your project is from PRIDE, make sure you run ``project_command.py`` first.

- If you want to know more, please read :doc:`project`.
- If your project is not from PRIDE, you can skip this step.

.. code:: python

   python project_command.py
      --project_accession PXD014414
      --sdrf PXD014414.sdrf.tsv
      --quantms_version 1.12
      --output_folder result

DE converter tool
--------------------
Differential expression file 
Store the differential express proteins between two contrasts, 
with the corresponding fold changes and p-values.It can be easily visualized using tools such as 
`Volcano Plot <https://en.wikipedia.org/wiki/Volcano_plot_(statistics)>`__ and 
easily integrated with other omics data resources.

- If you want to know more, please read :doc:`de`.

-  PRIDE projet (make sure you have run the ``project_command.py``)

.. code:: python

   python differential_expression_command.py
      --msstats_file PXD014414.sdrf_openms_design_msstats_in_comparisons.csv
      --project_file result/PXD014414.json
      --sdrf_file PXD014414.sdrf.tsv
      --output_folder result


-  Non-PRIDE project(Don't not need to run the ``project_command.py``)

.. code:: python

   python differential_expression_command.py
      --msstats_file PXD014414.sdrf_openms_design_msstats_in_comparisons.csv
      --generate_project False
      --sdrf_file PXD014414.sdrf.tsv
      --output_folder result

- Optional parameter


.. code:: python
   
   --fdr_threshold   FDR threshold to use to filter the results(default 0.05)
   --output_prefix_file   Prefix of the df expression file(like {prefix}-{uu.id}-{extension})
   --delete_existing   Delete existing files in the output folder(default True)

AE converter tool
--------------------
The absolute expression format aims to visualize absolute expression (AE) results using
iBAQ values and store the AE results of each protein on each sample.

- If you want to know more, please read :doc:`ae`.

-  PRIDE projet (make sure you have run the ``project_command.py``)
.. code:: python

   python absolute_expression_command.py
      --ibaq_file PXD004452-ibaq.csv
      --project_file result/PXD004452.json
      --output_folder result

-  Non-PRIDE project(Don't not need to run the ``project_command.py``)

.. code:: python

   python absolute_expression_command.py
     --ibaq_file PXD004452-ibaq.csv
     --generate_project False
     --output_folder result


- Optional parameter

.. code:: python

   --output_prefix_file    Prefix of the df expression file(like {prefix}-{uu.id}-{extension})
   --delete_existing    Delete existing files in the output folder(default True)


Feature converter tool
-------------------------
The Peptide table aims to cover detail on peptide level including peptide intensity. 
The most of content are from peptide part of mzTab. 
It store peptide intensity to perform down-stream analysis and integration.

- If you want to know more, please read :doc:`feature`.

In some projects, mzTab files can be very large, so we provide both ``diskcache`` and ``no-diskcache`` versions of the tool. 
You can choose the desired version according to your server configuration.

-  PRIDE projet (make sure you have run the ``project_command.py``)

.. code:: python

   python feature_command.py
      --sdrf_file PXD014414.sdrf.tsv
      --msstats_file PXD014414.sdrf_openms_design_msstats_in.csv
      --mztab_file PXD014414.sdrf_openms_design_openms.mzTab
      --output_folder result


-  Non-PRIDE project(Don't not need to run the ``project_command.py``)

.. code:: python

   python feature_command.py
     --sdrf_file PXD014414.sdrf.tsv
     --msstats_file PXD014414.sdrf_openms_design_msstats_in.csv
     --mztab_file PXD014414.sdrf_openms_design_openms.mzTab
     --generate_project False
     --output_folder result

- Optional parameter

.. code:: python

   --use_cache    Whether to use diskcache instead of memory(default True)
   --output_prefix_file   The prefix of the result file(like {prefix}-{uu.id}-{extension})
   --consensusxml_file   The consensusXML file used to retrieve the mz/rt(default None)


Psm converter tool
---------------------
The PSM table aims to cover detail on PSM level for AI/ML training and other use-cases.
It store details on PSM level including spectrum mz/intensity for specific use-cases such as AI/ML training.

- If you want to know more, please read :doc:`psm`.

-  PRIDE projet(make sure you have run the ``project_command.py``)
    
.. code:: python

   python psm_command.py convert-psm-file
      --mztab_file PXD014414.sdrf_openms_design_openms.mzTab
      --output_folder result

-  Non-PRIDE project(Don't not need to run the ``project_command.py``)

.. code:: python

   python feature_command.py convert-psm-file
      --mztab_file PXD014414.sdrf_openms_design_openms.mzTab
      --generate_project False
      --output_folder result

- Optional parameter

.. code:: python

   --output_prefix_file   The prefix of the result file(like {prefix}-{uu.id}-{extension})
   --verbose  Output debug information(default True)

Compare psm.parquet
-------------------
This tool is used to compare peptide information in result files obtained by different search engines.


.. code:: python

   python feature_command.py compare-set-of-psms
      --parquets PXD014414-comet.parquet PXD014414-sage.parquet PXD014414-msgf.parquet
      --tags comet sage msgf

Generate spectra message
-------------------------

generate_spectra_message support psm and feature. It can be used directly for spectral clustering.

- ``--label`` contains two options: ``psm`` and ``feature``.
- ``--partion`` contains two options: ``charge`` and ``reference_file_name``.
Since the result file is too large, you can specify –partition to split the result file.

.. code:: python

   python generate_spectra_message_command.py 
      --parquet_path PXD014414-f4fb88f6-0a45-451d-a8a6-b6d58fb83670.psm.parquet
      --mzml_directory mzmls
      --output_path psm/PXD014414.parquet
      --label psm
      --partition charge

Map proteins accessions
------------------------

get_unanimous_name support parquet and tsv. For parquet, map_parameter
have two option (``map_protein_name`` or ``map_protein_accession``), and the
label controls whether it is PSM or Feature.

-  parquet

.. code:: python

   python get_unanimous_command.py map-unanimous-for-parquet
      --parquet_path PXD014414-f4fb88f6-0a45-451d-a8a6-b6d58fb83670.psm.parquet
      --fasta Reference fasta database
      --output_path psm/PXD014414.psm.parquet
      --map_parameter map_protein_name
      --label psm

- tsv 
.. code:: python

   python get_unanimous_command.py get-unanimous-for-tsv
      --path PXD014414-c2a52d63-ea64-4a64-b241-f819a3157b77.differential.tsv
      --fasta Reference fasta database
      --output_path psm/PXD014414.de.tsv
      --map_parameter map_protein_name

Compare two parquet files
--------------------------
This tool is used to compare the feature.parquet file generated by two versions (``diskcache`` or ``no-diskcache``).

.. code:: python

   python parquet_command.py
      --parquet_path_one res_lfq2_discache.parquet
      --parquet_path_two res_lfq2_no_cache.parquet
      --report_path report.txt

