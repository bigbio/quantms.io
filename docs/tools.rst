quantms.io tools
=================================

quantms.io tools provides a standardized set of commands to generate different files for your project. You can generate separate files or complete project files depending on your needs. A completed project contains the following files:

- ``project.json`` -> contains descriptive information about the entire project.

Example:

.. code:: json

   {
    "project_accession": "PXD014414",
    "project_title": "",
    "project_sample_description": "",
    "project_data_description": "",
    "project_pubmed_id": 32265444,
    "organisms": [
        "Homo sapiens"
    ],
    "organism_parts": [
        "mammary gland",
        "adjacent normal tissue"
    ],
    "diseases": [
        "metaplastic breast carcinomas",
        "Triple-negative breast cancer",
        "Normal",
        "not applicable"
    ],
    "cell_lines": [
        "not applicable"
    ],
    "instruments": [
        "Orbitrap Fusion"
    ],
    "enzymes": [
        "Trypsin"
    ],
    "experiment_type": [
        "Triple-negative breast cancer",
        "Wisp3",
        "Tandem mass tag (tmt) labeling",
        "Ccn6",
        "Metaplastic breast carcinoma",
        "Precision therapy",
        "Lc-ms/ms shotgun proteomics"
    ],
    "acquisition_properties": [
        {"proteomics data acquisition method": "TMT"},
        {"proteomics data acquisition method": "Data-dependent acquisition"},
        {"dissociation method": "HCD"},
        {"precursor mass tolerance": "20 ppm"},
        {"fragment mass tolerance": "0.6 Da"}
    ],
    "quantms_files": [
        {"feature_file": "PXD014414-943a8f02-0527-4528-b1a3-b96de99ebe75.featrue.parquet"},
        {"sdrf_file": "PXD014414-f05eca35-9381-40d8-a7da-2fe57745afaf.sdrf.tsv"},
        {"psm_file": "PXD014414-f4fb88f6-0a45-451d-a8a6-b6d58fb83670.psm.parquet"},
        {"differential_file": "PXD014414-3026e5d5-fb0e-45e9-a4f0-c97d86536716.differential.tsv"}
    ],
    "quantms_version": "1.1.1",
    "comments": []
   }


- ``absolute_expression.tsv`` or ``differential_expression.tsv``
The differential expression format by quantms is based on the
`MSstats <https://msstats.org/wp-content/uploads/2017/01/MSstats_v3.7.3_manual.pdf>`__
output.

Example:

+---------+-------------------------+-----+----+---+----+-------+----+
| protein | label                   | log | se | d | pv | adj.p | i  |
|         |                         | 2fc |    | f | al | value | ss |
|         |                         |     |    |   | ue |       | ue |
+=========+=========================+=====+====+===+====+=======+====+
| LV86    | normal-squamous cell    | 0   | 0. | 8 | 0. | 0.62  | NA |
| 1_HUMAN | carcinoma               | .60 | 87 |   | 51 |       |    |
+---------+-------------------------+-----+----+---+----+-------+----+

The absolute expression format by quantms contains IBAQ message.

Example:

=========== ================ ========= ====== =====
protein     sample_accession condition ibaq   ribaq
=========== ================ ========= ====== =====
LV861_HUMAN Sample-1         heart     1234.1 12.34
=========== ================ ========= ====== =====

- ``feature.parquet``
The ``feature.parquet`` cover detail on peptide level.

Example: 

+------------+-------------------------------------------------------------------------+-------------------------+-----------------------+-----------------------+--------+--------------------------------+----------------+--------+--------------------+---------------------+-------------------------------+-----------------------------+---------------+----------+-----------+----------------+---------------------+-----------+----------+----------------------+--------------+--------------------+-------+---------+------------------------------------------------------------+---------------------+------------------------------+----------------------+----------+-----------------+-----------+-----------------+------------+
| sequence   | protein_accessions                                                      | protein_start_positions | protein_end_positions | protein_global_qvalue | unique | modifications                  | retention_time | charge | exp_mass_to_charge | calc_mass_to_charge | peptidoform                   | posterior_error_probability | global_qvalue | is_decoy | intensity | spectral_count | sample_accession    | condition | fraction | biological_replicate | fragment_ion | isotope_label_type | run   | channel | id_scores                                                  | reference_file_name | best_psm_reference_file_name | best_psm_scan_number | mz_array | intensity_array | num_peaks | gene_accessions | gene_names |
+============+=========================================================================+=========================+=======================+=======================+========+================================+================+========+====================+=====================+===============================+=============================+===============+==========+===========+================+=====================+===========+==========+======================+==============+====================+=======+=========+============================================================+=====================+==============================+======================+==========+=================+===========+=================+============+
| ASPDWGYDDK | ['sp|CONTAMINANT_P00915|CONTAMINANT_CAH1_HUMAN','sp|P00915|CAH1_HUMAN'] | [1 2]                   | [10 11]               | 0.001882796           | 0      | ['0-UNIMOD:1' '10-UNIMOD:737'] | 7522.223146    | 2      | 712.831298         | 712.8302134         | [Acetyl]-ASPDWGYDDK[TMT6plex] | 4.97E-05                    | 0             | 0        | 454585.3  | 1              | PXD014414-Sample-10 | Norm      | 1        | 10                   | None         | L                  | 1_1_1 | TMT131  | ["'OpenMS:Best PSM Score':0.0",'Best PSM PEP:4.96872e-05'] | UM_F_50cm_2019_0414 | UM_F_50cm_2019_0430          | 53434                |          |                 |           |                 |            |
+------------+-------------------------------------------------------------------------+-------------------------+-----------------------+-----------------------+--------+--------------------------------+----------------+--------+--------------------+---------------------+-------------------------------+-----------------------------+---------------+----------+-----------+----------------+---------------------+-----------+----------+----------------------+--------------+--------------------+-------+---------+------------------------------------------------------------+---------------------+------------------------------+----------------------+----------+-----------------+-----------+-----------------+------------+

- ``psm.parquet``

``psm.parquet`` store details on PSM level including spectrum mz/intensity for specific use-cases such as AI/ML training.

Example: 

+----------+--------------------------+-------------------------+-----------------------+-----------------------+--------+------------------+----------------+--------+--------------------+---------------------+------------------+-----------------------------+---------------+----------+-------------------------------------------------------------------------------------------------+-------------------+---------------------+-------------+----------+-----------------+-----------+-----------------+------------+
| sequence | protein_accessions       | protein_start_positions | protein_end_positions | protein_global_qvalue | unique | modifications    | retention_time | charge | exp_mass_to_charge | calc_mass_to_charge | peptidoform      | posterior_error_probability | global_qvalue | is_decoy | id_scores                                                                                       | consensus_support | reference_file_name | scan_number | mz_array | intensity_array | num_peaks | gene_accessions | gene_names |
+==========+==========================+=========================+=======================+=======================+========+==================+================+========+====================+=====================+==================+=============================+===============+==========+=================================================================================================+===================+=====================+=============+==========+=================+===========+=================+============+
| SSPGHR   | ['sp|P29692|EF1D_HUMAN'] | [118]                   | [123]                 | 0.001882796           | 1      | ['1-UNIMOD:737'] | 1258.2         | 2      | 435.2432855        | 435.2431809         | S[TMT6plex]SPGHR | 0.35875                     |               | 0        | ["'OpenMS:Target-decoy PSM q-value': 0.040626999360205",'Posterior error probability: 0.35875'] |                   | UM_F_50cm_2019_0428 | 2193        |          |                 |           |                 |            |
+----------+--------------------------+-------------------------+-----------------------+-----------------------+--------+------------------+----------------+--------+--------------------+---------------------+------------------+-----------------------------+---------------+----------+-------------------------------------------------------------------------------------------------+-------------------+---------------------+-------------+----------+-----------------+-----------+-----------------+------------+

- ``sdrf.tsv``

``sdrf.tsv`` is a file used by quantMS to search the library.

Example: 

+--------------------+---------------------------+--------------------------------+--------------------------------------+-------------------------------+-------------------------------------+----------------------+----------------------+----------------------------+----------------------------+---------------------------------------+-----------------------------+---------------+------------+------------------------------------------+----------------+-------------------------+----------------------------------------------------------------------------------------+------------------------------+------------------------------+---------------------------------+----------------------------------+--------------------------------------------+-----------------------------------------------------+----------------------------------------------+-----------------------------------------+---------------------------------------------------------+--------------------------------------------+------------------------------+---------------------------+-----------------------------------+----------------------------------+-------------------------------+
| source name        | characteristics[organism] | characteristics[organism part] | characteristics[developmental stage] | characteristics[disease]      | characteristics[histologic subtype] | characteristics[sex] | characteristics[age] | characteristics[cell type] | characteristics[cell line] | characteristics[biological replicate] | characteristics[individual] | Material Type | assay name | Technology Type                          | comment[label] | comment[data file]      | comment[file uri]                                                                      | comment[technical replicate] | comment[fraction identifier] | comment[cleavage agent details] | comment[instrument]              | comment[modification parameters]           | comment[modification parameters]                    | comment[modification parameters]             | comment[modification parameters]        | comment[modification parameters]                        | comment[modification parameters]           | comment[dissociation method] | comment[collision energy] | comment[precursor mass tolerance] | comment[fragment mass tolerance] | factor value[disease]         |
+====================+===========================+================================+======================================+===============================+=====================================+======================+======================+============================+============================+=======================================+=============================+===============+============+==========================================+================+=========================+========================================================================================+==============================+==============================+=================================+==================================+============================================+=====================================================+==============================================+=========================================+=========================================================+============================================+==============================+===========================+===================================+==================================+===============================+
| PXD014414-Sample-1 | Homo sapiens              | mammary gland                  | adult                                | metaplastic breast carcinomas | Chondroid                           | female               | 43Y                  | not applicable             | not applicable             | 1                                     | C1                          | tissue        | run 1      | proteomic profiling by mass spectrometry | TMT126         | UM_F_50cm_2019_0414.raw | ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/04/PXD014414/UM_F_50cm_2019_0414.raw | 1                            | 1                            | AC=MS:1001251;NT=Trypsin        | NT=Orbitrap Fusion;AC=MS:1002416 | NT=Oxidation;MT=Variable;TA=M;AC=UNIMOD:35 | NT=Acetyl;AC=UNIMOD:1;PP=Protein N-term;MT=variable | NT=Carbamidomethyl;TA=C;MT=fixed;AC=UNIMOD:4 | NT=TMT6plex;AC=UNIMOD:737;TA=K;MT=Fixed | NT=TMT6plex;AC=UNIMOD:737;PP=Protein N-term;MT=Variable | NT=TMT6plex;AC=UNIMOD:737;TA=S;MT=Variable | NT=HCD;AC=PRIDE:0000590      | 55 NCE                    | 20 ppm                            | 0.6 Da                           | metaplastic breast carcinomas |
+--------------------+---------------------------+--------------------------------+--------------------------------------+-------------------------------+-------------------------------------+----------------------+----------------------+----------------------------+----------------------------+---------------------------------------+-----------------------------+---------------+------------+------------------------------------------+----------------+-------------------------+----------------------------------------------------------------------------------------+------------------------------+------------------------------+---------------------------------+----------------------------------+--------------------------------------------+-----------------------------------------------------+----------------------------------------------+-----------------------------------------+---------------------------------------------------------+--------------------------------------------+------------------------------+---------------------------+-----------------------------------+----------------------------------+-------------------------------+

- If you want see a full example, please click `here <https://github.com/bigbio/quantms.io/tree/main/python/quantmsio/quantms_io/data>`__

Project converter tool
-------------------------
If your project comes from the PRIDE database, 
you can generate a ``project.json`` that contains 
descriptive information about the entire project.

- If you want to know more, please read :doc:`project`.
- If your project is not from PRIDE, you can skip this step.

.. code:: python

   python project_command.py generate_pride_project_json
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

- If you have generated project.json, you can use this parameter ``--project_file`` to add project information for DE files.
- If you want to know more, please read :doc:`de`.

Example: 

.. code:: python

   python differential_expression_command.py convert_msstats_differential
      --msstats_file PXD014414.sdrf_openms_design_msstats_in_comparisons.csv
      --sdrf_file PXD014414.sdrf.tsv
      --output_folder result

- Optional parameter

.. code:: python

   --project_file   Descriptive information from project.json(project json path)
   --fdr_threshold   FDR threshold to use to filter the results(default 0.05)
   --output_prefix_file   Prefix of the df expression file(like {prefix}-{uu.id}-{extension})
   --delete_existing   Delete existing files in the output folder(default True)

AE converter tool
--------------------
The absolute expression format aims to visualize absolute expression (AE) results using
iBAQ values and store the AE results of each protein on each sample.

- If you have generated project.json, you can use this parameter ``--project_file`` to add project information for AE files.
- If you want to know ibaq, please read `ibaqpy <https://github.com/bigbio/ibaqpy>`__
- If you want to know more, please read :doc:`ae`.

Example: 

.. code:: python

   python absolute_expression_command.py attach_file_to_json
      --ibaq_file PXD004452-ibaq.csv
      --sdrf_file PXD014414.sdrf.tsv
      --output_folder result

- Optional parameter

.. code:: python

   --project_file   Descriptive information from project.json(project json path)
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

Example: 

.. code:: python

   python feature_command.py convert_feature_file
      --sdrf_file PXD014414.sdrf.tsv
      --msstats_file PXD014414.sdrf_openms_design_msstats_in.csv
      --mztab_file PXD014414.sdrf_openms_design_openms.mzTab
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

Example: 
    
.. code:: python

   python psm_command.py convert_psm_file
      --mztab_file PXD014414.sdrf_openms_design_openms.mzTab
      --output_folder result


- Optional parameter

.. code:: python

   --output_prefix_file   The prefix of the result file(like {prefix}-{uu.id}-{extension})
   --verbose  Output debug information(default True)

DiaNN convert 
--------------------------
For DiaNN, the command supports generating ``feature.parquet`` and ``psm.parquet`` directly from diann_report files.

- ``--modifications`` is a list of 2 lengths containing both fixed and variable modifications. The different modifications in each modification are separated by ``,``.

Example: 

.. code:: python

   python diann_convert_command.py diann_convert_to_parquet
      --report_path diann_report.tsv
      --design_file PXD037682.sdrf_openms_design.tsv
      --modifications "Carbamidomethyl (C)" "null"
      --qvalue_threshold 0.05
      --mzml_info_folder mzml
      --sdrf_path PXD037682.sdrf.tsv
      --output_folder result
      --output_prefix_file PXD037682
      --threads 60


Compare psm.parquet
-------------------
This tool is used to compare peptide information in result files obtained by different search engines.

- ``--tags`` or ``-t`` are used to specify the tags of the PSM table.

Example: 

.. code:: python

   python feature_command.py compare_set_of_psms
      -p PXD014414-comet.parquet
      -p PXD014414-sage.parquet
      -p PXD014414-msgf.parquet
      -t comet
      -t sage
      -t msgf

Generate spectra message
-------------------------

generate_spectra_message support psm and feature. It can be used directly for spectral clustering.

- ``--label`` contains two options: ``psm`` and ``feature``.
- ``--partion`` contains two options: ``charge`` and ``reference_file_name``.
Since the result file is too large, you can specify ``–-partition`` to split the result file.

Example: 

.. code:: python

   python generate_spectra_message_command.py map_spectrum_message_to_parquet
      --parquet_path PXD014414-f4fb88f6-0a45-451d-a8a6-b6d58fb83670.psm.parquet
      --mzml_directory mzmls
      --output_path psm/PXD014414.parquet
      --label psm
      --chunksize(default 100000)
      --partition charge

Map proteins accessions
------------------------

get_unanimous_name support parquet and tsv. For parquet, map_parameter
have two option (``map_protein_name`` or ``map_protein_accession``), and the
label controls whether it is PSM or Feature.

-  parquet
- ``--label`` contains two options: ``psm`` and ``feature``

Example: 

.. code:: python

   python get_unanimous_command.py map_unanimous_for_parquet
      --parquet_path PXD014414-f4fb88f6-0a45-451d-a8a6-b6d58fb83670.psm.parquet
      --fasta Reference fasta database
      --output_path psm/PXD014414.psm.parquet
      --map_parameter map_protein_name
      --label psm

- tsv
  
Example: 

.. code:: python

   python get_unanimous_command.py get_unanimous_for_tsv
      --path PXD014414-c2a52d63-ea64-4a64-b241-f819a3157b77.differential.tsv
      --fasta Reference fasta database
      --output_path psm/PXD014414.de.tsv
      --map_parameter map_protein_name

Compare two parquet files
--------------------------
This tool is used to compare the feature.parquet file generated by two versions (``diskcache`` or ``no-diskcache``).

Example: 

.. code:: python

   python parquet_command.py compare_two_parquet
      --parquet_path_one res_lfq2_discache.parquet
      --parquet_path_two res_lfq2_no_cache.parquet
      --report_path report.txt

Generate report about files 
-----------------------------
This tool is used to generate report about all feature files or psm files.
You can build ``psm parquet`` or ``feature parquet`` multiple times for the same project and use this command to verify its consistency.

- ``--label`` contains two options: ``psm`` and ``feature``

Example: 

.. code:: python

   python generate_report_command.py generate_report_about_files
      --check_dir file_path
      --label psm

Register file 
--------------------------
This tool is used to register the file to ``project.json``.
If your project comes from the PRIDE database, You can use this command to add file information for ``project.json``.

- The parameter ``--category`` has three options: ``feature_file``, ``psm_file``, ``differential_file``, ``absolute_file``.You can add the above file types.
- The parameter ``--replace_existing`` is enable then we remove the old file and add this one. If not then we can have a list of files for a category.

Example: 

.. code:: python
   
   python attach_file_command.py attach_file_to_json
      --project_file PXD014414/project.json
      --attach_file PXD014414-943a8f02-0527-4528-b1a3-b96de99ebe75.featrue.parquet
      --category feature_file
      --replace_existing

Data preview
--------------------------
This tool is used to preview your feature files and AE files.
You can run ``streamlit run .\visualize_web_commond.py`` start a web service.
Then set up your working directory to preview the data.

.. image:: data_view.png
   :width: 800
   :align: center

* If you want to manipulate data on NoteBook, you can introduce the ``Statistic`` class.

.. code:: python

   from quantms_io.core.statistic import Statistic