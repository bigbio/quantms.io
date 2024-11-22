==== quantms.io tools

quantms.io tools provides a standardized set of commands to generate different files for your project. It is mainly used to consolidate the data of each file and generate a standardized file representation. 

image::map.png[width=80%]
   :width: 80%
   :align: center

You can generate separate files or complete project files depending on your needs.A completed project contains the following files:

Project converter tool
-------------------------
If your project comes from the PRIDE database, 
you can use the pride accession to generate a `project.json` that contains 
descriptive information about the entire project.
Or, customize a Project Accession to generate an entirely new project.
You can create a `protein.txt` to generate specific project of the protein information in it.
It's like this:

.. code:: shell

   #protein.txt
   Q07878
   O43660
   P63261

* If you want to know more, please read :doc:`project`.
* If your project is not from PRIDE, you can skip this step.

.. code:: shell

   quantmsioc generate-pride-project-json
      --project_accession PXD014414
      --sdrf PXD014414.sdrf.tsv
      --output_folder result

* Optional parameter

.. code:: shell
   --software_name   software name used to generate the data
   --software_version   software version used to generate the data
   --delete_existing   Delete existing files in the output folder(default False)

DE converter tool
--------------------
Differential expression file 
Store the differential express proteins between two contrasts, 
with the corresponding fold changes and p-values.It can be easily visualized using tools such as 
`Volcano Plot <https://en.wikipedia.org/wiki/Volcano_plot_(statistics)>`__ and 
easily integrated with other omics data resources.

* If you have generated project.json, you can use this parameter `--project_file` to add project information for DE files.
* If you want to know more, please read :doc:`de`.

Example: 

.. code:: shell

   quantmsioc convert-de
      --msstats_file PXD014414.sdrf_openms_design_msstats_in_comparisons.csv
      --sdrf_file PXD014414.sdrf.tsv
      --output_folder result

* Optional parameter

.. code:: shell

   --project_file   Descriptive information from project.json(project json path)
   --protein_file   Protein file that meets specific requirements(protein.txt)
   --fdr_threshold   FDR threshold to use to filter the results(default 0.05)
   --output_prefix_file   Prefix of the df expression file(like {prefix}-{uu.id}-{extension})
   --delete_existing   Delete existing files in the output folder(default True)

AE converter tool
--------------------
The absolute expression format aims to visualize absolute expression (AE) results using
iBAQ values and store the AE results of each protein on each sample.

* If you have generated project.json, you can use this parameter `--project_file` to add project information for AE files.
* If you want to know ibaq, please read `ibaqpy <https://github.com/bigbio/ibaqpy>`__
* If you want to know more, please read :doc:`ae`.

Example: 

.. code:: shell

   quantmsioc convert-ae
      --ibaq_file PXD004452-ibaq.csv
      --sdrf_file PXD014414.sdrf.tsv
      --output_folder result

* Optional parameter

.. code:: shell

   --project_file   Descriptive information from project.json(project json path)
   --protein_file   Protein file that meets specific requirements(protein.txt)
   --output_prefix_file    Prefix of the df expression file(like {prefix}-{uu.id}-{extension})
   --delete_existing    Delete existing files in the output folder(default True)


Feature converter tool
-------------------------
The Peptide table aims to cover detail on peptide level including peptide intensity. 
The most of content are from peptide part of mzTab. 
It store peptide intensity to perform down-stream analysis and integration.

* If you want to know more, please read :doc:`feature`.

Mztab
>>>>>>>

Example: 

.. code:: shell

   quantmsioc convert-feature
      --sdrf_file PXD014414.sdrf.tsv
      --msstats_file PXD014414.sdrf_openms_design_msstats_in.csv
      --mztab_file PXD014414.sdrf_openms_design_openms.mzTab
      --output_folder result

* Optional parameter

.. code:: shell

   --file_num Read batch size(default 50)
   --protein_file   Protein file that meets specific requirements(protein.txt)
   --partitions The field used for splitting files, multiple fields are separated by `,`
   --output_prefix_file   The prefix of the result file(like {prefix}-{uu.id}-{extension})
   --duckdb_max_memory  The maximum amount of memory allocated by the DuckDB engine (e.g 16GB)
   --duckdb_threads The number of threads for the DuckDB engine (e.g 4)

Maxquant
>>>>>>>>>

Example: 

.. code:: shell

   quantmsioc convert-maxquant-feature
      --evidence_file evidence.txt
      --sdrf_file PXD014414.sdrf.tsv
      --output_folder result

* Optional parameter

.. code:: shell
   --protein_file   Protein file that meets specific requirements(protein.txt)
   --partitions The field used for splitting files, multiple fields are separated by `,`
   --chunksize Read batch size
   --output_prefix_file Prefix of the parquet file needed to generate the file name


DiaNN
>>>>>>>

Example: 

.. code:: shell

   quantmsioc convert-diann
      --report_path diann_report.tsv
      --qvalue_threshold 0.05
      --mzml_info_folder mzml
      --sdrf_path PXD037682.sdrf.tsv
      --output_folder result

* Optional parameter

.. code:: shell
   --protein_file Protein file that meets specific requirements
   --partitions The field used for splitting files, multiple fields are separated by `,`
   --output_prefix_file Prefix of the Json file needed to generate the file name
   --duckdb_max_memory  The maximum amount of memory allocated by the DuckDB engine (e.g 16GB)
   --duckdb_threads The number of threads for the DuckDB engine (e.g 4)
   --file_num Read batch size(default 50)


Psm converter tool
---------------------

The PSM table aims to cover detail on PSM level for AI/ML training and other use-cases.
It store details on PSM level including spectrum mz/intensity for specific use-cases such as AI/ML training.

* If you want to know more, please read :doc:`psm`.

Mztab
>>>>>>>

Example: 
    
.. code:: shell

   quantmsioc convert-psm
      --mztab_file PXD014414.sdrf_openms_design_openms.mzTab
      --output_folder result

* Optional parameter

.. code:: shell

   --protein_file   Protein file that meets specific requirements(protein.txt)
   --chunksize Read batch size
   --output_prefix_file   The prefix of the result file(like {prefix}-{uu.id}-{extension})

Maxquant
>>>>>>>>>

Example: 
    
.. code:: shell

   quantmsioc convert-maxquant-psm
      --msms_file the msms.txt file, this will be used to extract the peptide information
      --output_folder result

* Optional parameter

.. code:: shell

   --chunksize Read batch size
   --output_prefix_file The prefix of the result file(like {prefix}-{uu.id}-{extension})

Compare psm.parquet
-------------------
This tool is used to compare peptide information in result files obtained by different search engines.

* `--tags` or `-t` are used to specify the tags of the PSM table.

Example: 

.. code:: shell

   quantmsioc compare-set-psms
      -p PXD014414-comet.parquet
      -p PXD014414-sage.parquet
      -p PXD014414-msgf.parquet
      -t comet
      -t sage
      -t msgf

Generate spectra message
-------------------------
generate_spectra_message support psm. It can be used directly for spectral clustering.


Since the result file is too large, you can specify `–-partitions` to split the result file.

Example: 

.. code:: shell

   quantmsioc map-spectrum-message-to-parquet
      --parquet_path PXD014414-f4fb88f6-0a45-451d-a8a6-b6d58fb83670.psm.parquet
      --mzml_directory mzmls
      --output_folder result

* Optional parameter

.. code:: shell

   --file_num The number of rows of parquet read using pandas streaming
   --partitions The field used for splitting files, multiple fields are separated by `,`

Generate gene message
-------------------------
generate_gene_message support feature. 

Example: 

.. code:: shell

   quantmsioc map-gene-msg-to-parquet 
   --parquet_path PXD000672-0beee055-ae78-4d97-b6ac-1f191e91bdd4.featrue.parquet
   --fasta_path Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta
   --output_folder result

* Optional parameter

.. code:: shell

   --file_num The number of rows of parquet read using pandas streaming
   --partitions The field used for splitting files, multiple fields are separated by `,`
   --species species type(default human)

* `species`
  
+-------------+-------------------------+
| Common name |       Genus name        |
+=============+=========================+
|    human    |      Homo sapiens       |
+-------------+-------------------------+
|    mouse    |      Mus musculus       |
+-------------+-------------------------+
|     rat     |    Rattus norvegicus    |
+-------------+-------------------------+
|  fruitfly   | Drosophila melanogaster |
+-------------+-------------------------+
|  nematode   | Caenorhabditis elegans  |
+-------------+-------------------------+
|  zebrafish  |       Danio rerio       |
+-------------+-------------------------+
| thale-cress |  Arabidopsis thaliana   |
+-------------+-------------------------+
|    frog     |   Xenopus tropicalis    |
+-------------+-------------------------+
|     pig     |       Sus scrofa        |
+-------------+-------------------------+



Register file 
--------------------------
This tool is used to register the file to `project.json`.
If your project comes from the PRIDE database, You can use this command to add file information for `project.json`.

* The parameter `--category` has three options: `sdrf_file`, `feature_file`, `psm_file`, `differential_file`, `absolute_file`.You can add the above file types.
* The parameter `--replace_existing` is enable then we remove the old file and add this one. If not then we can have a list of files for a category.

Example: 

.. code:: shell
   
   quantmsioc attach-file
      --project_file PXD014414/project.json
      --attach_file PXD014414-943a8f02-0527-4528-b1a3-b96de99ebe75.featrue.parquet
      --category feature_file


* Optional parameter

.. code:: shell 
   --is_folder A boolean value that indicates if the file is a folder or not
   --replace_existing Whether to delete old files
   --partitions The fields that are used to partition the data in the file. This is used to optimize the data retrieval and filtering of the data. This field is optional.

Statistics
-----------
This tool is used for statistics.
Example: 

.. code:: shell

   quantmsioc project-ae-statistics
      --absolute_path PXD010154-51b34353-227f-4d38-a181-6d42824de9f7.absolute.tsv
      --parquet_path PXD010154-51b34353-227f-4d38-a181-6d42824de9f7.featrue.parquet
      --save_path PXD014414.statistic.txt

.. code:: shell

   quantmsioc parquet-psm-statistics
      --parquet_path PXD010154-51b34353-227f-4d38-a181-6d42824de9f7.psm.parquet
      --save_path PXD014414.statistic.txt

Plots
-------
This tool is used for visualization.
* plot-psm-peptides
  
.. code:: shell

   quantmsioc plot plot-psm-peptides
      --psm_parquet_path PXD010154-51b34353-227f-4d38-a181-6d42824de9f7.psm.parquet
      --sdrf_path PXD010154.sdrf.tsv
      --save_path PXD014414_psm_peptides.svg

* plot-ibaq-distribution
  
.. code:: shell

   quantmsioc plot plot-ibaq-distribution
      --ibaq_path PXD010154-51b34353-227f-4d38-a181-6d42824de9f7.ibaq.tsv
      --select_column IbaqLog
      --save_path PXD014414_psm_peptides.svg

* plot-kde-intensity-distribution

.. code:: shell

      quantmsioc plot plot-kde-intensity-distribution
      --feature_path PXD010154-51b34353-227f-4d38-a181-6d42824de9f7.featrue.parquet
      --num_samples 10
      --save_path PXD014414_psm_peptides.svg

* plot-bar-peptide-distribution

.. code:: shell

      quantmsioc plot plot-bar-peptide-distribution
      --feature_path PXD010154-51b34353-227f-4d38-a181-6d42824de9f7.featrue.parquet
      --num_samples 10
      --save_path PXD014414_psm_peptides.svg

* plot-box-intensity-distribution

.. code:: shell

      quantmsioc plot plot-box-intensity-distribution
      --feature_path PXD010154-51b34353-227f-4d38-a181-6d42824de9f7.featrue.parquet
      --num_samples 10
      --save_path PXD014414_psm_peptides.svg


