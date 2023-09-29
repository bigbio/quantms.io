Generate Psm.parquet or compare psm.parquet
===========================================

Generate Psm.parquet
--------------------

Use cases
~~~~~~~~~

-  PRIDE projet(make sure you have run the ``project_command.py``)

::

   python psm_command.py convert-psm-file
   --mztab_file PXD014414.sdrf_openms_design_openms.mzTab
   --output_folder result

-  Non-PRIDE project(Don’t not need to run the ``project_command.py``)

::

   python feature_command.py convert-psm-file
   --mztab_file PXD014414.sdrf_openms_design_openms.mzTab
   --generate_project False
   --output_folder result

Optional parameter
~~~~~~~~~~~~~~~~~~

-  –-output_prefix_file The prefix of the result file.
-  –-verbose Output debug information.

Compare psm.parquet
-------------------

Use case
~~~~~~~~

::

   python feature_command.py compare-set-of-psms
   --parquets PXD014414-comet.parquet PXD014414-sage.parquet PXD014414-msgf.parquet
   --tags comet sage msgf
