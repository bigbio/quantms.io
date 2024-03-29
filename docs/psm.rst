PSM table format
================

Use cases
---------

-  The PSM table aims to cover detail on PSM level for AI/ML training
   and other use-cases.
-  Most of the content is similar to mzTab, a PSM would be a peptide
   identification in an specific msrun file.
-  The representation should be in a parquet file. \*
-  Store details on PSM level including spectrum mz/intensity for
   specific use-cases such as AI/ML training.
-  Fast and easy visualization and scanning on PSM level.

*\**\ `Parquet <https://github.com/apache/parquet-format>`__ file is a
columnar storage format that supports nested data. For these large-scale
analyses, Parquet has helped its users reduce storage requirements by at
least one-third on large datasets, in addition, it greatly improved scan
and deserialization time (web use-cases), hence the overall costs. The
following table compares the savings as well as the speedup obtained by
converting data into Parquet from CSV.

+----------+----------------------+------------------+----------------+
| Dataset  | Size on Amazon S3    | Query Run Time   | Data Scanned   |
+==========+======================+==================+================+
| Data     | 1 TB                 | 236 seconds      | 1.15 TB        |
| stored   |                      |                  |                |
| as CSV   |                      |                  |                |
| files    |                      |                  |                |
+----------+----------------------+------------------+----------------+
| Data     | 130 GB               | 6.78 seconds     | 2.51 GB        |
| stored   |                      |                  |                |
| in       |                      |                  |                |
| Apache   |                      |                  |                |
| Parquet  |                      |                  |                |
| Format   |                      |                  |                |
+----------+----------------------+------------------+----------------+

Format
------

The `Avro PSM schema <psm.avsc>`__ is used to define the PSM table
format. The following table describes the fields of the PSM table.

-  ``sequence``: The peptide’s sequence corresponding to the PSM -> ``string``
-  ``protein_accessions``: A list protein’s accessions -> ``list[string]``
-  ``protein_start_positions``: A list of protein’s start positions -> ``list[int]``
-  ``protein_end_positions``: A list of protein’s end positions -> ``list[int]``
-  ``protein_global_qvalue``: The global q-value of the associated protein or protein group -> ``double``
-  ``unique``: Indicates whether the peptide sequence (coming from the PSM) is unique for this protein in respect to the searched database -> ``boolean (0/1)``
-  ``modifications``: A list of modifications for a give peptide ``[modification1, modification2, ...]``. A modification should be recorded as string like `modification definition <README.md#modifications>`__-> ``list[string]``
-  ``retention_time``: The retention time of the spectrum -> ``float``
-  ``charge``: The charge assigned by the search engine/software -> ``integer``
-  ``exp_mass_to_charge``: The PSM’s experimental mass to charge (m/z) -> ``double``
-  ``calc_mass_to_charge``: The PSM’s calculated (theoretical) mass to charge (m/z) -> ``double``
-  ``reference_file_name``: The reference file name that contains the spectrum. -> ``string``
-  ``scan_number``: The scan number of the spectrum. The scan number or index of the spectrum in the file -> ``string``
-  ``peptidoform``: Peptidoform of the PSM. See more `documentation here <README.md#peptidoform>`__. -> ``string``
-  ``posterior_error_probability``: Posterior Error Probability score from quantms -> ``double``
-  ``global_qvalue``: Global q-value from quantms -> ``double``
-  ``is_decoy``: Indicates whether the peptide sequence (coming from the PSM) is decoy -> ``boolean (0/1)``

Optional fields:

-  ``gene_accessions``: A list of gene accessions -> ``list[string]``
-  ``gene_names``: A list of gene names -> ``list[string]``
-  ``consensus_support``: Global consensus support scores for multiple search engines -> ``float``
-  ``mz_array``: A list of mz values for the spectrum -> ``list[double]``
-  ``intensity_array``: A list of intensity values for the spectrum -> ``list[float]``
-  ``num_peaks``: The number of peaks in the spectrum, this is the size of previous lists intensity and mz -> ``integer``
-  ``id_scores``: A list of identification scores, search engine, percolator etc. Each search engine score will be a key/value pair ``(e.g. "MS-GF:RawScore": 78.9)`` -> ``list[string]``
