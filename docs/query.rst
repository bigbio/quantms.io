Query parquet
=============

The query module provides the ability to quickly search through Parquet files.

Basic query
-------------------------
Basic query operations allow you to query all ``samples``, ``peptides``, ``proteins``, ``genes``, and ``MZML`` files. The query results will be deduplicated and then returned in a list format.

.. code:: python

    from quantmsio.core.query import Parquet
    P = Parquet('PXD007683.feature.parquet')
    P.get_unique_samples()
    """
    ['PXD007683-Sample-3',
    'PXD007683-Sample-9',
    'PXD007683-Sample-4',
    'PXD007683-Sample-6',
    'PXD007683-Sample-11',
    'PXD007683-Sample-7',
    'PXD007683-Sample-1',
    'PXD007683-Sample-10',
    'PXD007683-Sample-8',
    'PXD007683-Sample-2',
    'PXD007683-Sample-5']
    """
    P.get_unique_peptides()
    P.get_unique_proteins()
    P.get_unique_genes()
    P.get_unique_references()


Specific query
----------------
Specific queries allow you to individually search for certain values based on specific conditions.
The results are returned in the form of a ``DataFrame``.

.. code:: python

    P.query_peptide('QPAYVSK') 
    """
        sequence    protein_accessions	    protein_start_positions ...
    	QPAYVSK	    [sp|P36016|LHS1_YEAST]	[739]
        QPAYVSK	    [sp|P36016|LHS1_YEAST]	[739]
        QPAYVSK	    [sp|P36016|LHS1_YEAST]	[739]
        ...
    """
    P.query_peptide('QPAYVSK',columns=['protein_start_positions','protein_end_positions']) 

.. code:: python

    P.query_peptides(['QPAYVSK','QPCPSQYSAIK'],columns=None) 
    """
        sequence    protein_accessions	    protein_start_positions ...
    	QPAYVSK	    [sp|P36016|LHS1_YEAST]	[739]
        QPAYVSK	    [sp|P36016|LHS1_YEAST]	[739]
        QPAYVSK	    [sp|P36016|LHS1_YEAST]	[739]
        ...
        QPCPSQYSAIK [sp|O95861|BPNT1_HUMAN]	[98]
        ...
    """

.. code:: python

    P.query_protein('P36016',columns=None)
    P.query_proteins(['P36016','O95861'],columns=None)
    """
        sequence    protein_accessions	    protein_start_positions ...
    	QPAYVSK	    [sp|P36016|LHS1_YEAST]	[739]
        QPAYVSK	    [sp|P36016|LHS1_YEAST]	[739]
        QPAYVSK	    [sp|P36016|LHS1_YEAST]	[739]
        ...
        QPCPSQYSAIK [sp|O95861|BPNT1_HUMAN]	[98]
        ...
    """

.. code:: python
    
    P.get_samples_from_database(['PXD007683-Sample-3','PXD007683-Sample-9'],columns=None)
    """
    sequence                protein_accessions          sample_accession
    AAAAAAAAAAAAAAAGAGAGAK  [sp|P55011|S12A2_HUMAN]	PXD007683-Sample-3
    AAAAAAAAAAAAAAAGAGAGAK  [sp|P55011|S12A2_HUMAN]	PXD007683-Sample-3
    AAAAAAAAAK	            [sp|Q99453|PHX2B_HUMAN]	PXD007683-Sample-3
    """
    P.get_report_from_database(['a05063','a05059'],columns=None) # mzml
    """
    sequence    protein_accessions      reference_file_name
    AAAAAAALQAK [sp|P36578|RL4_HUMAN]   a05063
    AAAAAAALQAK [sp|P36578|RL4_HUMAN]   a05063
    AAAAAAALQAK [sp|P36578|RL4_HUMAN]   a05063
    """

Iter bacth
----------------
You can use the following method to produce values in batches.

.. code-block:: python
    :linenos:

    for samples,df in P.iter_samples(file_num=10,columns=None):
        # A batch contains ten samples.
        print(samples,df)

    for df in P.iter_chunk(batch_size=500000,columns=None):
        # A batch contains 500,000 rows.
        print(df)
    
    for refs,df in P.iter_file(file_num=20,columns=None): # mzml
        # A batch contains 20 mzML files.
        print(refs,df)

Inject message
----------------
You can use the following method to fill in additional information.

.. code-block:: python
    :linenos:
    
    df = P.get_report_from_database(['a05063','a05059'],columns=None)
    df = P.inject_spectrum_msg(df, mzml_directory='./mzml')
    fasta = './Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta'
    protein_dict = P.get_protein_dict(fasta_path=fasta)
    df = P.inject_position_msg(df, protein_dict)
    df = P.inject_gene_msg(df,fasta,map_parameter = "map_protein_accession",species = "human")