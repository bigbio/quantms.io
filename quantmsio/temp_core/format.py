import pyarrow as pa
PEPTIDE_FIELDS = [
    pa.field(
        "sequence",
        pa.string(),
        metadata={"description": "The peptide’s sequence corresponding to the PSM"},
    ),
    pa.field(
        "peptidoform",
        pa.string(),
        metadata={"description": "Peptide sequence with modifications: Read the specification for more details"},
    ),
    pa.field(
        "modifications",
        pa.list_(pa.string()),
        metadata={"description": "List of modifications as string array, easy for search and filter"},
    ),
    pa.field(
        "modification_details",
        pa.list_(pa.string()),
        metadata={"description": "List of alternative site probabilities for the modification format: read the specification for more details"},
    ),
    pa.field(
        "posterior_error_probability",
        pa.float32(),
        metadata={"description": "Posterior error probability for the given peptide spectrum match"},
    ),
    pa.field("global_qvalue", pa.float32(), metadata={"description": "Global q-value of the peptide or psm at the level of the experiment"}),
    
    pa.field(
        "is_decoy",
        pa.int32(),
        metadata={"description": "Decoy indicator, 1 if the PSM is a decoy, 0 target"},
    ),
    pa.field(
        "calculated_mz",
        pa.float32(),
        metadata={"description": "Theoretical peptide mass-to-charge ratio based on identified sequence and modifications"},
    ),
    pa.field(
        "additional_scores",
        pa.list_(
            pa.struct([
                ("name", pa.string()),
                ("value", pa.float32()) 
            ])
        ),
        metadata={"description": "List of structures, each structure contains two fields: name and value"},
    ),

    pa.field(
        "pg_accessions",
        pa.list_(pa.string()),
        metadata={"description": "Protein group accessions of all the proteins that the peptide maps to"},
    ),
    pa.field(
        "pg_positions",
        pa.list_(pa.string()),
        metadata={"description": "Protein start and end positions written as start_post:end_post"},
    ),
    pa.field(
        "unique",
        pa.int32(),
        metadata={"description": "Unique peptide indicator, if the peptide maps to a single protein, the value is 1, otherwise 0"},
    ),
    pa.field(
        "protein_global_qvalue",
        pa.float32(),
        metadata={"description": "Global q-value of the protein group at the experiment level"},
    ),
    pa.field(
        "gg_accessions",
        pa.list_(pa.string()),
        metadata={"description": "Gene accessions, as string array"},
    ),
    pa.field(
        "gg_names",
        pa.list_(pa.string()),
        metadata={"description": "Gene names, as string array"},
    ),

    pa.field(
        "precursor_charge",
        pa.int32(),
        metadata={"description": "Precursor charge"},
    ),
    pa.field(
        "observed_mz",
        pa.float32(),
        metadata={"description": "Experimental peptide mass-to-charge ratio of identified peptide (in Da)"},
    ),
    pa.field(
        "rt",
        pa.float32(),
        metadata={"description": "MS2 scan’s precursor retention time (in seconds)"},
    ),
    pa.field(
        "predicted_rt",
        pa.float32(),
        metadata={"description": "Predicted retention time of the peptide (in seconds)"},
    ),
    pa.field(
        "quantmsio_version",
        pa.string(),
        metadata={"description": "The version of quantms.io"},
    )
]

PSM_UNIQUE_FIELDS = [
    pa.field(
        "reference_file_name",
        pa.string(),
        metadata={"description": "Spectrum file name with no path information and not including the file extension"},
    ),
    pa.field(
        "scan_number",
        pa.string(),
        metadata={"description": "Scan number of the spectrum"},
    ),
    pa.field(
        "ion_mobility",
        pa.float32(),
        metadata={"description": "Ion mobility value for the precursor ion"},
    ),
    pa.field("num_peaks", pa.int32(), metadata={"description": "Number of peaks in the spectrum used for the peptide spectrum match"}),
    pa.field(
        "mz_array",
        pa.list_(pa.float32()),
        metadata={"description": "Array of m/z values for the spectrum used for the peptide spectrum match"},
    ),
    pa.field(
        "intensity_array",
        pa.list_(pa.float32()),
        metadata={"description": "Array of intensity values for the spectrum used for the peptide spectrum match"},
    ),
    pa.field("rank", pa.int32(), metadata={"description": "Rank of the peptide spectrum match in the search engine output"}),

    pa.field(
        "cv_params",
        pa.list_(
            pa.struct([
                ("name", pa.string()),
                ("value", pa.string()) 
            ])
        ),
        metadata={"description": "Optional list of CV parameters for additional metadata"},
    ),
]
















PSM_FIELDS = PEPTIDE_FIELDS + PSM_UNIQUE_FIELDS