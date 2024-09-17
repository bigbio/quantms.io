PSM_MAP = {
    "sequence": "sequence",
    "modifications": "modifications",
    "opt_global_Posterior_Error_Probability_score": "posterior_error_probability",
    "opt_global_q-value": "global_qvalue",
    "opt_global_cv_MS:1002217_decoy_peptide": "is_decoy",
    "calc_mass_to_charge": "calculated_mz",
    "accession": "pg_accessions",
    "unique": "unique",
    "charge": "precursor_charge",
    "exp_mass_to_charge": "observed_mz",
    "retention_time": "rt",
}
PSM_USECOLS = list(PSM_MAP.keys()) + [
    "spectra_ref",
    "start",
    "end",
]

MSSTATS_MAP = {
    "ProteinName": "pg_accessions",
    "Reference": "reference_file_name",
    "Intensity": "intensity",
    "Channel": "channel",
    "RetentionTime": "rt",
    "PeptideSequence": "peptidoform"
}
MSSTATS_USECOLS = set(MSSTATS_MAP.keys())