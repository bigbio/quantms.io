from quantmsio import __version__
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

SDRF_MAP = {
    "comment[data file]": "reference",
    "comment[label]": "label",
    "source name": "sample_accession",
    "comment[fraction identifier]": "fraction",
    "characteristics[biological replicate]": "biological_replicate",
}

SDRF_USECOLS = set( list(SDRF_MAP.keys()) + [
    "comment[technical replicate]"
])

DIANN_MAP = {
    'File.Name': "reference_file_name",
    "Precursor.Normalised": "intensity",
    "RT.Start": "rt_start",
    'RT.Stop': 'rt_stop',
    "RT": 'rt',
    'Predicted.RT': "predicted_rt",
    'Protein.Ids': 'pg_accessions',
    "PEP": "posterior_error_probability",
    "Global.Q.Value": "global_qvalue",
    "Global.PG.Q.Value": "protein_global_qvalue",
    "Precursor.Charge": "precursor_charge",
    "Stripped.Sequence": "sequence",
    "Modified.Sequence": "peptidoform",
    "Genes":"gg_names",
    "opt_global_spectrum_reference": "scan_number",
    "Calculate.Precursor.Mz": "calculated_mz",
    "exp_mass_to_charge": "observed_mz"
}

QUANTMSIO_VERSION = __version__
          


