from quantmsio import __version__
from quantmsio.core.format import PSM_FIELDS, FEATURE_FIELDS, IBAQ_FIELDS
import pyarrow as pa

PSM_MAP = {
    "sequence": "sequence",
    "modifications": "modifications",
    "opt_global_cv_MS:1000889_peptidoform_sequence": "peptidoform",
    "opt_global_q-value": "global_qvalue",
    "opt_global_consensus_support": "consensus_support",
    "opt_global_cv_MS:1002217_decoy_peptide": "is_decoy",
    "calc_mass_to_charge": "calculated_mz",
    "accession": "mp_accessions",
    "charge": "precursor_charge",
    "exp_mass_to_charge": "observed_mz",
    "retention_time": "rt",
}
PSM_USECOLS = list(PSM_MAP.keys()) + ["spectra_ref"]
PEP = [
    "opt_global_Posterior_Error_Probability_score",
    "opt_global_Posterior_Error_Probability",
    "opt_global_MS:1001493_score",
]

MSSTATS_MAP = {
    "ProteinName": "pg_accessions",
    "Reference": "reference_file_name",
    "Intensity": "intensity",
    "Channel": "channel",
    "RetentionTime": "rt",
    "PeptideSequence": "peptidoform",
}

MSSTATS_USECOLS = set(MSSTATS_MAP.keys())

SDRF_MAP = {
    "comment[data file]": "reference",
    "comment[label]": "label",
    "source name": "sample_accession",
    "comment[fraction identifier]": "fraction",
    "characteristics[biological replicate]": "biological_replicate",
}

SDRF_USECOLS = set(list(SDRF_MAP.keys()) + ["comment[technical replicate]"])

DIANN_MAP = {
    "File.Name": "reference_file_name",
    "Precursor.Quantity": "intensity",
    "RT.Start": "rt_start",
    "RT.Stop": "rt_stop",
    "RT": "rt",
    "Predicted.RT": "predicted_rt",
    "Protein.Group": "pg_accessions",
    "Protein.Ids": "mp_accessions",
    "PEP": "posterior_error_probability",
    "Global.Q.Value": "global_qvalue",
    "Global.PG.Q.Value": "pg_global_qvalue",
    "Q.Value": "qvalue",
    "PG.Q.Value": "pg_qvalue",
    "Precursor.Normalised": "normalize_intensity",
    "PG.MaxLFQ": "lfq",
    "Quantity.Quality": "quantitative_score",
    "Precursor.Charge": "precursor_charge",
    "Stripped.Sequence": "sequence",
    "Modified.Sequence": "peptidoform",
    "Genes": "gg_names",
    "Run": "run"
}
DIANN_USECOLS = list(DIANN_MAP.keys())

MAXQUANT_PSM_MAP = {
    "Sequence": "sequence",
    "Proteins": "mp_accessions",
    "PEP": "posterior_error_probability",
    "Modified sequence": "peptidoform",
    "Reverse": "is_decoy",
    "m/z": "observed_mz",
    "Scan number": "scan",
    "Retention time": "rt",
    "Charge": "precursor_charge",
    "Modified sequence": "peptidoform",
    "Raw file": "reference_file_name",
    "Score": "andromeda_score",
    "Delta score": "andromeda_delta_score",
    "PIF": "parent_ion_score",
}

MAXQUANT_FEATURE_MAP = {
    "Sequence": "sequence",
    "Proteins": "mp_accessions",
    "Leading proteins": "pg_accessions",
    "Gene names": "gg_names",
    "PEP": "posterior_error_probability",
    "Modified sequence": "peptidoform",
    "Charge": "precursor_charge",
    "Modified sequence": "peptidoform",
    "Raw file": "reference_file_name",
    "Score": "andromeda_score",
    "Delta score": "andromeda_delta_score",
    "PIF": "parent_ion_score",
    "Reverse": "is_decoy",
    "m/z": "observed_mz",
    "MS/MS scan number": "scan",
    "Calibrated retention time": "rt",
    "Calibrated retention time start": "rt_start",
    "Calibrated retention time finish": "rt_stop",
}

IBAQ_USECOLS = [
    "pg_accessions",
    "peptidoform",
    "sequence",
    "precursor_charge",
    "intensities",
    "reference_file_name",
    "unique",
]

MAXQUANT_PSM_USECOLS = list(MAXQUANT_PSM_MAP.keys())

MAXQUANT_FEATURE_USECOLS = list(MAXQUANT_FEATURE_MAP.keys())

QUANTMSIO_VERSION = __version__

PSM_SCHEMA = pa.schema(
    PSM_FIELDS,
    metadata={"description": "psm file in quantms.io format"},
)
FEATURE_SCHEMA = pa.schema(
    FEATURE_FIELDS,
    metadata={"description": "feature file in quantms.io format"},
)

IBAQ_SCHEMA = pa.schema(
    IBAQ_FIELDS,
    metadata={"description": "ibaq file in quantms.io format"},
)
