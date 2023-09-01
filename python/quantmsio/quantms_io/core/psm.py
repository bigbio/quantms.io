import numpy as np
import pandas as pd

from quantms_io.core.mztab import MztabHandler
import pyarrow as pa

from quantms_io.core.parquet_handler import ParquetHandler
from quantms_io.utils.pride_utils import standardize_protein_list_accession, get_quantmsio_modifications


class PSMHandler(ParquetHandler):

    PSM_FIELDS = [pa.field("sequence", pa.string(),
                               metadata={"description": "Peptide sequence of the feature"}),
                      pa.field("protein_accessions", pa.list_(pa.string()),
                               metadata={"description": "accessions of associated proteins"}),
                      pa.field("protein_start_positions", pa.list_(pa.int32()),
                               metadata={"description": "start positions in the associated proteins"}),
                      pa.field("protein_end_positions", pa.list_(pa.int32()),
                               metadata={"description": "end positions in the associated proteins"}),
                      pa.field("protein_global_qvalue", pa.float64(),
                               metadata={"description": "global q-value of the associated protein or protein group"}),
                      pa.field("unique", pa.int32(),
                               metadata={"description": "if the peptide is unique to a particular protein"}),
                      pa.field("modifications", pa.list_(pa.string()),
                               metadata={"description": "peptide modifications"}),
                      pa.field("retention_time", pa.float64(),
                               metadata={"description": "retention time"}),
                      pa.field("charge", pa.int32(),
                               metadata={"description": "charge state of the feature"}),
                      pa.field("exp_mass_to_charge", pa.float64(),
                               metadata={"description": "experimentally measured mass-to-charge ratio"}),
                      pa.field("calc_mass_to_charge", pa.float64(),
                               metadata={"description": "calculated mass-to-charge ratio"}),
                      pa.field("peptidoform", pa.string(),
                               metadata={"description": "peptidoform in proforma notation"}),
                      pa.field("posterior_error_probability", pa.float64(),
                               metadata={"description": "posterior error probability"}),
                      pa.field("global_qvalue", pa.float64(),
                               metadata={"description": "global q-value"}),
                      pa.field("is_decoy", pa.int32(),
                               metadata={
                                   "description": "flag indicating if the feature is a decoy (1 is decoy, 0 is not decoy)"}),
                      pa.field("id_scores", pa.list_(pa.string()),
                               metadata={"description": "identification scores as key value pairs"}),
                      pa.field("consensus_support", pa.float32(),
                                metadata={"description": "consensus support value"}),
                      pa.field("reference_file_name", pa.string(),
                               metadata={"description": "file name of the reference file for the feature"}),
                      pa.field("scan_number", pa.string(),
                               metadata={"description": "scan number of the best PSM"}),
                      pa.field("mz_array", pa.list_(pa.float64()),
                               metadata={"description": "mass-to-charge ratio values"}),
                      pa.field("intensity_array", pa.list_(pa.float64()),
                               metadata={"description": "intensity array values"}),
                      pa.field("num_peaks", pa.int32(),
                               metadata={"description": "number of peaks"}),
                      pa.field("gene_accessions", pa.list_(pa.string()),
                               metadata={"description": "accessions of associated genes"}),
                      pa.field("gene_names", pa.list_(pa.string()),
                               metadata={"description": "names of associated genes"}),
                      ]

    def __init__(self) -> None:
        super().__init__()
        self.parquet_path = None

    def _create_schema(self):
        """
        Create the schema for the psm file. The schema is defined in the docs folder of this repository.
        (https://github.com/bigbio/quantms.io/blob/main/docs/PSM.md)
        """
        return pa.schema(PSMHandler.PSM_FIELDS, metadata={"description": "PSM file in quantms.io format"})

    def _create_psm_table(self, psm_list: list) -> pa.Table:
        return pa.Table.from_pandas(pd.DataFrame(psm_list), schema=self.schema)

    def convert_mztab_to_feature(self, mztab_path: str, parquet_path: str = None):
        """
        Convert a mzTab file to a feature file
        :param mztab_path: path to the mzTab file
        :return: path to the feature file
        """
        if parquet_path is not None:
            self.parquet_path = parquet_path

        mztab_handler = MztabHandler(mztab_file=mztab_path)
        mztab_handler.create_mztab_psm_iterator(mztab_path)
        psm_list = []
        for it in iter(mztab_handler.read_next_psm, None):
            psm_list.append(self._transform_psm_from_mztab(psm=it, mztab_handler=mztab_handler))
            print(it["sequence"] + "---" + it["accession"])

        feature_table = self._create_psm_table(psm_list)
        self.write_single_file_parquet(feature_table, parquet_output=self.parquet_path, write_metadata=True)

    def _transform_psm_from_mztab(self, psm, mztab_handler) -> dict:
        """
        Transform a an mztab psm to quantms io psm.
        :param psm mztab psm
        :param mztab_handler the mztab handler with all information about scores, ms_runs, modification
        :return: dictionary of psm following the schema defined in PSMHandler.PSM_FIELDS
        """

        sequence = psm["sequence"]
        protein_accessions = standardize_protein_list_accession(psm["accession"])
        protein_start_positions = [int(i) for i in psm["start"].split(",")]
        protein_end_positions = [int(i) for i in psm["end"].split(",")]
        unique = 1 if len(protein_accessions) == 1 else 0

        protein_accession_nredundant = list(dict.fromkeys(protein_accessions))
        protein_qvalue = mztab_handler.get_protein_qvalue_from_index_list(protein_accession_list=protein_accession_nredundant)
        protein_qvalue = None if protein_qvalue is None else np.float64(protein_qvalue[0])

        retention_time = None if ("retention_time" not in psm or psm["retention_time"] is None) else np.float64(psm["retention_time"])
        charge = int(psm["charge"])
        calc_mass_to_charge = None if ("calc_mass_to_charge" not in psm or psm["calc_mass_to_charge"] is None) \
            else np.float64(psm["calc_mass_to_charge"])
        exp_mass_to_charge = None if ("exp_mass_to_charge" not in psm or psm["exp_mass_to_charge"] is None) \
            else np.float64(psm["exp_mass_to_charge"])

        modifications_string = psm["modifications"]  # Mods
        modifications = get_quantmsio_modifications(modifications_string=modifications_string,
                                                    modification_definition=mztab_handler.get_modifications_definition())
        modifications_string = ""
        for key, value in modifications.items():
            modifications_string += "|".join(map(str, value["position"]))
            modifications_string = modifications_string + "-" + value["unimod_accession"] + ","
        modifications_string = None if len(modifications_string) == 0 else modifications_string[
                                                                           :-1]  # Remove last comma
        modification_list = None if modifications_string is None else modifications_string.split(",")
        posterior_error_probability = None if ("posterior_error_probability" not in psm
                                               or psm["posterior_error_probability"] is None) \
            else np.float64(psm["posterior_error_probability"])


        global_qvalue = None if ("global_qvalue" not in psm
                                               or psm["global_qvalue"] is None) \
            else np.float64(psm["global_qvalue"])

        consensus_support = None if psm["consensus_support"] else np.float32(psm["consensus_support"])

        psm_score = np.float64(psm["score"])
        peptide_score_name = mztab_handler.get_search_engine_scores()["psm_score"]

        return {
            "sequence": sequence,
            "protein_accessions": protein_accessions,
            "protein_start_positions": protein_start_positions,
            "protein_end_positions": protein_end_positions,
            "protein_global_qvalue": protein_qvalue,
            "unique": unique,
            "modifications": modification_list,
            "retention_time": retention_time,
            "charge": charge,
            "calc_mass_to_charge": calc_mass_to_charge,
            "peptidoform": psm["proforma_peptidoform"],  # Peptidoform in proforma notation
            "posterior_error_probability": posterior_error_probability,
            "global_qvalue": global_qvalue,
            "is_decoy": int(psm["is_decoy"]),
            "consensus_support": consensus_support,
            "id_scores": [f"{peptide_score_name}: {psm_score}", f"Posterior error probability: {posterior_error_probability}"],
            "reference_file_name": psm["ms_run"],
            "scan_number": psm["scan_number"],
            "exp_mass_to_charge": exp_mass_to_charge,
            "mz_array": None,
            "intensity_array": None,
            "num_peaks": None,
            "gene_accessions": None,
            "gene_names": None
        }


