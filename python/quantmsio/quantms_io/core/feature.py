"""
The feature file handle and manage the feature table in column format. The main serialization format is Apache Parquet.
The feature file is defined in the docs folder of this repository.
(https://github.com/bigbio/quantms.io/blob/main/docs/FEATURE.md). Among other information, the feature file contains:
    - Peptide sequence
    - Peptide modifications
    - Peptide charge
    - Peptide retention time
    - Peptide intensity
    - Sample accession
The feature file is a column format that defines the peptide quantification/identification and its relation with each
sample in the experiment.
"""
import random

import pyarrow as pa
import pyarrow.parquet as pq
import pandas as pd
from scipy.linalg._solve_toeplitz import float64

from quantms_io.core.mztab import MztabHandler
from quantms_io.core.parquet_handler import ParquetHandler
from quantms_io.core.sdrf import SDRFHandler
from quantms_io.utils.pride_utils import clean_peptidoform_sequence, compare_protein_lists, \
    standardize_protein_list_accession, standardize_protein_string_accession, get_modifications_object_from_mztab_line


def get_quantmsio_modifications(modifications_string: str, modification_definition: dict) -> dict:
    """
    Get the modifications in quantms.io format from a string of modifications in mztab format.
    :param modifications_string: modifications string in mztab format
    :return: modifications in quantms.io format
    """
    if modifications_string is None or modifications_string == "null":
        return {}

    return get_modifications_object_from_mztab_line(modification_string = modifications_string,
                                                             modifications_definition=modification_definition)


class FeatureHandler(ParquetHandler):
    """
    This class handle protein tables in column format. The main serialization format is Apache Parquet.
    """

    FEATURE_FIELDS = [pa.field("sequence", pa.string(),
                               metadata={"description": "Peptide sequence of the feature"}),
                      pa.field("protein_accessions", pa.list_(pa.string()),
                               metadata={"description": "accessions of associated proteins"}),
                      pa.field("protein_start_positions", pa.list_(pa.int32()),
                               metadata={"description": "start positions in the associated proteins"}),
                      pa.field("protein_end_positions", pa.list_(pa.int32()),
                               metadata={"description": "end positions in the associated proteins"}),
                      pa.field("protein_global_qvalue", pa.float64(),
                               metadata={"description": "global q-value of the associated protein or protein group"}),
                      pa.field("protein_best_id_score", pa.string(),
                               metadata={"description": "best identification score of the associated protein or protein group"}),
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
                               metadata={"description": "flag indicating if the feature is a decoy (1 is decoy, 0 is not decoy)"}),
                      pa.field("best_id_score", pa.string(),
                               metadata={"description": "best identification score as key value pair"}),
                      pa.field("intensity", pa.float64(),
                               metadata={"description": "intensity value"}),
                      pa.field("spectral_count", pa.int32(),
                               metadata={"description": "number of spectral counts"}),
                      pa.field("sample_accession", pa.string(),
                               metadata={"description": "accession of the associated sample"}),
                      pa.field("condition", pa.string(),
                               metadata={"description": "experimental condition, value of the experimental factor"}),
                      pa.field("fraction", pa.string(),
                               metadata={"description": "fraction information"}),
                      pa.field("biological_replicate", pa.string(),
                               metadata={"description": "biological replicate information"}),
                      pa.field("fragment_ion", pa.string(),
                               metadata={"description": "fragment ion information"}),
                      pa.field("isotope_label_type", pa.string(),
                               metadata={"description": "type of isotope label"}),
                      pa.field("run", pa.string(),
                               metadata={"description": "experimental run information"}),
                      pa.field("channel", pa.string(),
                               metadata={"description": "experimental channel information"}),
                      pa.field("id_scores", pa.list_(pa.string()),
                               metadata={"description": "identification scores as key value pairs"}),
                      pa.field("consensus_support", pa.float64(),
                               metadata={"description": "consensus support value"}),
                      pa.field("reference_file_name", pa.string(),
                               metadata={"description": "file name of the reference file"}),
                      pa.field("scan_number", pa.string(),
                               metadata={"description": "scan number of the best PSM"}),
                      pa.field("mz", pa.list_(pa.float64()),
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

    def __init__(self, parquet_path: str = None):
        self.schema = self._create_schema()
        self.parquet_path = parquet_path
        self.dataset = None

    def _create_schema(self):
        """
        Create the schema for the protein file. The schema is defined in the docs folder of this repository.
        (https://github.com/bigbio/quantms.io/blob/main/docs/FEATURE.md)
        """
        return pa.schema(FeatureHandler.FEATURE_FIELDS, metadata={"description": "Feature file in quantms.io format"})

    def read_feature_table(self) -> pa.Table:
        table = pq.ParquetDataset(
            self.parquet_path, use_legacy_dataset=False, schema=self.schema).read()  # type: pa.Table
        return table

    def create_feature_table(self, feature_list: list):
        return pa.Table.from_pandas(pd.DataFrame(feature_list), schema=self.schema)

    def convert_mztab_msstats_to_feature(self, msstats_file: str, sdrf_file: str, mztab_file: str,
                                         use_cache: bool = False):
        """
        Convert a MSstats input file and mztab into a quantms.io file format.
        :param msstats_file: MSstats input file
        :param sdrf_file: SDRF file
        :param mztab_file: mztab file
        :param use_cache: use cache to store the mztab file
        :return: none
        """
        sdrf_handler = SDRFHandler(sdrf_file)
        mztab_handler = MztabHandler(mztab_file, use_cache=use_cache)
        mztab_handler.load_mztab_file(use_cache=use_cache)
        feature_list = []

        # Read the MSstats file line by line and convert it to a quantms.io file feature format
        with open(msstats_file, "r") as msstats_file_handler:
            line = msstats_file_handler.readline()
            if line.startswith("Protein"):
                # Skip the header
                msstats_columns = line.rstrip().split(",")
            else:
                raise Exception("The MSstats file does not have the expected header")
            line = msstats_file_handler.readline()
            while line.rstrip() != "":
                line = line.rstrip()
                msstats_values = line.split(",")
                # Create a dictionary with the values
                feature_dict = dict(zip(msstats_columns, msstats_values))
                msstats_feature = self._fetch_msstats_feature(feature_dict, sdrf_handler, mztab_handler)
                if msstats_feature is not None:
                    feature_list.append(msstats_feature)
                #feature_table = self.create_feature_table([msstats_feature])
                print(msstats_feature)
                line = msstats_file_handler.readline()
                # Create a feature table
        feature_table = self.create_feature_table(feature_list)
        # Write the feature table to a parquet file

        self.parquet_path = "feature.parquet"
        self.write_single_file_parquet(feature_table, write_metadata=True)
        mztab_handler.close()


    def describe_schema(self):
        schema_description = []
        for field in self.schema:
            field_description = {
                "name": field.name,
                "type": str(field.type),
                "description": field.metadata.get("description", "")
            }
            schema_description.append(field_description)
        return schema_description

    def _fetch_msstats_feature(self, feature_dict: dict, sdrf_handler: SDRFHandler, mztab_handler: MztabHandler):
        """
        Fetch a feature from a MSstats dictionary and convert to a quantms.io format row
        :param feature_dict: MSstats feature dictionary
        :param sdrf_handler: SDRF handler
        :param mztab_handler: mztab handler
        :return: quantms.io format row
        """
        protein_accessions_list = standardize_protein_list_accession(feature_dict["ProteinName"])
        protein_accessions_string = standardize_protein_string_accession(feature_dict["ProteinName"])

        peptidoform = feature_dict["PeptideSequence"] # Peptidoform is the Msstats form .e.g. EM(Oxidation)QDLGGGER
        peptide_sequence = clean_peptidoform_sequence(peptidoform)  # Get sequence .e.g. EMQDLGGGER
        charge = feature_dict["PrecursorCharge"] # Charge is the Msstats form .e.g. 2
        reference_file = feature_dict["Reference"].replace("\"", "").split(".")[0] # Reference is the Msstats form .e.g. HeLa_1to1_01

        peptide_mztab_qvalue = mztab_handler.get_peptide_qvalue_from_index(msstats_peptidoform=peptidoform,
                                                                           charge=charge)
        peptide_score_name = mztab_handler.get_search_engine_scores()["peptide_score"]

        peptide_mztab_qvalue_accession = standardize_protein_list_accession(peptide_mztab_qvalue[0])
        peptide_qvalue = peptide_mztab_qvalue[1] # Peptide q-value index 1
        peptide_mztab_qvalue_decoy = peptide_mztab_qvalue[2] # Peptide q-value decoy index 2

        # Mods in quantms.io format
        modifications_string = peptide_mztab_qvalue[4] # Mods
        modifications = get_quantmsio_modifications(modifications_string= modifications_string,
                                                    modification_definition=mztab_handler.get_modifications_definition())
        modifications_string = ""
        for key, value in modifications.items():
            modifications_string += "|".join(map(str, value["position"]))
            modifications_string = modifications_string + "-" + value["unimod_accession"] + ","
        modifications_string = None if len(modifications_string)==0 else modifications_string[:-1] # Remove last comma
        try:
            peptide_count = mztab_handler.get_number_psm_from_index(msstats_peptidoform=peptidoform,
                                                                charge=charge, reference_file=reference_file)
        except:
            print(f"MBR peptide: {peptidoform}, {charge}, {reference_file}")
            return None

        start_positions = peptide_count[2].split(",") # Start positions in the protein
        start_positions = [int(i) for i in start_positions]

        end_positions = peptide_count[3].split(",")   # End positions in the protein
        end_positions = [int(i) for i in end_positions]
        peptide_count_accession = standardize_protein_list_accession(peptide_count[1]) # Accession of the protein
        spectral_count = peptide_count[0] # Spectral count

        # find calculated mass and experimental mass in peptide_count. Here we are using the scan number and the
        # reference file name to find the calculated mass and the experimental mass.
        scan_number = peptide_mztab_qvalue[6]
        mztab_reference_file = peptide_mztab_qvalue[5]
        calculated_mass = None
        experimental_mass = None
        for row in peptide_count[4]:
            if row[3] == scan_number and row[2] == mztab_reference_file:
                calculated_mass = row[0]
                experimental_mass = row[1]
                break

        try:
            protein_qvalue_object = mztab_handler.get_protein_qvalue_from_index(protein_accession=protein_accessions_string)
            protein_qvalue = protein_qvalue_object[0] # Protein q-value index 0
        except:
            print("Error in line: {}".format(feature_dict))
            return None

        # TODO: Importantly, the protein accessions in the peptide section of the mzTab files peptide_mztab_qvalue_accession
        #  are not the same as the protein accessions in the protein section of the mzTab file protein_accessions_list, that is
        #  why we are using the protein_accessions_list to compare with the protein accessions in the MSstats file.
        if not compare_protein_lists(protein_accessions_list, peptide_count_accession):
            print("Error in line: {}".format(feature_dict))
            print("Protein accessions: {}-{}-{}".format(protein_accessions_list, peptide_mztab_qvalue_accession, peptide_count_accession))
            raise Exception("The protein accessions in the MSstats file do not match with the mztab file")

        # Unique is different in PSM section, Peptide, We are using the msstats number of accessions.
        unique = 1 if len(protein_accessions_list) == 1 else 0

        # TODO: get retention time from mztab file. The rentention time in the mzTab peptide level is not clear how is
        # related to the retention time in the MSstats file. The retention time in the MSstats file is the retention time.
        rt = None
        return {
            "sequence": peptide_sequence,
            "protein_accessions":protein_accessions_list,
            "protein_start_positions": start_positions,
            "protein_end_positions": end_positions,
            "protein_global_qvalue": float64(protein_qvalue),
            "protein_best_id_score": None,
            "unique": unique,
            "modifications": modifications_string,
            "retention_time": rt,
            "charge": int(charge),
            "calc_mass_to_charge": float64(calculated_mass),
            "peptidoform": peptide_mztab_qvalue[0],
            "posterior_error_probability": None,
            "global_qvalue": float64(peptide_qvalue),
            "is_decoy": int(peptide_mztab_qvalue_decoy),
            "best_id_score": f"{peptide_score_name}: {peptide_qvalue}",
            "intensity": float64(feature_dict["Intensity"]),
            "spectral_count": spectral_count,
            "sample_accession": "S1",
            "condition": "ConditionA",
            "fraction": "F1",
            "biological_replicate": "BioRep1",
            "fragment_ion": "b2",
            "isotope_label_type": "LabelA",
            "run": "Run1",
            "channel": "ChannelA",
            "id_scores": [f"{peptide_score_name}: {peptide_mztab_qvalue}"],
            "consensus_support": None,
            "reference_file_name": reference_file,
            "scan_number": scan_number,
            "exp_mass_to_charge": float64(experimental_mass),
            "mz": None,
            "intensity_array": None,
            "num_peaks": None,
            "gene_accessions": None,
            "gene_names": None
        }



