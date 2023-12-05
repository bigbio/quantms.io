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
import numpy as np
import pandas as pd

import pyarrow as pa
import pyarrow.parquet as pq

from quantms_io.core.feature_in_memory import FeatureInMemory
from quantms_io.core.mztab import MztabHandler
from quantms_io.core.openms import OpenMSHandler
from quantms_io.core.parquet_handler import ParquetHandler
from quantms_io.core.sdrf import SDRFHandler
from quantms_io.utils.constants import ITRAQ_CHANNEL, TMT_CHANNELS
from quantms_io.utils.pride_utils import (clean_peptidoform_sequence,
                                          compare_protein_lists,
                                          get_quantmsio_modifications,
                                          standardize_protein_list_accession,
                                          standardize_protein_string_accession)

import logging
logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__)

def get_msstats_in_batches(msstats_file: str, batch_size: int) -> int:
    """
    :param msstats_file: MSstats input file
    :param batch_size: batch size
    :return: number of batches
    """
    total_len = sum(1 for _ in open(msstats_file)) - 1  # neglect header
    if total_len <= batch_size:
        return 1
    elif total_len % batch_size != 0:
        return total_len // batch_size + 1
    else:
        return total_len // batch_size


def get_additional_properties_from_sdrf(
    feature_dict: dict, experiment_type: str, sdrf_samples: dict
) -> dict:
    if (
        "FragmentIon" not in feature_dict
    ):  # FragmentIon is not in the MSstats if the experiment is label free
        feature_dict["FragmentIon"] = None

    if "IsotopeLabelType" not in feature_dict:
        feature_dict["IsotopeLabelType"] = "L"

    if "TMT" in experiment_type:  # TMT experiment
        feature_dict["Channel"] = TMT_CHANNELS[experiment_type][
            int(feature_dict["Channel"]) - 1
        ]
    elif "ITRAQ" in experiment_type:  # ITRAQ experiment
        feature_dict["Channel"] = ITRAQ_CHANNEL[experiment_type][
            int(feature_dict["Channel"]) - 1
        ]
    else:  # Label free experiment
        feature_dict["Channel"] = "LABEL FREE SAMPLE"

    # Get the sample accession from the SDRF file.
    sample_id = sdrf_samples[
        feature_dict["Reference"] + ":_:" + feature_dict["Channel"]
    ]
    feature_dict["SampleName"] = sample_id

    return feature_dict


def _fetch_msstats_feature(
    feature_dict: dict,
    experiment_type: str,
    sdrf_samples: dict,
    mztab_handler: MztabHandler = None,
    intensity_map: dict = None,
):
    """
    fetch a feature from a MSstats dictionary and convert to a quantms.io format row
    :param feature_dict: MSstats feature dictionary
    :param experiment_type: experiment type
    :param sdrf_samples: SDRF samples
    :param mztab_handler: mztab handler
    :param intensity_map: intensity map
    :return: quantms.io format row
    """

    # Reference is the Msstats form .e.g. HeLa_1to1_01
    feature_dict["Reference"] = feature_dict["Reference"].replace('"', "").split(".")[0]
    feature_dict = get_additional_properties_from_sdrf(
        feature_dict, experiment_type, sdrf_samples
    )

    protein_accessions_list = standardize_protein_list_accession(
        feature_dict["ProteinName"]
    )
    protein_accessions_string = standardize_protein_string_accession(
        feature_dict["ProteinName"]
    )

    peptidoform = feature_dict[
        "PeptideSequence"
    ]  # Peptidoform is the Msstats form .e.g. EM(Oxidation)QDLGGGER
    peptide_sequence = clean_peptidoform_sequence(
        peptidoform
    )  # Get sequence .e.g. EMQDLGGGER
    charge = None
    if "PrecursorCharge" in feature_dict:
        charge = feature_dict["PrecursorCharge"]  # Charge is the Msstats form .e.g. 2
    elif "Charge" in feature_dict:
        charge = feature_dict["Charge"]  # Charge is the Msstats form .e.g. 2

    peptide_indexed = mztab_handler.get_peptide_index(
        msstats_peptidoform=peptidoform, charge=charge
    )
    peptide_score_name = mztab_handler.get_search_engine_scores()["peptide_score"]

    peptide_mztab_qvalue_accession = standardize_protein_list_accession(
        peptide_indexed["protein_accession"]
    )
    peptide_qvalue = peptide_indexed["peptide_qvalue"]  # Peptide q-value index 1

    posterior_error_probability = peptide_indexed[
        "posterior_error_probability"
    ]  # Get posterior error probability
    if posterior_error_probability is not None:
        posterior_error_probability = np.float64(posterior_error_probability)

    # Get if decoy or not
    peptide_mztab_qvalue_decoy = peptide_indexed[
        "is_decoy"
    ]  # Peptide q-value decoy index 2

    # Mods in quantms.io format
    modifications_string = peptide_indexed["modifications"]  # Mods
    modifications = get_quantmsio_modifications(
        modifications_string=modifications_string,
        modification_definition=mztab_handler.get_modifications_definition(),
    )
    modifications_string = ""
    for key, value in modifications.items():
        modifications_string += "|".join(map(str, value["position"]))
        modifications_string = (
            modifications_string + "-" + value["unimod_accession"] + ","
        )
    modifications_string = (
        None if len(modifications_string) == 0 else modifications_string[:-1]
    )  # Remove last comma
    modification_list = (
        None if modifications_string is None else modifications_string.split(",")
    )

    start_positions = peptide_indexed["psm_protein_start"].split(
        ","
    )  # Start positions in the protein
    start_positions = [int(i) for i in start_positions]
    end_positions = peptide_indexed["psm_protein_end"].split(
        ","
    )  # End positions in the protein
    end_positions = [int(i) for i in end_positions]

    # The spectral count is cero for a given file if a peptide does not appear in the psm section because for example is
    # match between runs. However, for TMT/ITRAQ experiments, the spectral count can't be provided.
    spectral_count = None if "LABEL FREE" not in feature_dict["Channel"] else 0  #

    spectral_count_list = peptide_indexed["spectral_count"]  # Spectral count
    if feature_dict["Reference"] in spectral_count_list:
        spectral_count = spectral_count_list[feature_dict["Reference"]]
    else:
        reference_file = feature_dict["Reference"]
        logger.debug(f"MBR peptide: {peptidoform}, {charge}, {reference_file}")

    # Find calculated mass and experimental mass in peptide_count. Here we are using the scan number and the
    # reference file name to find the calculated mass and the experimental mass.
    peptide_ms_run = peptide_indexed["peptide_ms_run"]
    peptide_scan_number = peptide_indexed["peptide_scan_number"]
    calculated_mass = None
    experimental_mass = None
    for row in peptide_indexed["list_psms"]:
        calculated_mass = row[0]
        if row[3] == peptide_scan_number and row[2] == peptide_ms_run:
            experimental_mass = row[1]
            break

    # If the experimental mass is not found, in the psm list this means that the peptide is MBR product.
    if experimental_mass is None:
        peptide_ms_run = None
        peptide_scan_number = None

    try:
        protein_qvalue_object = mztab_handler.get_protein_qvalue_from_index(
            protein_accession=protein_accessions_string
        )
        protein_qvalue = protein_qvalue_object[0]  # Protein q-value index 0
    except:
        logger.error("Error in line: {}".format(feature_dict))
        return None

    # TODO: Importantly, the protein accessions in the peptide section of the mzTab files peptide_mztab_qvalue_accession
    #  are not the same as the protein accessions in the protein section of the mzTab file protein_accessions_list, that is
    #  why we are using the protein_accessions_list to compare with the protein accessions in the MSstats file.
    if not compare_protein_lists(
        protein_accessions_list, peptide_mztab_qvalue_accession
    ):
        logger.error("Error in line: {}".format(feature_dict))
        logger.error(
            "Protein accessions: {}-{}".format(
                protein_accessions_list, peptide_mztab_qvalue_accession
            )
        )
        raise Exception(
            "The protein accessions in the MSstats file do not match with the mztab file"
        )

    # Unique is different in PSM section, Peptide, We are using the msstats number of accessions.
    unique = 1 if len(protein_accessions_list) == 1 else 0

    # TODO: get retention time from mztab file. The retention time in the mzTab peptide level is not clear how is
    # related to the retention time in the MSstats file. If the consensusxml file is provided, we can get the retention
    # time from the consensusxml file.
    rt = (
        None
        if "RetentionTime" not in feature_dict
        else np.float64(feature_dict["RetentionTime"])
    )
    if intensity_map is not None:
        key = None
        if "LABEL FREE" in feature_dict["Channel"]:
            key = peptidoform + ":_:" + str(charge) + ":_:" + feature_dict["Reference"]
        elif "TMT" in feature_dict["Channel"] or "ITRAQ" in feature_dict["Channel"]:
            key = (
                peptidoform
                + ":_:"
                + str(charge)
                + ":_:"
                + feature_dict["Reference"]
                + ":_:"
                + feature_dict["Channel"]
            )
        if key is not None and key in intensity_map:
            consensus_intensity = intensity_map[key]["intensity"]
            if abs(consensus_intensity - np.float64(feature_dict["Intensity"])) < 0.1:
                rt = rt if rt is not None else intensity_map[key]["rt"]
                experimental_mass = intensity_map[key]["mz"]

    return {
        "sequence": peptide_sequence,
        "protein_accessions": protein_accessions_list,
        "protein_start_positions": start_positions,
        "protein_end_positions": end_positions,
        "protein_global_qvalue": np.float64(protein_qvalue),
        "unique": unique,
        "modifications": modification_list,
        "retention_time": rt,
        "charge": int(charge),
        "calc_mass_to_charge": np.float64(calculated_mass),
        "peptidoform": peptide_indexed[
            "peptidoform"
        ],  # Peptidoform in proforma notation
        "posterior_error_probability": posterior_error_probability,
        "global_qvalue": np.float64(peptide_qvalue),
        "is_decoy": int(peptide_mztab_qvalue_decoy),
        "intensity": np.float64(feature_dict["Intensity"]),
        "spectral_count": spectral_count,
        "sample_accession": feature_dict["SampleName"],
        "condition": feature_dict["Condition"],
        "fraction": feature_dict["Fraction"],
        "biological_replicate": feature_dict["BioReplicate"],
        "fragment_ion": feature_dict["FragmentIon"],
        "isotope_label_type": feature_dict["IsotopeLabelType"],
        "run": feature_dict["Run"],
        "channel": feature_dict["Channel"],
        "id_scores": [
            f"{peptide_score_name}: {peptide_qvalue}",
            f"Best PSM PEP: {posterior_error_probability}",
        ],
        "reference_file_name": feature_dict["Reference"],
        "best_psm_reference_file_name": peptide_ms_run,
        "best_psm_scan_number": peptide_scan_number,
        "exp_mass_to_charge": np.float64(experimental_mass),
        "mz_array": None,
        "intensity_array": None,
        "num_peaks": None,
        "gene_accessions": None,
        "gene_names": None,
    }


class FeatureHandler(ParquetHandler):
    """
    this class handle protein tables in column format.
    the main serialization format is Apache Parquet.
    """

    FEATURE_FIELDS = [
        pa.field(
            "sequence",
            pa.string(),
            metadata={"description": "Peptide sequence of the feature"},
        ),
        pa.field(
            "protein_accessions",
            pa.list_(pa.string()),
            metadata={"description": "accessions of associated proteins"},
        ),
        pa.field(
            "protein_start_positions",
            pa.list_(pa.int32()),
            metadata={"description": "start positions in the associated proteins"},
        ),
        pa.field(
            "protein_end_positions",
            pa.list_(pa.int32()),
            metadata={"description": "end positions in the associated proteins"},
        ),
        pa.field(
            "protein_global_qvalue",
            pa.float64(),
            metadata={
                "description": "global q-value of the associated protein or protein group"
            },
        ),
        pa.field(
            "unique",
            pa.int32(),
            metadata={
                "description": "if the peptide is unique to a particular protein"
            },
        ),
        pa.field(
            "modifications",
            pa.list_(pa.string()),
            metadata={"description": "peptide modifications"},
        ),
        pa.field(
            "retention_time", pa.float64(), metadata={"description": "retention time"}
        ),
        pa.field(
            "charge",
            pa.int32(),
            metadata={"description": "charge state of the feature"},
        ),
        pa.field(
            "exp_mass_to_charge",
            pa.float64(),
            metadata={"description": "experimentally measured mass-to-charge ratio"},
        ),
        pa.field(
            "calc_mass_to_charge",
            pa.float64(),
            metadata={"description": "calculated mass-to-charge ratio"},
        ),
        pa.field(
            "peptidoform",
            pa.string(),
            metadata={"description": "peptidoform in proforma notation"},
        ),
        pa.field(
            "posterior_error_probability",
            pa.float64(),
            metadata={"description": "posterior error probability"},
        ),
        pa.field(
            "global_qvalue", pa.float64(), metadata={"description": "global q-value"}
        ),
        pa.field(
            "is_decoy",
            pa.int32(),
            metadata={
                "description": "flag indicating if the feature is a decoy (1 is decoy, 0 is not decoy)"
            },
        ),
        # pa.field("best_id_score", pa.string(),
        #          metadata={"description": "best identification score as key value pair"}),
        pa.field(
            "intensity", pa.float64(), metadata={"description": "intensity value"}
        ),
        pa.field(
            "spectral_count",
            pa.int32(),
            metadata={"description": "number of spectral counts"},
        ),
        pa.field(
            "sample_accession",
            pa.string(),
            metadata={"description": "accession of the associated sample"},
        ),
        pa.field(
            "condition",
            pa.string(),
            metadata={
                "description": "experimental condition, value of the experimental factor"
            },
        ),
        pa.field(
            "fraction", pa.string(), metadata={"description": "fraction information"}
        ),
        pa.field(
            "biological_replicate",
            pa.string(),
            metadata={"description": "biological replicate information"},
        ),
        pa.field(
            "fragment_ion",
            pa.string(),
            metadata={"description": "fragment ion information"},
        ),
        pa.field(
            "isotope_label_type",
            pa.string(),
            metadata={"description": "type of isotope label"},
        ),
        pa.field(
            "run", pa.string(), metadata={"description": "experimental run information"}
        ),
        pa.field(
            "channel",
            pa.string(),
            metadata={"description": "experimental channel information"},
        ),
        pa.field(
            "id_scores",
            pa.list_(pa.string()),
            metadata={"description": "identification scores as key value pairs"},
        ),
        # pa.field("consensus_support", pa.float64(),
        #          metadata={"description": "consensus support value"}),
        pa.field(
            "reference_file_name",
            pa.string(),
            metadata={"description": "file name of the reference file for the feature"},
        ),
        pa.field(
            "best_psm_reference_file_name",
            pa.string(),
            metadata={
                "description": "file name of the reference file for the best PSM"
            },
        ),
        pa.field(
            "best_psm_scan_number",
            pa.string(),
            metadata={"description": "scan number of the best PSM"},
        ),
        pa.field(
            "mz_array",
            pa.list_(pa.float64()),
            metadata={"description": "mass-to-charge ratio values"},
        ),
        pa.field(
            "intensity_array",
            pa.list_(pa.float64()),
            metadata={"description": "intensity array values"},
        ),
        pa.field("num_peaks", pa.int32(), metadata={"description": "number of peaks"}),
        pa.field(
            "gene_accessions",
            pa.list_(pa.string()),
            metadata={"description": "accessions of associated genes"},
        ),
        pa.field(
            "gene_names",
            pa.list_(pa.string()),
            metadata={"description": "names of associated genes"},
        ),
    ]

    def __init__(self, parquet_path: str = None):
        super().__init__(parquet_path)
        self.schema = self._create_schema()
        self.parquet_path = parquet_path
        self.dataset = None

    def _create_schema(self):
        """
        Create the schema for the feature file. The schema is defined in the docs folder of this repository.
        (https://github.com/bigbio/quantms.io/blob/main/docs/FEATURE.md)
        """
        return pa.schema(
            FeatureHandler.FEATURE_FIELDS,
            metadata={"description": "Feature file in quantms.io format"},
        )

    def read_feature_table(self) -> pa.Table:
        table = pq.ParquetDataset(
            self.parquet_path, use_legacy_dataset=False, schema=self.schema
        ).read()  # type: pa.Table
        return table

    def create_feature_table(self, feature_list: list) -> pa.Table:
        return pa.Table.from_pandas(pd.DataFrame(feature_list), schema=self.schema)

    def convert_mztab_msstats_to_feature(
        self,
        msstats_file: str,
        sdrf_file: str,
        mztab_file: str,
        consesusxml_file: str = None,
        batch_size: int = 1000000,
        use_cache: bool = False,
    ):
        """
        convert a MSstats input file and mztab into a quantms.io file format.
        :param msstats_file: MSstats input file
        :param sdrf_file: SDRF file
        :param mztab_file: mztab file
        :param consesusxml_file: consensusxml file
        :param batch_size: batch size
        :param use_cache: use cache to store the mztab file
        :return: none
        """
        sdrf_handler = SDRFHandler(sdrf_file)
        experiment_type = sdrf_handler.get_experiment_type_from_sdrf()
        # Get the intensity map from the consensusxml file
        intensity_map = {}
        if consesusxml_file is not None:
            consensus_handler = OpenMSHandler()
            intensity_map = consensus_handler.get_intensity_map(
                    consensusxml_path=consesusxml_file, experiment_type=experiment_type
            )
        if use_cache:
            sdrf_samples = sdrf_handler.get_sample_map()

            mztab_handler = MztabHandler(mztab_file, use_cache=use_cache)
            mztab_handler.load_mztab_file(use_cache=use_cache)

            feature_list = []

            # Read the MSstats file line by line and convert it to a quantms.io file feature format
            batches = get_msstats_in_batches(msstats_file, batch_size)
            batch_count = 1
            with open(msstats_file, "r") as msstats_file_handler:
                line = msstats_file_handler.readline()
                if "Protein" in line:
                    # Skip the header
                    msstats_columns = line.rstrip().split(",")
                else:
                    raise Exception(
                        "The MSstats file does not have the expected header"
                    )
                line = msstats_file_handler.readline()
                pqwriter = None
                while line.rstrip() != "":
                    line = line.rstrip()
                    msstats_values = line.split(",")
                    # Create a dictionary with the values
                    feature_dict = dict(zip(msstats_columns, msstats_values))
                    msstats_feature = _fetch_msstats_feature(
                        feature_dict,
                        experiment_type,
                        sdrf_samples,
                        mztab_handler,
                        intensity_map,
                    )
                    if msstats_feature is not None:
                        feature_list.append(msstats_feature)
                    # batches > 1
                    if len(feature_list) == batch_size and batch_count < batches:
                        feature_table = self.create_feature_table(feature_list)
                        feature_list = []
                        batch_count += 1
                        if not pqwriter:
                            pqwriter = pq.ParquetWriter(
                                self.parquet_path, feature_table.schema
                            )
                        pqwriter.write_table(feature_table)
                    line = msstats_file_handler.readline()
                # batches = 1
                if batch_count == 1:
                    feature_table = self.create_feature_table(feature_list)
                    pqwriter = pq.ParquetWriter(self.parquet_path, feature_table.schema)
                    pqwriter.write_table(feature_table)
                # final batch
                else:
                    feature_table = self.create_feature_table(feature_list)
                    pqwriter.write_table(feature_table)
                if pqwriter:
                    pqwriter.close()
            mztab_handler.close()
        else:
            convert = FeatureInMemory(experiment_type, self.schema)
            convert.merge_mztab_and_sdrf_to_msstats_in(
                mztab_file,
                msstats_file,
                sdrf_file,
                self.parquet_path,
                msstats_chunksize=batch_size,
            )

    def describe_schema(self):
        """
        describe the schema of the feature file
        """
        schema_description = []
        for field in self.schema:
            field_description = {
                "name": field.name,
                "type": str(field.type),
                "description": field.metadata.get("description", ""),
            }
            schema_description.append(field_description)
        return schema_description
