import numpy as np

from quantms_io.core.core import DiskCache
from quantms_io.utils.constants import PROTEIN_DETAILS
from quantms_io.utils.file_utils import calculate_buffer_size
from quantms_io.utils.pride_utils import (fetch_modifications_from_mztab_line,
                                          fetch_ms_runs_from_mztab_line,
                                          fetch_peptide_from_mztab_line,
                                          fetch_protein_from_mztab_line,
                                          fetch_psm_from_mztab_line,
                                          get_key_peptide_combination,
                                          get_permutations_of_original_list,
                                          parse_score_name_in_mztab,
                                          standardize_protein_string_accession)

import logging
logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__)


def create_peptide_for_index(
    peptide_qvalue: str,
    is_decoy: str,
    peptidoform: str,
    modifications: str,
    peptide_ms_run: str,
    peptide_scan_number: str,
    peptide_position_mztab: str,
    psm_protein_start: str = None,
    psm_protein_end: str = None,
    posterior_error_probability: str = None,
    spectral_count: dict = None,
    list_psms: list = None,
) -> dict:
    """
    Create a peptide to be store in the index. An index peptide have the following structure:
      # - protein_accession: protein accession
      # - peptide_qvalue: peptide qvalue
      # - is_decoy: is decoy
      # - peptidoform: peptidoform in proforma notation
      # - modifications:  modifications
      # - posterior_error_probability: posterior error probability
      # - peptide_ms_run: ms run
      # - peptide_scan_number: scan number
      # - peptide_position_mztab: position in the mztab file
      # - psm_protein_start: protein start
      # - psm_protein_end: protein end
      # - spectral_count:
           # - reference file: spectral count
      # - list_psms:
           # - 0 - calculated mass to charge
           # - 1 - experimental mass to charge
           # - 2 - reference file
           # - 3 - scan number
           # - 4 - position in the mztab file
     :param peptide_qvalue: peptide qvalue
     :param is_decoy: is decoy
     :param peptidoform: peptidoform in proforma
     :param modifications: modifications
     :param posterior_error_probability: posterior error probability
     :param peptide_ms_run: ms run
     :param peptide_scan_number: scan number
     :param peptide_position_mztab: position in the mztab file
     :param psm_protein_start: protein start
     :param psm_protein_end: protein end
     :param spectral_count: spectral count
     :param list_psms: list of psms
     :return: peptide dictionary
    """
    if spectral_count is None:
        spectral_count = {}
    if list_psms is None:
        list_psms = []
    peptide = {
        "peptide_qvalue": peptide_qvalue,
        "is_decoy": is_decoy,
        "peptidoform": peptidoform,
        "modifications": modifications,
        "peptide_ms_run": peptide_ms_run,
        "peptide_scan_number": peptide_scan_number,
        "peptide_position_mztab": peptide_position_mztab,
        "psm_protein_start": psm_protein_start,
        "psm_protein_end": psm_protein_end,
        "spectral_count": spectral_count,
        "list_psms": list_psms,
        "posterior_error_probability": posterior_error_probability,
    }
    return peptide


class MztabHandler:
    """
    mztab handler class for quantms.io. this class is used to load mztab files and create indexes for the mztab files.
    the mztab handler has multiple indexes to access the mztab content: protein, peptide and psm.
    - Protein index: contains the position and the qvalue of each protein in the mztab file.
    - Peptide index: contains the position and the qvalue/modifications of each peptide in the mztab file.
    """

    def __init__(self, mztab_file, use_cache: bool = False):
        self._search_engine_scores = {}  # search engine scores names

        # The key is the combination of the msstats_peptidoform notations and the charge state. The value is a list
        # with the following information:
        # - protein_accession: protein accession
        # - peptide_qvalue: peptide qvalue
        # - is_decoy: is decoy
        # - peptidoform: peptidoform in proforma notation
        # - modifications:  modifications
        # - peptide_ms_run: ms run
        # - peptide_scan_number: scan number
        # - peptide_position_mztab: position in the mztab file
        # - psm_protein_start: protein start
        # - psm_protein_end: protein end
        # - spectral_count:
        # - reference file: spectral count
        # - list_psms:
        #  - 0 - calculated mass to charge
        #  - 1 - experimental mass to charge
        #  - 2 - reference file
        #  - 3 - scan number
        #  - 4 - position in the mztab file
        self._peptide_index = None  # peptide index details

        # The dictionary contains a key value pair of the information from the protein, the key is the accessions of the
        # proteins, the value is the qvalue score and the position in the mztab file.
        self._protein_details = None  # protein index details

        self._ms_runs = {}  # ms runs index
        self._modifications = {}
        self._mztab_file = mztab_file  # mztab file
        self._use_cache = use_cache  # use cache to store peptide information

    def create_peptide_index(self):
        """
        A peptide index is a dictionary where the key is the combination of the msstats_peptidoform notations and the
        the charge state. The value can be read in self._peptide_index.
        """
        if self._use_cache:
            self._peptide_index = DiskCache("peptide_details")
        else:
            self._peptide_index = {}

    def add_peptide_to_index(self, peptide_key: str, peptide_value: dict):
        """
        add a peptide to the index if the file is big then use cache.
        :param peptide_key: peptide key
        :param peptide_value: peptide value
        :return: None
        """
        if self._use_cache:
            if self._peptide_index.contains(peptide_key):
                pep = self._peptide_index.get_item(peptide_key)
                if float(peptide_value['peptide_qvalue']) < float(pep['peptide_qvalue']):
                    self._peptide_index.add_item(peptide_key, peptide_value)
            else:
                self._peptide_index.add_item(peptide_key, peptide_value)
        else:
            self._peptide_index[peptide_key] = peptide_value

    def add_psm_to_count(
        self,
        psm_key: str,
        protein_accession: str,
        protein_start: str,
        protein_end: str,
        calc_mass_to_charge: str,
        exp_mass_to_charge: str,
        reference_file: str,
        scan_number: str,
        mztab_postion: str,
        posterior_error_probability: str = None,
    ):
        """
        psm count is used to count the number of psms for a given peptidoform in a file.
        :param psm_key: psm key
        :param protein_accession: protein accession
        :param protein_start: protein start
        :param protein_end: protein end
        :param calc_mass_to_charge: calculated mass to charge
        :param exp_mass_to_charge: experimental mass to charge
        :param reference_file: reference file
        :param scan_number: scan number
        :param mztab_postion: position in the mztab file
        :param posterior_error_probability: posterior error probability
        :return: None
        """
        if self._use_cache:
            if self._peptide_index.contains(psm_key):
                psm = self._peptide_index.get_item(psm_key)
                peptide_protein_accession = (
                    psm["protein_accession"] if "protein_accession" in psm else None
                )
                peptide_protein_start = (
                    psm["psm_protein_start"] if "psm_protein_start" in psm else None
                )
                peptide_protein_end = (
                    psm["psm_protein_end"] if "psm_protein_end" in psm else None
                )
                if posterior_error_probability is not None:
                    if psm["posterior_error_probability"] is None or (
                        np.float64(posterior_error_probability)
                        < psm["posterior_error_probability"]
                    ):
                        psm["posterior_error_probability"] = np.float64(
                            posterior_error_probability
                        )

                if (
                    peptide_protein_start is not None
                    and peptide_protein_end is not None
                ) and [
                    peptide_protein_accession,
                    peptide_protein_start,
                    peptide_protein_end,
                ] != [
                    protein_accession,
                    protein_start,
                    protein_end,
                ]:
                    raise Exception(
                        "The protein information is different for the same psm key"
                    )
                else:
                    psm["psm_protein_start"] = protein_start
                    psm["psm_protein_end"] = protein_end
                    psm["protein_accession"] = protein_accession

                spectral_count = psm[
                    "spectral_count"
                ]  # dictionary with the spectral count reference file: spectral count
                if spectral_count is None:
                    spectral_count = {}
                if reference_file not in spectral_count:
                    spectral_count[reference_file] = 1
                else:
                    spectral_count[reference_file] = spectral_count[reference_file] + 1
                psm["spectral_count"] = spectral_count

                list_psms = psm["list_psms"]  # list of psms information
                if list_psms is None:
                    list_psms = []
                list_psms.append(
                    [
                        calc_mass_to_charge,
                        exp_mass_to_charge,
                        reference_file,
                        scan_number,
                        mztab_postion,
                    ]
                )
                psm["list_psms"] = list_psms

                self._peptide_index.add_item(psm_key, psm)
        else:
            if psm_key in self._peptide_index:
                psm = self._peptide_index[psm_key]
                peptide_protein_accession = psm["protein_accession"]
                peptide_protein_start = psm["psm_protein_start"]
                peptide_protein_end = psm["psm_protein_end"]
                if (
                    posterior_error_probability is not None
                    and posterior_error_probability < psm["posterior_error_probability"]
                ):
                    psm["posterior_error_probability"] = posterior_error_probability
                if (
                    peptide_protein_start is not None
                    and peptide_protein_end is not None
                ) and [
                    peptide_protein_accession,
                    peptide_protein_start,
                    peptide_protein_end,
                ] != [
                    protein_accession,
                    protein_start,
                    protein_end,
                ]:
                    raise Exception(
                        "The protein information is different for the same psm key"
                    )
                else:
                    psm["psm_protein_start"] = protein_start
                    psm["psm_protein_end"] = protein_end

                spectral_count = psm[
                    "spectral_count"
                ]  # dictionary with the spectral count reference file: spectral count
                if spectral_count is None:
                    spectral_count = {}
                if reference_file not in spectral_count:
                    spectral_count[reference_file] = 1
                else:
                    spectral_count[reference_file] = spectral_count[reference_file] + 1
                psm["spectral_count"] = spectral_count

                list_psms = psm["list_psms"]  # list of psms information
                if list_psms is None:
                    list_psms = []
                list_psms.append(
                    [
                        calc_mass_to_charge,
                        exp_mass_to_charge,
                        reference_file,
                        scan_number,
                    ]
                )
                psm["list_psms"] = list_psms

                self._peptide_index[psm_key] = psm

    def create_mztab_index(
        self, mztab_file: str, qvalue_index: bool = True, psm_count_index: bool = True
    ):
        """
        create an index for a mztab file; the index contains a structure with the position of each psm, peptide and
        protein in the file.
        :param mztab_file: mztab file
        :param qvalue_index: create a qvalue index
        :param psm_count_index: create a psm count index
        """
        self._protein_details = {}
        self.create_peptide_index()

        with open(mztab_file, "r", buffering=calculate_buffer_size(mztab_file)) as f:
            line = f.readline()
            while line != "":
                line = line.strip()
                if line.startswith("MTD") and "location" in line:
                    self._ms_runs = fetch_ms_runs_from_mztab_line(line, self._ms_runs)
                if line.startswith("MTD") and "search_engine_score" in line:
                    self._get_search_engine_scores(line)
                if line.startswith("MTD") and "_mod[":
                    self._modifications = fetch_modifications_from_mztab_line(
                        line, self._modifications
                    )
                elif line.startswith("PRH"):
                    logger.info("-- End of the Metadata section of the mztab file -- ")
                    protein_columns = line.split("\t")
                elif line.startswith("PRT"):
                    protein_info = line.split("\t")
                    if PROTEIN_DETAILS not in line:
                        es = dict(zip(protein_columns, protein_info))
                        protein = fetch_protein_from_mztab_line(es)
                        self._protein_details[protein["accession"]] = [
                            protein["score"]
                        ]
                elif line.startswith("PEH"):
                    logger.info("-- All proteins have been read, starting peptide section -- ")
                    peptide_columns = line.split("\t")
                elif line.startswith("PEP"):
                    peptide_info = line.split("\t")
                    es = dict(zip(peptide_columns, peptide_info))
                    peptide = fetch_peptide_from_mztab_line(
                        es,
                        ms_runs=self._ms_runs,
                        modification_definition=self._modifications,
                    )
                    peptide_key = get_key_peptide_combination(
                        msstats_peptidoform=peptide["peptidoform"],
                        charge=peptide["charge"],
                    )

                    peptide_value = create_peptide_for_index(
                        peptide_qvalue=peptide["score"],
                        is_decoy=peptide["is_decoy"],
                        peptidoform=peptide["proforma_peptidoform"],
                        modifications=peptide["modifications"],
                        peptide_ms_run=peptide["ms_run"],
                        peptide_scan_number=peptide["scan_number"],
                        peptide_position_mztab=peptide["pos"],
                    )
                    self.add_peptide_to_index(peptide_key, peptide_value)
                elif line.startswith("PSH"):
                    logger.info("-- All peptides have been read, starting psm section -- ")
                    psm_columns = line.split("\t")
                elif line.startswith("PSM"):
                    psm_info = line.split("\t")
                    es = dict(zip(psm_columns, psm_info))
                    psm = fetch_psm_from_mztab_line(
                        es,
                        ms_runs=self._ms_runs,
                        modifications_definition=self._modifications,
                    )
                    psm_key = get_key_peptide_combination(
                        psm["peptidoform"], psm["charge"]
                    )
                    self.add_psm_to_count(
                        psm_key=psm_key,
                        protein_accession=psm["accession"],
                        protein_start=psm["start"],
                        protein_end=psm["end"],
                        posterior_error_probability=psm["posterior_error_probability"],
                        calc_mass_to_charge=psm["calc_mass_to_charge"],
                        exp_mass_to_charge=psm["exp_mass_to_charge"],
                        reference_file=psm["ms_run"],
                        scan_number=psm["scan_number"],
                        mztab_postion=psm["pos"],
                    )

                line = f.readline()

    def load_mztab_file(self, use_cache: bool = False):
        """
        load a mztab file it can be in memory or in cache. if the file is in cache, it will be loaded from there.
        :param use_cache: use cache
        """
        if self._mztab_file is None:
            raise Exception("Mztab file is None")
        self.create_mztab_index(self._mztab_file)

    def print_all_peptides(self):
        """
        print all the peptides in the mztab file.
        """
        if self._use_cache:
            for key in self._peptide_index.get_all_keys():
                logger.debug("Key {} Value {}".format(key, self._peptide_index.get_item(key)))
        else:
            for peptide in self._peptide_index.keys():
                logger.debug("Key {} Value {}".format(peptide, self._peptide_index[peptide]))

    def print_mztab_stats(self):
        """
        print the mztab stats.
        """
        logger.info("Mztab stats")
        logger.info("Number of proteins {}".format(len(self._protein_details)))
        logger.info("Number of peptides {}".format(self._peptide_index.length()))
        logger.info("Number of ms runs {}".format(len(self._ms_runs)))

        if self._use_cache:
            logger.info("stats {}".format(self._peptide_index.get_stats()))

    def close(self):
        if self._use_cache:
            self._peptide_index.close()

    def get_peptide_index(self, msstats_peptidoform: str, charge: str) -> dict:
        """
        get the peptide qvalue from the peptide index.
        :param msstats_peptidoform: msstats peptidoform
        :param charge: charge
        :return: peptide qvalue object
        """
        key = get_key_peptide_combination(msstats_peptidoform, charge)
        if self._peptide_index is None:
            raise Exception("Peptide qvalues index is None")

        if self._use_cache:
            return self._peptide_index.get_item(key)
        else:
            return self._peptide_index[key]

    def get_protein_qvalue_from_index(self, protein_accession: str):
        """
        get the protein qvalue from the protein index.
        :param protein_accession: protein accession
        :return: protein qvalue object
        """
        if (
            protein_accession is None
            or self._protein_details is None
            or protein_accession not in self._protein_details
        ):
            raise Exception(
                "Protein accession to be search is None or the protein accession is not in the index"
            )
        return self._protein_details[protein_accession]

    def get_protein_qvalue_from_index_list(self, protein_accession_list: list):
        """
        get protein qvalue from the protein index using the protein acessions in a list structure.
        :param protein_accession_list list of protein accessions
        :return: qvalue structure.
        """

        if protein_accession_list is None or self._protein_details is None:
            raise Exception(
                "Protein accession to be search is None or the protein accession is not in the index"
            )

        # # If the number of proteins where the peptide maps is to big, we will have an error in the permutations.
        # if len(protein_accession_list) > 4:
        #     return None

        protein_accession_string = ";".join(protein_accession_list)
        protein_accession_string = standardize_protein_string_accession(
            protein_string=protein_accession_string, sorted=True
        )
        if protein_accession_string in self._protein_details:
            return self._protein_details[protein_accession_string]

        # # If the number of proteins where the peptide maps is to big, we will have an error in the permutations.
        if len(protein_accession_list) > 10:
            return None

        solutions = get_permutations_of_original_list(protein_accession_list)
        for accession in solutions:
            accession_string = ";".join(accession)
            if accession_string in self._protein_details:
                return self._protein_details[accession_string]
        return None

    def _get_search_engine_scores(self, line: str):
        if line.lower().__contains__("protein_search_engine_score"):
            self._search_engine_scores["protein_score"] = parse_score_name_in_mztab(
                line
            )
        elif line.lower().__contains__("peptide_search_engine_score"):
            self._search_engine_scores["peptide_score"] = parse_score_name_in_mztab(
                line
            )
        elif line.lower().__contains__("psm_search_engine_score"):
            self._search_engine_scores["psm_score"] = parse_score_name_in_mztab(line)

    def get_search_engine_scores(self):
        return self._search_engine_scores

    def get_modifications_definition(self):
        return self._modifications

    def create_mztab_psm_iterator(self, mztab_file: str):
        """
        create an iterator for the mztab file. the iterator is used to read the psm in the mztab file.
        :param mztab_file: mztab file
        :return: None
        """
        self._protein_details = {}
        self._psm_iterator = open(mztab_file, "r", buffering=calculate_buffer_size(mztab_file))
        line = self._psm_iterator.readline()
        while line != "":
            line = line.strip()
            if line.startswith("MTD") and "location" in line:
                self._ms_runs = fetch_ms_runs_from_mztab_line(line, self._ms_runs)
            if line.startswith("MTD") and "search_engine_score" in line:
                self._get_search_engine_scores(line)
            if line.startswith("MTD") and "_mod[":
                self._modifications = fetch_modifications_from_mztab_line(
                    line, self._modifications
                )
            elif line.startswith("PRH"):
                logger.info("-- End of the Metadata section of the mztab file -- ")
                protein_columns = line.split("\t")
            elif line.startswith("PRT"):
                protein_info = line.split("\t")
                logger.info("Protein info {}".format(protein_info))
                if PROTEIN_DETAILS not in line:
                    es = dict(zip(protein_columns, protein_info))
                    protein = fetch_protein_from_mztab_line(es)
                    protein["accession"] = standardize_protein_string_accession(
                        protein["accession"], sorted=True
                    )
                    self._protein_details[protein["accession"]] = [ protein["score"]]
                    logger.info("Added protein to the protein index -- {}".format(protein["accession"]))
            elif line.startswith("PSH"):
                logger.info("-- All peptides have been read, starting psm section -- ")
                self._psm_columns = line.split("\t")
                return
            line = self._psm_iterator.readline()

    def read_next_psm(self):
        """
        read the next psm in the mztab file. a precondition of the method is that the mztab file has been loaded and
        the iterator has been created using the method create_mztab_psm_iterator.
        :return: psm dictionary
        """
        if self._psm_iterator is None:
            raise Exception(
                "The mztab file has not been loaded or the iterator has not been created"
            )

        line = self._psm_iterator.readline()
        if line != "" and line.startswith("PSM"):
            psm_info = line.split("\t")
            es = dict(zip(self._psm_columns, psm_info))
            psm = fetch_psm_from_mztab_line(
                es,
                ms_runs=self._ms_runs,
                modifications_definition=self._modifications,
            )
            return psm
        elif line == "":
            self._psm_iterator.close()  # close the iterator

        return None
