import re

from scipy.linalg._solve_toeplitz import float64

from quantms_io.core.core import DiskCache
from quantms_io.utils.constants import PROTEIN_DETAILS
from quantms_io.utils.pride_utils import get_key_peptide_combination, standardize_protein_string_accession, \
    parse_score_name_in_mztab, get_modifications_object_from_mztab_line


def fetch_peptide_spectra_ref(peptide_spectra_ref: str):
    """
    Get the ms run and the scan number from a spectra ref. The spectra ref is in the format:
    ms_run[1]:controllerType=0 controllerNumber=1 scan=1
    :param peptide_spectra_ref: spectra ref
    :return: ms run and scan number
    """
    ms_run = peptide_spectra_ref.split(":")[0]
    scan_number = peptide_spectra_ref.split(":")[1].split(" ")[-1].split("=")[-1]
    return ms_run, scan_number


def get_peptidoform_proforma_version_in_mztab(peptide_sequence: str, modification_string: str,
                                              modifications_definition: dict) -> str:
    """
    Get the peptidoform in propoforma notation from a mztab line definition, it could be for a PEP definition
    or for a PSM. The peptidoform in proforma notation is defined as:
     - EM[Oxidation]EVEES[Phospho]PEK
     - [iTRAQ4plex]-EMEVNESPEK-[Methyl]
    The modification string is defined as:
     - 3|4|8-UNIMOD:21, 2-UNIMO:35
     - 2[MS,MS:1001876, modification probability, 0.8]-UNIMOD:21, 3[MS,MS:1001876, modification probability, 0.8]-UNIMOD:31
     This modification means that Phospho can be in position 3, 4 or 8 and Oxidation is in position 2.
    :param peptide_sequence: Peptide sequence
    :param modification_string: modification string
    :param modifications_definition: dictionary modifications definition
    :return: peptidoform in proforma
    """
    if modification_string == "null" or modification_string is None or modification_string == "":
        return peptide_sequence

    modifications = get_modifications_object_from_mztab_line(modification_string=modification_string,
                                                             modifications_definition=modifications_definition)

    aa_index = 0
    result_peptide = ""
    peptide_sequence = list(peptide_sequence)
    # Add n-term modification if it is present
    for key_index, value_index in modifications.items():
        if aa_index in value_index["position"]:
            result_peptide = "[" + key_index + "]" + result_peptide
    if len(result_peptide) > 0:
        result_peptide = result_peptide + "-"

    aa_index += 1
    for aa in peptide_sequence:
        add_mod = ""
        for key_index, value_index in modifications.items():
            if aa_index in value_index["position"]:
                add_mod = add_mod + "[" + key_index + "]"
        result_peptide = result_peptide + aa + add_mod
        aa_index += 1
    # Add c-term modification if it is present
    for key_index, value_index in modifications.items():
        if aa_index in value_index["position"]:
            result_peptide = result_peptide + "-[" + key_index + "]"
    return result_peptide


def compare_msstats_peptidoform_with_proforma(msstats_peptidoform: str, proforma: str) -> bool:
    """
    Compare a msstats peptidoform with a proforma representation.
    - [Acetyl]-EM[Oxidation]EVEES[Phospho]PEK vs .(Acetyl)EM(Oxidation)EVEES(Phospho)PEK
    :param msstats_peptidoform: msstats peptidoform
    :param proforma: proforma peptidoform
    :return: True if the peptidoforms are the same, False otherwise
    """
    proforma = proforma.replace("[", "(").replace("]", ")")
    if "-" in proforma:
        proforma_parts = proforma.split("-")
        if len(proforma_parts) > 3:
            raise Exception("The proforma peptidoform is not valid has more than 3 termini modifications")
        result_proforma = ""
        par_index = 0
        for part in proforma_parts:
            if part.startswith("(") and par_index == 0:  # n-term modification
                result_proforma = result_proforma + "." + part
            elif part.startswith("(") and par_index > 1:
                result_proforma = result_proforma + "." + part
            else:
                result_proforma = result_proforma + part
            par_index += 1
        proforma = result_proforma
    return msstats_peptidoform == proforma


def get_petidoform_msstats_notation(peptide_sequence: str, modification_string: str,
                                    modifications_definition: dict) -> str:
    """
    Get the peptidoform in msstats notation from a mztab line definition, it could be for a PEP definition
    or for a PSM. The peptidoform in msstats notation is defined as:
        - EM(Oxidation)EVEES(Phospho)PEK
        - .(iTRAQ4plex)EMEVNESPEK.(Methyl)
    :param peptide_sequence: peptide sequence
    :param modifications_string: modification string
    :param modifications_definition: dictionary modifications definition
    :return: peptidoform in msstats
    """
    if modification_string == "null" or modification_string is None or modification_string == "":
        return peptide_sequence

    modifications = get_modifications_object_from_mztab_line(modification_string=modification_string,
                                                             modifications_definition=modifications_definition)

    aa_index = 0
    result_peptide = ""
    peptide_sequence = list(peptide_sequence)
    # Add n-term modification if it is present
    for key_index, value_index in modifications.items():
        if aa_index in value_index["position"]:
            result_peptide = "." + "(" + key_index + ")" + result_peptide

    aa_index += 1
    for aa in peptide_sequence:
        add_mod = ""
        for key_index, value_index in modifications.items():
            if aa_index in value_index["position"]:
                add_mod = add_mod + "(" + key_index + ")"
        result_peptide = result_peptide + aa + add_mod
        aa_index += 1
    # Add c-term modification if it is present
    for key_index, value_index in modifications.items():
        if aa_index in value_index["position"]:
            result_peptide = result_peptide + ".(" + key_index + ")"
    return result_peptide


def fetch_peptide_from_mztab_line(pos: int, peptide_dict: dict, ms_runs: dict = None,
                                  modification_definition: dict = None) -> dict:
    """
    Get the peptide from a mztab line include the post.
    :param pos: position of the peptide in the mztab file
    :param peptide_dict: dictionary with the peptide information
    :param ms_runs: ms runs dictionary
    :param modification_definition: modification definition
    :return: peptide dictionary
    """
    keys = ["sequence", "modifications", "charge", "retention_time", "accession"]
    peptide = dict(zip(keys, [peptide_dict[k] for k in keys]))

    peptide["accession"] = standardize_protein_string_accession(peptide["accession"])
    peptide["pos"] = pos
    peptide["score"] = peptide_dict["best_search_engine_score[1]"]

    if "opt_global_cv_MS:1002217_decoy_peptide" in peptide_dict:  # check if the peptide is a decoy
        peptide["is_decoy"] = peptide_dict["opt_global_cv_MS:1002217_decoy_peptide"]
    else:
        peptide["is_decoy"] = None  # if the information is not available

    if "opt_global_cv_MS:1000889_peptidoform_sequence" in peptide_dict:
        peptide["peptidoform"] = peptide_dict["opt_global_cv_MS:1000889_peptidoform_sequence"]
    else:
        peptide["peptidoform"] = get_petidoform_msstats_notation(peptide["sequence"], peptide["modifications"],
                                                                 modification_definition)  # if the information is not available

    peptide["proforma_peptidoform"] = get_peptidoform_proforma_version_in_mztab(peptide_sequence=peptide["sequence"],
                                                                                modification_string=peptide[
                                                                                    "modifications"],
                                                                                modifications_definition=modification_definition)

    if "spectra_ref" in peptide_dict:
        ms_run, scan_number = fetch_peptide_spectra_ref(peptide_dict["spectra_ref"])
        if ms_runs is not None and ms_run in ms_runs:
            peptide["ms_run"] = ms_runs[ms_run]
        elif ms_runs is not None and ms_run not in ms_runs:
            raise Exception("The ms run {} is not in the ms runs index".format(ms_run))
        else:
            peptide["ms_run"] = ms_run
        peptide["scan_number"] = scan_number

    else:
        peptide["ms_run"] = None  # if the information is not available
        peptide["scan_number"] = None  # if the information is not available

    comparison_representation = compare_msstats_peptidoform_with_proforma(msstats_peptidoform=peptide["peptidoform"],
                                                                          proforma=peptide["proforma_peptidoform"])

    if not comparison_representation:
        raise Exception("The proforma peptidoform {} & the msstats peptidoform {} are not equal".format(
            peptide["proforma_peptidoform"],
            peptide["peptidoform"]))
    return peptide


def fetch_protein_from_mztab_line(pos: int, protein_dict: dict):
    """
    Get the protein from a mztab line include the post.
    :param pos: position of the protein in the mztab file
    :param protein_dict: dictionary with the protein information
    :return: protein dictionary
    """
    accession_string = standardize_protein_string_accession(protein_dict["ambiguity_members"])
    return {"accession": accession_string, "score": protein_dict["best_search_engine_score[1]"], "pos": pos}


def fetch_ms_runs_from_mztab_line(mztab_line: str, ms_runs: dict) -> dict:
    """
    Get the ms runs from a mztab line. A mztab line contains the ms runs information. The structure of an msrun line
    in a mztab file is:
       MTD  ms_run[1]-location file:///C:/Users/alexis/Downloads/FileRAW.mzML
    :param mztab_line: mztab line
    :param ms_runs: ms runs dictionary
    :return: ms runs dictionary
    """
    mztab_line = mztab_line.strip()
    line_parts = mztab_line.split("\t")
    if line_parts[0] == 'MTD' and line_parts[1].split("-")[-1] == 'location':
        ms_runs[line_parts[1].split("-")[0]] = line_parts[2].split("//")[-1].split(".")[0]
    return ms_runs


def fetch_psm_from_mztab_line(pos: int, es: dict, ms_runs: dict = None, modifications_definition: dict = None) -> dict:
    """
    Get the psm from a mztab line include the post.
    :param pos: Position of the psm in the mztab file
    :param es: dictionary with the psm information
    :param ms_runs: ms runs dictionary
    :param modifications_definition: modifications definition
    :return: psm dictionary
    """
    keys = ["sequence", "modifications", "charge", "retention_time", "accession", "start", "end", "calc_mass_to_charge",
            "exp_mass_to_charge"]
    psm = dict(zip(keys, [es[k] for k in keys]))

    psm["accession"] = standardize_protein_string_accession(psm["accession"])
    psm["pos"] = pos
    psm["score"] = es["search_engine_score[1]"]

    if "opt_global_cv_MS:1002217_decoy_peptide" in es:  # check if the peptide is a decoy
        psm["is_decoy"] = es["opt_global_cv_MS:1002217_decoy_peptide"]
    else:
        psm["is_decoy"] = None  # if the information is not available

    if "opt_global_Posterior_Error_Probability_score" in es:
        psm["posterior_error_probability"] = es["opt_global_Posterior_Error_Probability_score"]
    else:
        psm["posterior_error_probability"] = None

    if "opt_global_q-value" in es:
        psm["global_qvalue"] = es["opt_global_q-value"]
    else:
        psm["global_qvalue"] = None

    if "opt_global_cv_MS:1000889_peptidoform_sequence" in es:
        psm["peptidoform"] = es["opt_global_cv_MS:1000889_peptidoform_sequence"]
    else:
        psm["peptidoform"] = get_petidoform_msstats_notation(psm["sequence"], psm["modifications"],
                                                             modifications_definition)  # if the information is not available

    psm["proforma_peptidoform"] = get_peptidoform_proforma_version_in_mztab(psm["sequence"], psm["modifications"],
                                                                            modifications_definition)
    if "spectra_ref" in es:
        ms_run, scan_number = fetch_peptide_spectra_ref(es["spectra_ref"])
        if ms_runs is not None and ms_run in ms_runs:
            psm["ms_run"] = ms_runs[ms_run]
        elif ms_runs is not None and ms_run not in ms_runs:
            raise Exception("The ms run {} is not in the ms runs index".format(ms_run))
        else:
            psm["ms_run"] = ms_run
        psm["scan_number"] = scan_number
    else:
        psm["ms_run"] = None  # if the information is not available
        psm["scan_number"] = None  # if the information is not available

    return psm


def fetch_modifications_from_mztab_line(line: str, _modifications: dict) -> dict:
    """
    Get the modifications from a mztab line. An mzTab modification could be a fixed or variable modification.
    The structure of a fixed is the following:
      MTD	fixed_mod[1]	[UNIMOD, UNIMOD:4, Carbamidomethyl, ]
      MTD	fixed_mod[1]-site	C
      MTD	fixed_mod[1]-position	Anywhere
    while the structure of a variable modification is the following:
      MTD	var_mod[1]	[UNIMOD, UNIMOD:21, Phospho, ]
      MTD	var_mod[1]-site	S
      MTD   var_mod[1]-position	Anywhere

    :param line: mztab line
    :param _modifications: modifications dictionary
    :return: modification dictionary
    """
    line = line.strip()
    line_parts = line.split("\t")
    if line_parts[0] == 'MTD' and "_mod[" in line_parts[1]:
        if "site" not in line_parts[1] and "position" not in line_parts[1]:
            values = line_parts[2].replace("[", "").replace("]", "").split(",")
            accession = values[1].strip()
            name = values[2].strip()
            index = line_parts[1].split("[")[1].split("]")[0]
            _modifications[accession] = [name, index, None, None]
        elif "site" in line_parts[1]:
            index = line_parts[1].split("[")[1].split("]")[0]
            accession = None
            for key, value in _modifications.items():  # for name, age in dictionary.iteritems():  (for Python 2.x)
                if value[1] == index:
                    accession = key
            if accession is None:
                raise Exception("The accession for the modification is None")
            _modifications[accession][2] = line_parts[2]
        elif "position" in line_parts[1]:
            index = line_parts[1].split("[")[1].split("]")[0]
            accession = None
            for key, value in _modifications.items():
                if value[1] == index:
                    accession = key
            if accession is None:
                raise Exception("The accession for the modification is None")
            _modifications[accession][3] = line_parts[2]
    return _modifications


def create_peptide_for_index(peptide_qvalue: str, is_decoy: str, peptidoform: str,
                             modifications: str, peptide_ms_run: str, peptide_scan_number: str,
                             peptide_position_mztab: str, psm_protein_start: str = None, psm_protein_end: str = None,
                             posterior_error_probability: str = None, spectral_count: dict = None,
                             list_psms: list = None) -> dict:
    """
    Create a peptide to be store in the index. A index peptide have the following structure:
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
     :param protein_accession: protein accession
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
    peptide = {"peptide_qvalue": peptide_qvalue, "is_decoy": is_decoy,
               "peptidoform": peptidoform, "modifications": modifications, "peptide_ms_run": peptide_ms_run,
               "peptide_scan_number": peptide_scan_number, "peptide_position_mztab": peptide_position_mztab,
               "psm_protein_start": psm_protein_start, "psm_protein_end": psm_protein_end,
               "spectral_count": spectral_count, "list_psms": list_psms,
               "posterior_error_probability": posterior_error_probability}
    return peptide


class MztabHandler:
    """
    Mztab handler class for quantms.io. This class is used to load mztab files and create indexes for the mztab files.
    The mztab handler has multiple indexes to access the mztab content: protein, peptide and psm.
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
        self._index = {}  # index for each section of the mztab file

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
        Add a peptide to the index if the file is big then use cache.
        :param peptide_key: peptide key
        :param peptide_value: peptide value
        :return: None
        """
        if self._use_cache:
            self._peptide_index.add_item(peptide_key, peptide_value)
        else:
            self._peptide_index[peptide_key] = peptide_value

    def add_psm_to_count(self, psm_key: str, protein_accession: str, protein_start: str, protein_end: str,
                         calc_mass_to_charge: str, exp_mass_to_charge: str, reference_file: str, scan_number: str,
                         mztab_postion: str, posterior_error_probability: str = None):
        """
        PSM count is used to count the number of psms for a given peptidoform in a file.
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
        :param global_qvalue: global qvalue
        :return: None
        """
        if self._use_cache:
            if self._peptide_index.contains(psm_key):
                psm = self._peptide_index.get_item(psm_key)
                peptide_protein_accession = psm["protein_accession"] if "protein_accession" in psm else None
                peptide_protein_start = psm["psm_protein_start"] if "psm_protein_start" in psm else None
                peptide_protein_end = psm["psm_protein_end"] if "psm_protein_end" in psm else None
                if posterior_error_probability is not None:
                    if (psm["posterior_error_probability"] is None or
                            (float64(posterior_error_probability) < psm["posterior_error_probability"])):
                           psm["posterior_error_probability"] = float64(posterior_error_probability)

                if ((peptide_protein_start is not None and peptide_protein_end is not None)
                        and [peptide_protein_accession, peptide_protein_start, peptide_protein_end] != [
                            protein_accession, protein_start, protein_end]):
                    raise Exception("The protein information is different for the same psm key")
                else:
                    psm["psm_protein_start"] = protein_start
                    psm["psm_protein_end"] = protein_end
                    psm["protein_accession"] = protein_accession

                spectral_count = psm["spectral_count"]  # dictionary with the spectral count reference file: spectral count
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
                list_psms.append([calc_mass_to_charge, exp_mass_to_charge, reference_file, scan_number, mztab_postion])
                psm["list_psms"] = list_psms

                self._peptide_index.add_item(psm_key, psm)
            # else:
                # print("The psm key {} is not in the peptide index".format(psm_key))
        else:
            if psm_key in self._peptide_index:
                psm = self._peptide_index[psm_key]
                peptide_protein_accession = psm["protein_accession"] 
                peptide_protein_start = psm["psm_protein_start"] 
                peptide_protein_end = psm["psm_protein_end"] 
                if posterior_error_probability is not None and posterior_error_probability < psm["posterior_error_probability"]:
                    psm["posterior_error_probability"] = posterior_error_probability
                if ((peptide_protein_start is not None and peptide_protein_end is not None)
                        and [peptide_protein_accession, peptide_protein_start, peptide_protein_end] != [
                            protein_accession, protein_start,
                            protein_end]):
                    raise Exception("The protein information is different for the same psm key")
                else:
                    psm["psm_protein_start"] = protein_start
                    psm["psm_protein_end"] = protein_end

                spectral_count = psm[
                    "spectral_count"]  # dictionary with the spectral count reference file: spectral count
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
                list_psms.append([calc_mass_to_charge, exp_mass_to_charge, reference_file, scan_number])
                psm["list_psms"] = list_psms

                self._peptide_index[psm_key] = psm
            # else:
            #     print("The psm key {} is not in the peptide index".format(psm_key))

    def create_mztab_index(self, mztab_file: str, qvalue_index: bool = True, psm_count_index: bool = True):
        """
        Create an index for a mztab file; the index contains a structure with the position of each psm, peptide and
        protein in the file.
        :param mztab_file: mztab file
        """
        self._protein_details = {}
        self.create_peptide_index()

        with open(mztab_file) as f:
            pos = 0
            line = f.readline()
            while line != "":
                line = line.strip()
                if line.startswith("MTD") and "location" in line:
                    self._ms_runs = fetch_ms_runs_from_mztab_line(line, self._ms_runs)
                if line.startswith("MTD") and "search_engine_score" in line:
                    self._get_search_engine_scores(line)
                if line.startswith("MTD") and "_mod[":
                    self._modifications = fetch_modifications_from_mztab_line(line, self._modifications)
                elif line.startswith("PRH"):
                    print("-- End of the Metadata section of the mztab file -- ")
                    protein_columns = line.split("\t")
                    self._index["PRH"] = pos
                elif line.startswith("PRT"):
                    self._index["PRT"] = pos
                    protein_info = line.split("\t")
                    if PROTEIN_DETAILS not in line:
                        es = dict(zip(protein_columns, protein_info))
                        protein = fetch_protein_from_mztab_line(pos, es)
                        self._protein_details[protein["accession"]] = [protein["score"], pos]
                elif line.startswith("PEH"):
                    print("-- All proteins have been read, starting peptide section -- ")
                    self._index["PEH"] = pos
                    peptide_columns = line.split("\t")
                elif line.startswith("PEP"):
                    self._index["PEP"] = pos
                    peptide_info = line.split("\t")
                    es = dict(zip(peptide_columns, peptide_info))
                    peptide = fetch_peptide_from_mztab_line(pos, es, ms_runs=self._ms_runs,
                                                            modification_definition=self._modifications)
                    peptide_key = get_key_peptide_combination(msstats_peptidoform=peptide["peptidoform"],
                                                              charge=peptide["charge"])

                    peptide_value = create_peptide_for_index(peptide_qvalue=peptide["score"],
                                                             is_decoy=peptide["is_decoy"],
                                                             peptidoform=peptide["proforma_peptidoform"],
                                                             modifications=peptide["modifications"],
                                                             peptide_ms_run=peptide["ms_run"],
                                                             peptide_scan_number=peptide["scan_number"],
                                                             peptide_position_mztab=peptide["pos"])
                    self.add_peptide_to_index(peptide_key, peptide_value)
                elif line.startswith("PSH"):
                    print("-- All peptides have been read, starting psm section -- ")
                    self._index["PSH"] = pos
                    psm_columns = line.split("\t")
                elif line.startswith("PSM"):
                    self._index["PSM"] = pos
                    psm_info = line.split("\t")
                    es = dict(zip(psm_columns, psm_info))
                    psm = fetch_psm_from_mztab_line(pos, es, ms_runs=self._ms_runs,
                                                    modifications_definition=self._modifications)
                    psm_key = get_key_peptide_combination(psm["peptidoform"], psm["charge"])
                    self.add_psm_to_count(psm_key=psm_key, protein_accession=psm["accession"],
                                          protein_start=psm["start"], protein_end=psm["end"],
                                          posterior_error_probability=psm["posterior_error_probability"],
                                          calc_mass_to_charge=psm["calc_mass_to_charge"],
                                          exp_mass_to_charge=psm["exp_mass_to_charge"],
                                          reference_file=psm["ms_run"],
                                          scan_number=psm["scan_number"],
                                          mztab_postion=psm["pos"]
                                          )

                pos = f.tell()
                line = f.readline()
        return self._index

    def load_mztab_file(self, use_cache: bool = False):
        """
        Load a mztab file it can be in memory or in cache. If the file is in cache, it will be loaded from there.
        """
        if self._mztab_file is None:
            raise Exception("Mztab file is None")
        self.create_mztab_index(self._mztab_file)

    def print_all_peptides(self):
        """
        Print all the peptides in the mztab file.
        """
        if self._use_cache:
            for key in self._peptide_index.get_all_keys():
                print("Key {} Value {}".format(key, self._peptide_index.get_item(key)))
        else:
            for peptide in self._peptide_index.keys():
                print("Key {} Value {}".format(peptide, self._peptide_index[peptide]))

    def print_mztab_stats(self):
        """
        Print the mztab stats.
        """
        print("Mztab stats")
        print("Number of proteins {}".format(len(self._protein_details)))
        print("Number of peptides {}".format(self._peptide_index.length()))
        print("Number of ms runs {}".format(len(self._ms_runs)))

        if self._use_cache:
            print("stats {}".format(self._peptide_index.get_stats()))

    def close(self):
        if self._use_cache:
            self._peptide_index.close()

    def get_peptide_index(self, msstats_peptidoform: str, charge: str) -> dict:
        """
        Get the peptide qvalue from the peptide index.
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
        Get the protein qvalue from the protein index.
        :param protein_accession: protein accession
        :return: protein qvalue object
        """
        if protein_accession is None or self._protein_details is None or protein_accession not in self._protein_details:
            raise Exception("Protein accession to be search is None or the protein accession is not in the index")
        return self._protein_details[protein_accession]

    def _get_search_engine_scores(self, line: str):
        if line.lower().__contains__("protein_search_engine_score"):
            self._search_engine_scores["protein_score"] = parse_score_name_in_mztab(line)
        elif line.lower().__contains__("peptide_search_engine_score"):
            self._search_engine_scores["peptide_score"] = parse_score_name_in_mztab(line)
        elif line.lower().__contains__("psm_search_engine_score"):
            self._search_engine_scores["psm_score"] = parse_score_name_in_mztab(line)

    def get_search_engine_scores(self):
        return self._search_engine_scores

    def get_modifications_definition(self):
        return self._modifications
