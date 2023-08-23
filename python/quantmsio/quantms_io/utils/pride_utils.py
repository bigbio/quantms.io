"""
This file contains utility functions for parsing PRIDE JSON files
"""
import re

def get_pubmed_id_pride_json(pride_json: dict) -> str:
    """
    Parse the PubMed ID from the PRIDE JSON file
    :param pride_json: PRIDE JSON file
    :return: PubMed ID
    """
    pubmed_id = None
    if "references" in pride_json:
        for reference in pride_json["references"]:
            if "pubmedId" in reference:
                pubmed_id = reference["pubmedId"]
    return pubmed_id


def get_set_of_experiment_keywords(pride_json: dict) -> list:
    """
    Parse the experiment type from the PRIDE JSON file
    :param pride_json: PRIDE JSON file
    :return: Set of experiment types
    """
    experiment_types = set()
    if "projectTags" in pride_json:
        for experiment_type in pride_json["projectTags"]:
            experiment_types.add(experiment_type)
    if "keywords" in pride_json:
        for experiment_type in pride_json["keywords"]:
            experiment_types.add(experiment_type)
    return list(experiment_types)


def get_key_peptide_combination(msstats_peptidoform: str, charge: str, reference_file: str = None):
    """
    Get the key for a peptide or psm in mztab. For peptides the key is:
     - sequence peptidoform in msstats notation (e.g. PEPTIDE(Oxidation))
     - charge (e.g. 2)
    For psms the key is:
     - sequence peptidoform in msstats notation (e.g. PEPTIDE(Oxidation))
     - charge (e.g. 2)
     - reference_file (e.g. Run1-RAW)
    :param msstats_peptidoform: sequence peptidoform in msstats notation (e.g. PEPTIDE(Oxidation))
    :param charge: charge (e.g. 2)
    :param reference_file: reference file (e.g. Run1-RAW)
    :return: key
    """
    if reference_file is not None:
        # We use a different separator for the key to avoid conflicts with the separator of the msruns
        # which can sometimes use the separator "_".
        return msstats_peptidoform + ":_:" + charge + ":_:" + reference_file
    return msstats_peptidoform + "_" + charge


def decompose_key_peptide_combination(key: str):
    """
    Decompose the key for a peptide in mztab. The key is a combination of msstats_peptidoform,
    charge and in some cases the reference file.
    :param key: key
    :return: sequence, modification, charge and retention time
    """
    if ":_:" in key:
        return key.split(":_:")
    return key.split("_")

def clean_peptidoform_sequence(sequence: str) -> str:
    """
    Clean any peptidoform a proforma sequence or a msstats sequence. This function removes any modification from the
    sequence.
    :param sequence: sequence
    """
    if "(" not in sequence and "[" not in sequence: # no modification
        return sequence

    # remove any modification if a proforma sequence is provided
    sequence = re.sub("[\(\[].*?[\)\]]", "", sequence)
    sequence = sequence.replace(".", "").replace(" ", "").replace("-", "")
    return sequence

def compare_protein_lists(protein_list_1: list, protein_list_2: list) -> bool:
    """
    Compare two protein lists
    :param protein_list_1: protein list 1
    :param protein_list_2: protein list 2
    :return: True if the three protein lists are the same
    """
    protein_list_1 = set([x.strip() for x in protein_list_1])
    protein_list_2 = set([x.strip() for x in protein_list_2])
    return protein_list_1 == protein_list_2

def standardize_protein_string_accession(protein_string: str) -> str:
    """
    Standardize the protein string accession, in some cases the protein string accession is decided by commas
    instead of semicolons.
    :param protein_string: protein string
    :return: standardized protein string
    """
    return protein_string.replace(",", ";").strip()

def standardize_protein_list_accession(protein_string: str) -> list:
    """
    Get the list of protein accessions from a protein string join by semicolons.
    :param protein_string: protein string
    :return: list of protein accessions
    """
    return [x.strip() for x in protein_string.split(";")]

def parse_score_name_in_mztab(score_name_mztab_line: str) -> str:
    """
    Parse the score name in mztab. The score name in mztab is a combination of the score name and the score type.
    :param score_name_mztab_line: score name in mztab
    :return: score name
    """
    lines = score_name_mztab_line.split("\t")
    score_values = lines[2].replace("[", "").replace("]", "").split(",")
    score_name = score_values[2].strip()
    if ":" in score_name:
        score_name = "'{}'".format(score_name) # add quotes to the score name if it contains a colon like
        # "OpenMS:Target-decoy protein q-value"
    return score_name

def get_modifications_object_from_mztab_line(modification_string: str, modifications_definition: dict) -> dict:
    """
    Get the modifications from a mztab line. This method is used to transform peptide + modification strings to
    proteoform notations, for msstats notation and for proforma notation.
    :param modification_string: modification string
    :param modifications_definition: dictionary modifications definition
    :return: modifications dictionary
    """
    modifications = {}
    modification_values = re.split(r',(?![^\[]*\])', modification_string)
    for modification in modification_values:
        modification = modification.strip()
        accession = modification.split("-")[1]
        unimod_accession = accession
        if accession not in modifications_definition:
            raise Exception("The modification {} is not in the modifications definition".format(accession))
        accession = modifications_definition[accession][0]  # get the name of the modification
        position = []
        position_probability_string = modification.split("-")[0]
        if "[" not in position_probability_string and "|" not in position_probability_string:  # only one position
            position = [position_probability_string]
        elif "[" not in position_probability_string and "|" in position_probability_string:  # multiple positions not probability
            position = position_probability_string.split("|")  # multiple positions not probability
        else:
            positions_probabilities = position_probability_string.split("|")
            for position_probability in positions_probabilities:
                if "[" not in position_probability:
                    position.append(position_probability)
                else:
                    position_with_probability = position_probability.split("[")[0]
                    position.append(position_with_probability)
        position = [int(i) for i in position]
        if accession in modifications:  # Avoid error in OpenMS that do not write 1|4-UNIMOD:35 but 1-UNIMOD:35, 4-UNIMOD:35
            position = modifications[accession]["position"] + position
            modifications[accession] = {"position":position, "unimod_accession": unimod_accession}
        else:
            modifications[accession] = {"position": position, "unimod_accession":unimod_accession}
    return modifications

