"""
This file contains utility functions for parsing PRIDE JSON files
"""


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
