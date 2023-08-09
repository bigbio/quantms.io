import json
import uuid

import requests

from quantms_io.core.sdrf import SDRFHandler
from quantms_io.utils.pride_utils import get_pubmed_id_pride_json, get_set_of_experiment_keywords


class ProjectHandler:
    def __init__(self, project_accession: str):
        self.project_accession = project_accession
        self.project = ProjectDefinition()

    def _load_project_info(self, project_json_file: str = None):
        """
        Load the project information from a JSON file
        :param project_json_file: JSON file with the project information
        """
        try:
            with open(f"{project_json_file}", "r") as json_file:
                json_dict = json.load(json_file)
                self.project = ProjectDefinition()
                self.project.set_project_info(json_dict)
        except FileNotFoundError:
            self.project = ProjectDefinition()

    def populate_from_pride_archive(self):
        # Simulate API request to PRIDE Archive
        api_url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{self.project_accession}"
        response = requests.get(api_url)

        if response.status_code == 200:
            pride_data = response.json()
            self.project.project_info["project_accession"] = self.project_accession
            self.project.project_info["project_title"] = pride_data["title"]
            self.project.project_info["project_description"] = pride_data["projectDescription"]
            self.project.project_info["project_sample_description"] = pride_data["sampleProcessingProtocol"]
            self.project.project_info["project_data_description"] = pride_data["dataProcessingProtocol"]
            self.project.project_info["project_pubmed_id"] = get_pubmed_id_pride_json(pride_data)
            self.project.project_info["experiment_type"] = get_set_of_experiment_keywords(pride_data)
        else:
            print(f"Error retrieving data from PRIDE Archive API. Status code: {response.status_code}")

    def add_quantms_version(self, quantms_version: str):
        """
        Add the quantms version to the project information
        :param quantms_version: QuantMS version
        """
        self.project.project_info["quantms_version"] = quantms_version

    def add_sdrf_project_properties(self, sdrf: SDRFHandler):
        """
        Add the project properties from the SDRF file to the project information.
        :param sdrf: SDRFHandler object
        """
        self.project.project_info["organisms"] = sdrf.get_organisms()
        self.project.project_info["organism_parts"] = sdrf.get_organism_parts()
        self.project.project_info["diseases"] = sdrf.get_diseases()
        self.project.project_info["cell_lines"] = sdrf.get_cell_lines()
        self.project.project_info["instruments"] = sdrf.get_instruments()
        self.project.project_info["enzymes"] = sdrf.get_enzymes()
        self.project.project_info["acquisition_properties"] = sdrf.get_acquisition_properties()

    def add_quantms_file(self, file_section: str, file_extension: str):
        """
        Add a quantms file to the project information. The file name will be generated automatically. Read more about the
        quantms file naming convention in the docs folder of this repository
        (https://github.com/bigbio/quantms.io/blob/main/docs/PROJECT.md)
        :param file_section: Section of the quantms file (protein, peptide, psm, feature, absolute, differential, sdrf, etc.)
        :param file_extension: File extension of the quantms file (json, csv, tsv, etc.)
        """
        file_name = f"{self.project_accession}-{str(uuid.uuid4())}.{file_section}.{file_extension}"
        if "quantms_files" not in self.project.project_info:
            self.project.project_info["quantms_files"] = []
        self.project.project_info["quantms_files"].append({file_section + "_file": file_name})

    def save_project_info(self):
        """
        Save the updated project information to a JSON file. The file name will be generated automatically.
        """
        # Save the updated project info to a JSON file
        output_filename = f"{self.project_accession}-{str(uuid.uuid4())}.project.json"
        with open(output_filename, "w") as json_file:
            json.dump(self.project.project_info, json_file, indent=4)
        print(f"Updated project information saved to {output_filename}")

    def populate_from_sdrf(self, sdrf_file: str):
        """
        Populate the project information from an SDRF file using the SDRFHandler class.
        :param sdrf_file: SDRF file
        """
        sdrf = SDRFHandler(sdrf_file)
        self.add_sdrf_project_properties(sdrf)


class ProjectDefinition:
    """
    Class to define the project information. This class will be used to generate the project JSON file. Read more about
    the project JSON file in the docs folder of this repository
    (https://github.com/bigbio/quantms.io/blob/main/docs/PROJECT.md).
    """

    def __init__(self, project_accession: str = None):
        self.project_info = {
            "project_accession": "",
            "project_title": "",
            "project_description": "",
            "project_sample_description": "",
            "project_data_description": "",
            "project_pubmed_id": "",
            "organisms": [],
            "organism_parts": [],
            "diseases": [],
            "cell_lines": [],
            "instruments": [],
            "enzymes": [],
            "experiment_type": [],
            "acquisition_properties": [],
            "quantms_files": [],
            "quantms_version": "",
            "comments": []
        }
        if project_accession:
            self.set_basic_info(project_accession)

    def set_basic_info(self, project_accession: str):
        self.project_info["project_accession"] = project_accession

    # Define additional methods to set specific fields as needed

    def get_project_info(self):
        return self.project_info

    def set_project_info(self, project_info: dict):
        self.project_info = project_info
