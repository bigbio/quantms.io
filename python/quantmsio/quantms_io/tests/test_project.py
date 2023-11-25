from unittest import TestCase
from unittest.mock import patch

import requests

from quantms_io.core.project import ProjectHandler


class TestProjectHandler(TestCase):
    @patch("requests.get")
    def test_populate_from_pride_archive_successful(self, mock_get):
        # Mock the API response for a successful request
        mock_response = requests.models.Response()
        mock_response.status_code = 200
        mock_response.json = lambda: {
            "title": "Test Project",
            "projectDescription": "Test description",
            "sampleProcessingProtocol": "Test sample processing protocol",
            "dataProcessingProtocol": "Test data processing protocol",
            "references": [{"pubmedId": "12345678"}],
            "keywords": ["keyword1", "keyword2"],
            "projectTags": ["tag1", "tag2"],
            # Add other mock data as needed
        }
        mock_get.return_value = mock_response

        project_accession = "PXD123456"
        project_manager = ProjectHandler(project_accession)
        project_manager.populate_from_pride_archive()

        # Assert that the project_info has been updated
        self.assertEqual(
            project_manager.project.project_info["project_title"], "Test Project"
        )
        self.assertEqual(
            project_manager.project.project_info["project_description"],
            "Test description",
        )

    def test_populate_from_pride_archive_api(self):
        project_accession = "PXD020453"
        project_manager = ProjectHandler(project_accession)
        project_manager.populate_from_pride_archive()

        # Assert that the project_info has been updated
        self.assertEqual(
            project_manager.project.project_info["project_title"],
            "Structural insights into Cullin4-RING ubiquitin ligase remodelling by Vpr from simian immunodeficiency viruses",
        )
        self.assertEqual(
            project_manager.project.project_info["project_description"],
            "crosslinking mass spectrometry results for sulfo-SDA crosslinking of human CUL4-NEDD8/ROC1/DDB1/DCAF1-CtD in complex with SAMHD1 and Vpr protein from simian immunodeficiency virus infecting Cercopithecus cephus (SIVmus Vpr)",
        )
        print(project_manager.project.project_info)

    def test_save_project_info(self):
        project_accession = "PXD020187"
        sdrf_file = "data/PXD020187.sdrf.tsv"

        project_manager = ProjectHandler(project_accession)
        project_manager.populate_from_pride_archive()
        project_manager.populate_from_sdrf(sdrf_file)

        project_manager.save_project_info()  # Save the project information to a JSON file

    def test_load_project_from_json(self):
        project_file = (
            "data/PXD020187-934d85ce-6a7a-4417-8330-a21b750fd9e4.project.json"
        )
        project_manager = ProjectHandler(project_json_file=project_file)
        project_manager.save_project_info()  # Save the project information to a JSON file
