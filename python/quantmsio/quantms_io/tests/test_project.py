from unittest import TestCase
from unittest.mock import patch

import requests
from ddt import data
from ddt import ddt

from quantms_io.core.project import ProjectHandler


@ddt
class TestProjectHandler(TestCase):
    global test_datas
    test_datas = [
        ("PXD040438", "/examples/DDA-lfq/PXD040438.sdrf.tsv", "/examples/output/DDA-lfq/"),
        ("MSV000079033", "/examples/DDA-plex/MSV000079033-Blood-Plasma-iTRAQ.sdrf.tsv", "/examples/output/DDA-plex/"),
        ("PXD037682", "/examples/DIA-lfq/PXD037682.sdrf.tsv", "/examples/output/DIA-lfq/"),
    ]

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
        self.assertEqual(project_manager.project.project_info["project_title"], "Test Project")
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

    @data(*test_datas)
    def test_save_project_info(self, test_data):
        project_accession = test_data[0]
        sdrf_file = __package__ + test_data[1]
        output_folder = __package__ + test_data[2]

        project_manager = ProjectHandler(project_accession)
        project_manager.populate_from_pride_archive()
        project_manager.populate_from_sdrf(sdrf_file)
        project_manager.add_sdrf_file(
            sdrf_file_path=sdrf_file,
            output_folder=output_folder,
            delete_existing=False,
        )
        project_manager.save_updated_project_info(output_file_name=output_folder + "peoject.json")
