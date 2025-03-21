import pytest
from unittest.mock import patch
import requests
from quantmsio.core.project import ProjectHandler
from .common import datafile


test_data = [
    ("MSV000079033", "DDA-plex/MSV000079033-Blood-Plasma-iTRAQ.sdrf.tsv"),
]


@patch("requests.get")
def test_populate_from_pride_archive_successful(mock_get):
    """Test populating project info from PRIDE Archive with mocked response."""
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

    # Initialize ProjectHandler and populate from PRIDE Archive
    project_accession = "PXD123456"
    project_manager = ProjectHandler(project_accession)
    project_manager.populate_from_pride_archive()

    # Assert that the project_info has been updated
    assert project_manager.project.project_info["project_title"] == "Test Project"
    assert project_manager.project.project_info["project_description"] == "Test description"
    assert "sample_processing_protocol" in project_manager.project.project_info
    assert "data_processing_protocol" in project_manager.project.project_info


@pytest.mark.skip(reason="This test makes an actual API call to PRIDE Archive")
def test_populate_from_pride_archive_api():
    """Test populating project info from PRIDE Archive with real API call."""
    # Initialize ProjectHandler and populate from PRIDE Archive
    project_accession = "PXD020453"
    project_manager = ProjectHandler(project_accession)
    project_manager.populate_from_pride_archive()

    # Assert that the project_info has been updated
    assert project_manager.project.project_info["project_title"] == (
        "Structural insights into Cullin4-RING ubiquitin ligase remodelling by Vpr from simian immunodeficiency viruses"
    )
    assert "project_description" in project_manager.project.project_info
    assert "sample_processing_protocol" in project_manager.project.project_info
    assert "data_processing_protocol" in project_manager.project.project_info


@pytest.mark.parametrize("project_accession,sdrf_path", test_data)
def test_save_project_info(project_accession, sdrf_path):
    """Test saving project info."""
    # Resolve file path
    sdrf_file = datafile(sdrf_path)

    # Initialize ProjectHandler
    project_manager = ProjectHandler(project_accession)

    # Mock populate_from_pride_archive to avoid API call
    with patch.object(ProjectHandler, "populate_from_pride_archive") as mock_populate:
        mock_populate.return_value = None

        # Populate from SDRF
        project_manager.populate_from_sdrf(sdrf_file)

        # Assert that the project_info has been updated
        assert project_manager.project.project_info["project_accession"] == project_accession
        assert "samples" in project_manager.project.project_info
        assert len(project_manager.project.project_info["samples"]) > 0
