from unittest import TestCase

from quantms_io.core.openms import OpenMSHandler


class TestConsensusXMLHandler(TestCase):
    def test_load_consensusxml(self):
        consensusxml_path = (
            "data/raw_ae_example/PXD009219.sdrf_openms_design_openms.consensusXML"
        )
        consensus_xml_handler = OpenMSHandler()
        consensus_xml_handler.get_intensity_map(consensusxml_path)

    def test_get_spectrum_from_scan(self):
        openms_handler = OpenMSHandler()
        spectrum = openms_handler.get_spectrum_from_scan("data/test.mzML", 1)
        print(spectrum)
