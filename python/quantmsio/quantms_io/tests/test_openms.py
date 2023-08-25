from unittest import TestCase

from quantms_io.core.openms import ConsensusXMLHandler


class TestConsensusXMLHandler(TestCase):
    def test_load_consensusxml(self):
        consensusxml_path = "data/raw_ae_example/PXD009219.sdrf_openms_design_openms.consensusXML"
        consensus_xml_handler = ConsensusXMLHandler()
        consensus_xml_handler.get_intensity_map(consensusxml_path)
