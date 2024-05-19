from unittest import TestCase

from core.sdrf import SDRFHandler


class TestSDRFHandler(TestCase):
    def test__load_sdrf_info(self):
        file = __package__ + "/examples/DDA-lfq/PXD040438.sdrf.tsv"
        sdrf_handler = SDRFHandler(file)
        sdrf_handler._load_sdrf_info(file)

        self.assertEqual(sdrf_handler.get_organisms(), ["Homo sapiens"])
        self.assertEqual(sdrf_handler.get_instruments(), ["Electrospray ionization"])
        self.assertEqual(sdrf_handler.get_diseases(), ["COVID-19", "Normal"])
        self.assertEqual(sdrf_handler.get_enzymes(), ["Trypsin"])
        self.assertEqual(sdrf_handler.get_cell_lines(), [])
        self.assertEqual(
            sdrf_handler.get_acquisition_properties(),
            [
                {"proteomics data acquisition method": "Label free"},
                {"proteomics data acquisition method": "Data-dependent acquisition"},
                {"dissociation method": "HCD"},
                {"precursor mass tolerance": "20 ppm"},
                {"fragment mass tolerance": "0.6 Da"},
            ],
        )

        values_df = sdrf_handler.extract_feature_properties()
        print(values_df)

        print(sdrf_handler.get_experiment_type_from_sdrf())

    def test_get_labels(self):
        file = __package__ + "/examples/DDA-plex/MSV000079033-Blood-Plasma-iTRAQ.sdrf.tsv"
        sdrf_handler = SDRFHandler(file)
        self.assertEqual(len(sdrf_handler.get_sample_labels()), 4)

        experiment_type = sdrf_handler.get_experiment_type_from_sdrf()
        print(experiment_type)
