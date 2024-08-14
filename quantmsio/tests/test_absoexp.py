from unittest import TestCase

from quantmsio.tests.common import datafile
from quantmsio.core.ae import AbsoluteExpressionHander


class TestAEHandler(TestCase):

    def test_ae(self):
        (project_path, ibaq_path, sdrf_path) = map(
            datafile,
            (
                "AE/project.json",
                "AE/PXD016999.1-ibaq.tsv",
                "AE/PXD016999-first-instrument.sdrf.tsv",
            ),
        )
        ae_handler = AbsoluteExpressionHander()
        ae_handler.load_project_file(project_path)
        ae_handler.load_ibaq_file(ibaq_path)
        ae_handler.load_sdrf_file(sdrf_path)
