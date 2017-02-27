import unittest

from mock import patch, Mock
from autoseq.pipeline.alascca import AlasccaPipeline


 class TestAlascca(unittest.TestCase):
    def setUp(self):
        self.dummy_path_to_big_vcf = "dummy_path_to_big_vcf"
        self.dummy_path_to_V3V4_vcf = "dummy_path_to_V3V4_vcf"
        self.dummy_path_to_V3V4big_vcf = "dummy_path_to_V3V4big_vcf"

        dummy_sample_data = {}
        dummy_reference_data = {"contest_vcfs":{"big": self.dummy_path_to_big_vcf,
                                                "clinseqV3V4": self.dummy_path_to_V3V4_vcf,
                                                "clinseqV3V4big": self.dummy_path_to_V3V4big_vcf}}

        # Construct a AlasccaPipeline that is set up well enough for unit
        # testing to be run on it's instance methods:
        self.dummy_alascca_pipeline = AlasccaPipeline(sampledata=dummy_sample_data,
                                                      refdata=dummy_reference_data,
                                                      outdir="/tmp/autoseq-test",
                                                      libdir="/tmp",
                                                      maxcores=1,
                                                      analysis_id="1",
                                                      runner="shellrunner",
                                                      jobdb=None,
                                                      dot_file=None,
                                                      scratch="/tmp"
                                                      )

    @patch('autoseq.util.library.get_libdict')
    def test_get_pop_af_vcf_both_cb(self, get_libdict):
        get_libdict = Mock()
        get_libdict.side_effect = ["CB", "CB"]

        self.assertEquals(self.dummy_alascca_pipeline.get_pop_af_vcf(), self.dummy_path_to_big_vcf)
#
#    def test_get_pop_af_vcf_tpanel_cs_npanel_cs(self):
