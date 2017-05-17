import unittest

from autoseq.util.clinseq_barcode import *


class TestClinseqBarcode(unittest.TestCase):
    def setUp(self):
        self.test_capture1 = \
            UniqueCapture("LB", "P-00000001", "CFDNA", "01234567", "TP", "CM")
        self.test_capture2 = \
            UniqueCapture("LB", "P-00000001", "CFDNA", "01234567", "TP", "WG")

    def test_extract_unique_capture_valid1(self):
        self.assertEquals(self.test_capture1,
                          extract_unique_capture("LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"))

    def test_extract_unique_capture_valid2(self):
        self.assertEquals(self.test_capture1,
                          extract_unique_capture("LB-P-00000001-CFDNA-01234567-TP1-CM1"))

    def test_extract_unique_capture_valid3(self):
        self.assertEquals(self.test_capture2,
                          extract_unique_capture("LB-P-00000001-CFDNA-01234567-TP1-WGS"))
