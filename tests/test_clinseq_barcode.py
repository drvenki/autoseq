import unittest
from mock import mock_open, patch
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

    def test_extract_unique_capture_invalid(self):
        self.assertRaises(ValueError, lambda: extract_unique_capture("an_invalid_barcode"))

    def test_clinseq_barcode_is_valid_valid1(self):
        self.assertTrue(clinseq_barcode_is_valid("LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"))

    def test_clinseq_barcode_is_valid_invalid_too_few(self):
        self.assertFalse(clinseq_barcode_is_valid("LB-00000001-CFDNA-01234567-TP-CM2017001022000"))

    def test_clinseq_barcode_is_valid_invalid_project(self):
        self.assertFalse(clinseq_barcode_is_valid("01-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"))

    def test_sdid_valid_valid1(self):
        self.assertTrue(sdid_valid("P-A"))

    def test_sdid_valid_invalid1(self):
        self.assertFalse(sdid_valid("P-$"))

    def test_prep_id_valid_valid1(self):
        self.assertTrue(prep_id_valid("AA1"))

    def test_prep_id_valid_invalid_too_short(self):
        self.assertFalse(prep_id_valid("AA"))

    def test_prep_id_valid_invalid_no_number(self):
        self.assertFalse(prep_id_valid("AAZ"))

    def test_capture_id_valid_valid_wgs(self):
        self.assertTrue(capture_id_valid("WGS"))

    def test_capture_id_valid_invalid_no_number(self):
        self.assertFalse(capture_id_valid("AAZ"))

    def test_capture_id_valid_invalid_too_short(self):
        self.assertFalse(capture_id_valid("AA"))

    def test_extract_kit_id_valid(self):
        self.assertEquals(extract_kit_id("SomeString"), "So")

    def test_extract_kit_id_invalid(self):
        self.assertRaises(ValueError, lambda: extract_kit_id("C"))

    def test_extract_clinseq_barcodes_invalid_file_extension(self):
        self.assertRaises(ValueError, lambda: extract_clinseq_barcodes("some_file.invalid_extension"))

    @patch('autoseq.util.clinseq_barcode.parse_orderform')
    def test_extract_clinseq_barcodes_xlsx(self, mock_parse_orderform):
        mock_parse_orderform.return_value = ["a_mock_barcode", "another_mock_barcode"]
        self.assertEquals(len(extract_clinseq_barcodes("test.xlsx")), 2)

    def test_extract_clinseq_barcodes_txt(self):
        mocked_open = mock_open(read_data='a_mock_barcode\nanother_mock_barcode\n')
        with patch('autoseq.util.clinseq_barcode.open', mocked_open, create=True):
            observed_len = len(extract_clinseq_barcodes("test.txt"))
            self.assertEquals(observed_len, 2)

    @patch('autoseq.util.clinseq_barcode.os.path.exists')
    @patch('autoseq.util.clinseq_barcode.find_fastqs')
    def test_data_available_for_clinseq_barcode_file_exists(self, mock_find_fastqs, mock_os_path_exists):
        mock_find_fastqs.return_value = True
        mock_os_path_exists.return_value = True
        self.assertTrue(
            data_available_for_clinseq_barcode("test_libdir",
                                               "LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"))

    @patch('autoseq.util.clinseq_barcode.os.path.exists')
    @patch('autoseq.util.clinseq_barcode.find_fastqs')
    def test_data_available_for_clinseq_barcode_invalid_barcode(self, mock_find_fastqs, mock_os_path_exists):
        mock_find_fastqs.return_value = True
        mock_os_path_exists.return_value = True
        self.assertRaises(ValueError, lambda: data_available_for_clinseq_barcode("test_libdir",
                                                                                 "an_invalid_barcode"))

    @patch('autoseq.util.clinseq_barcode.os.path.exists')
    @patch('autoseq.util.clinseq_barcode.find_fastqs')
    def test_data_available_for_clinseq_barcode_no_dir(self, mock_find_fastqs, mock_os_path_exists):
        mock_find_fastqs.return_value = True
        mock_os_path_exists.return_value = False
        self.assertFalse(
            data_available_for_clinseq_barcode("test_libdir",
                                               "LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"))

    @patch('autoseq.util.clinseq_barcode.os.path.exists')
    @patch('autoseq.util.clinseq_barcode.find_fastqs')
    def test_data_available_for_clinseq_barcode_no_fastq(self, mock_find_fastqs, mock_os_path_exists):
        mock_find_fastqs.return_value = (None, None)
        mock_os_path_exists.return_value = True
        self.assertFalse(
            data_available_for_clinseq_barcode("test_libdir",
                                               "LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"))
