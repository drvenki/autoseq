import unittest

from autoseq.util.orderform import *


class TestOrderform(unittest.TestCase):
    def test_parse_orderform_block_valid(self):
        fields_to_parse = ["<SAMPLE ENTRIES>",
                           "LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000",
                           "</SAMPLE ENTRIES>"]
        parsed_clinseq_barcodes = parse_orderform_block(fields_to_parse)
        self.assertEquals(parsed_clinseq_barcodes,
                          ["LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"])

    def test_parse_orderform_block_empty(self):
        fields_to_parse = ["<SAMPLE ENTRIES>",
                           "</SAMPLE ENTRIES>"]
        parsed_clinseq_barcodes = parse_orderform_block(fields_to_parse)
        self.assertEquals(parsed_clinseq_barcodes, [])

    def test_parse_orderform_block_no_start(self):
        fields_to_parse = ["LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000",
                           "</SAMPLE ENTRIES>"]
        parsed_clinseq_barcodes = parse_orderform_block(fields_to_parse)
        self.assertEquals(parsed_clinseq_barcodes, [])

    def test_parse_orderform_block_no_end(self):
        fields_to_parse = ["<SAMPLE ENTRIES>",
                           "LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"]
        parsed_clinseq_barcodes = parse_orderform_block(fields_to_parse)
        self.assertEquals(parsed_clinseq_barcodes,
                          ["LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"])
