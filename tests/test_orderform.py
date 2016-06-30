import unittest

from autoseq.cli.liqbio import make_sample_dicts
from autoseq.util.orderform import parse_orderform


class TestOrderform(unittest.TestCase):
    orderform = 'tests/liqbio_test_orderform.xlsx'

    def test_parse_libraries(self):
        """
        test that an orderform is parsed correctly
        """

        libs = parse_orderform(self.orderform)
        # first item in orderform
        first_item = {'library_id': 'NA12877-T-03098849-TD1-TT1',
                      'capture_kit_name': 'test-regions',
                      'capture_id': 'TT1',
                      'type': 'T',
                      'sample_id': '03098849',
                      'prep_kit_name': 'THRUPLEX_DNASEQ',
                      'sdid': 'NA12877',
                      'prep_id': 'TD1',
                      'project_id': None,
                      'project_name': None}
        self.assertIn(first_item, libs,
                      "NA12877-T-03098849-TD1-TT1 could not be parsed from orderform")

        last_item = {'library_id': 'NA12877-P-03098850-TD1-WGS',
                     'capture_kit_name': 'lowpass_wgs',
                     'capture_id': 'WGS',
                     'type': 'P',
                     'sample_id': '03098850',
                     'prep_kit_name': 'THRUPLEX_DNASEQ',
                     'sdid': 'NA12877',
                     'prep_id': 'TD1',
                     'project_id': None,
                     'project_name': None
                     }
        # last item in orderform
        self.assertIn(last_item, libs,
                      "NA12877-P-03098850-TD1-WGS could not be parsed from orderform")

    def test_make_sample_dict(self):
        libs = parse_orderform(self.orderform)
        sample_dicts = make_sample_dicts(libs)

        self.assertIn('NA12877', sample_dicts)

        sample_dict_NA12877 = {
            "sdid": "NA12877",
            "panel": {
                "T": "NA12877-T-03098849-TD1-TT1",
                "N": "NA12877-N-03098121-TD1-TT1",
                "P": ["NA12877-P-03098850-TD1-TT1", "NA12877-P-03098850-TD1-TT2"]
            },
            "wgs": {
                "T": "NA12877-T-03098849-TD1-WGS",
                "N": "NA12877-N-03098121-TD1-WGS",
                "P": ["NA12877-P-03098850-TD1-WGS"]
            }
        }

        self.assertEqual(sample_dict_NA12877, sample_dicts['NA12877'])
