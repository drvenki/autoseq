import unittest
from mock import patch
from autoseq.pipeline.alascca import *
from autoseq.util.clinseq_barcode import UniqueCapture


class TestClinseq(unittest.TestCase):
    def setUp(self):
        self.sample_data_invalid = {
            "sdid": "P-NA12877",
            "T": ["AL-P-NA12877-T-03098849-TD1-TT1", "AL-P-NA12877-T-03098849-TD1-WGS"],
            "N": ["AL-P-NA12877-N-03098121-TD1-TT1", "AL-P-NA12877-N-03098121-TD1-WGS"],
            "CFDNA": ["LB-P-NA12877-CFDNA-03098850-TD1-TT1", "LB-P-NA12877-CFDNA-03098850-TD1-TT2",
                      "LB-P-NA12877-CFDNA-03098850-TD1-WGS"]
        }

        self.sample_data_valid = {
            "sdid": "P-NA12877",
            "T": ["AL-P-NA12877-T-03098849-TD1-TT1"],
            "N": ["AL-P-NA12877-N-03098121-TD1-TT1"],
            "CFDNA": []
        }

        self.ref_data = {
            "bwaIndex": "bwa/test-genome-masked.fasta",
            "chrsizes": "genome/test-genome-masked.chrsizes.txt",
            "clinvar": "variants/clinvar_20160203.vcf.gz",
            "cosmic": "variants/CosmicCodingMuts_v71.vcf.gz",
            "dbSNP": "variants/dbsnp142-germline-only.vcf.gz",
            "exac": "variants/ExAC.r0.3.1.sites.vep.vcf.gz",
            "icgc": "variants/icgc_release_20_simple_somatic_mutation.aggregated.vcf.gz",
            "reference_dict": "genome/test-genome-masked.dict",
            "reference_genome": "genome/test-genome-masked.fasta",
            "swegene_common": "variants/swegen_common.vcf.gz",
            "targets": {
                "test-regions": {
                    "cnvkit-ref": None,
                    "msisites": "intervals/targets/test-regions.msisites.tsv",
                    "targets-bed-slopped20": "intervals/targets/test-regions-GRCh37.slopped20.bed",
                    "targets-interval_list": "intervals/targets/test-regions-GRCh37.slopped20.interval_list",
                    "targets-interval_list-slopped20": "intervals/targets/test-regions-GRCh37.slopped20.interval_list"
                }
            },
            "contest_vcfs": {
                "test-regions": "test_contest.vcf"
            },
            "vep_dir": None
        }

        self.test_tumor_capture = UniqueCapture("AL", "P-NA12877", "T", "03098849", "TD", "TT")
        self.test_normal_capture = UniqueCapture("AL", "P-NA12877", "N", "03098121", "TD", "TT")

    # XXX CONTINUE HERE: TEST CONSTRUCTOR VALIDITY. THEN, WRITE TESTS FOR THE OTHER ALASCCA METHODS.

    @patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
    @patch('autoseq.pipeline.clinseq.find_fastqs')
    def test_get_normal_and_tumor_captures(self, mock_find_fastqs, mock_data_available_for_clinseq_barcode):
        mock_find_fastqs.return_value = ["test1.fq.gz", "test2.fq.gz"]
        mock_data_available_for_clinseq_barcode.return_value = True
        test_alascca_pipeline = AlasccaPipeline(self.sample_data_valid, self.ref_data, {}, "/tmp",
                                                "/nfs/LIQBIO/INBOX/exomes")
        normal_capture, tumor_capture = test_alascca_pipeline.get_normal_and_tumor_captures()
        self.assertEquals(normal_capture, self.test_normal_capture)
        self.assertEquals(tumor_capture, self.test_tumor_capture)
