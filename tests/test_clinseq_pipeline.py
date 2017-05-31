import unittest
from mock import patch
from autoseq.pipeline.clinseq import *
from autoseq.util.clinseq_barcode import UniqueCapture


class TestClinseq(unittest.TestCase):
    def setUp(self):
        self.test_single_panel_results = SinglePanelResults()
        self.test_cancer_vs_normal_results = CancerVsNormalPanelResults()
        sample_data = {
            "sdid": "P-NA12877",
            "T": ["AL-P-NA12877-T-03098849-TD1-TT1", "AL-P-NA12877-T-03098849-TD1-WGS"],
            "N": ["AL-P-NA12877-N-03098121-TD1-TT1", "AL-P-NA12877-N-03098121-TD1-WGS"],
            "CFDNA": ["LB-P-NA12877-CFDNA-03098850-TD1-TT1", "LB-P-NA12877-CFDNA-03098850-TD1-TT2",
                      "LB-P-NA12877-CFDNA-03098850-TD1-WGS"]
        }
        ref_data = {
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
        self.test_unique_capture = UniqueCapture("AL", "P-NA12877", "CFDNA", "03098850", "TD", "TT")
        self.test_clinseq_pipeline = ClinseqPipeline(sample_data, ref_data, {"cov-low-thresh-fraction": 0.8}, "/tmp", "/nfs/LIQBIO/INBOX/exomes")

    def test_single_panel_results(self):
        self.assertEquals(self.test_single_panel_results.merged_bamfile, None)

    def test_cancer_vs_normal_results(self):
        self.assertEquals(self.test_cancer_vs_normal_results.somatic_vcf, None)

    def test_pipeline_constructor(self):
        self.assertEquals(self.test_clinseq_pipeline.job_params, {'cov-low-thresh-fraction': 0.8})
        self.assertEquals(self.test_clinseq_pipeline.sampledata["sdid"], "P-NA12877")
        self.assertEquals(self.test_clinseq_pipeline.refdata["icgc"], "variants/icgc_release_20_simple_somatic_mutation.aggregated.vcf.gz")

    def test_get_job_param_set(self):
        self.assertEquals(self.test_clinseq_pipeline.get_job_param("cov-low-thresh-fraction"), 0.8)

    def test_get_job_param_default(self):
        self.assertEquals(self.test_clinseq_pipeline.get_job_param("cov-high-thresh-fold-cov"), 100)

    def test_set_germline_vcf(self):
        self.test_clinseq_pipeline.set_germline_vcf(self.test_unique_capture, "test.vcf")
        self.assertEquals(self.test_clinseq_pipeline.normal_capture_to_vcf[self.test_unique_capture], "test.vcf")

    def test_get_germline_vcf_exists(self):
        self.test_clinseq_pipeline.set_germline_vcf(self.test_unique_capture, "test.vcf")
        self.assertEquals(self.test_clinseq_pipeline.get_germline_vcf(self.test_unique_capture), "test.vcf")

    def test_get_germline_vcf_none(self):
        self.assertEquals(self.test_clinseq_pipeline.get_germline_vcf(self.test_unique_capture), None)

    def test_set_capture_bam(self):
        self.test_clinseq_pipeline.set_capture_bam(self.test_unique_capture, "test.bam")
        self.assertEquals(self.test_clinseq_pipeline.capture_to_results[self.test_unique_capture].merged_bamfile,
                          "test.bam")

    def test_set_capture_cnr(self):
        self.test_clinseq_pipeline.set_capture_cnr(self.test_unique_capture, "test.cnr")
        self.assertEquals(self.test_clinseq_pipeline.capture_to_results[self.test_unique_capture].cnr,
                          "test.cnr")

    def test_set_capture_cns(self):
        self.test_clinseq_pipeline.set_capture_cns(self.test_unique_capture, "test.cns")
        self.assertEquals(self.test_clinseq_pipeline.capture_to_results[self.test_unique_capture].cns,
                          "test.cns")

    def test_get_capture_bam_exists(self):
        self.test_clinseq_pipeline.set_capture_bam(self.test_unique_capture, "test.bam")
        self.assertEquals(self.test_clinseq_pipeline.get_capture_bam(self.test_unique_capture),
                          "test.bam")

    def test_get_capture_bam_none(self):
        self.assertEquals(self.test_clinseq_pipeline.get_capture_bam(self.test_unique_capture),
                          None)

    @patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
    def test_check_sampledata_all_available(self, mock_data_available_for_clinseq_barcode):
        mock_data_available_for_clinseq_barcode.return_value = True
        # Test that no changes occur to the sample data if the data is all available:
        self.test_clinseq_pipeline.check_sampledata()
        self.assertEquals(self.test_clinseq_pipeline.sampledata,
                          {
                              "sdid": "P-NA12877",
                              "T": ["AL-P-NA12877-T-03098849-TD1-TT1", "AL-P-NA12877-T-03098849-TD1-WGS"],
                              "N": ["AL-P-NA12877-N-03098121-TD1-TT1", "AL-P-NA12877-N-03098121-TD1-WGS"],
                              "CFDNA": ["LB-P-NA12877-CFDNA-03098850-TD1-TT1", "LB-P-NA12877-CFDNA-03098850-TD1-TT2",
                                        "LB-P-NA12877-CFDNA-03098850-TD1-WGS"]
                          })

    @patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
    def test_check_sampledata_none_available(self, mock_data_available_for_clinseq_barcode):
        mock_data_available_for_clinseq_barcode.return_value = False
        # Test that no changes occur to the sample data if the data is all available:
        self.test_clinseq_pipeline.check_sampledata()
        self.assertEquals(self.test_clinseq_pipeline.sampledata,
                          {
                              "sdid": "P-NA12877",
                              "T": [],
                              "N": [],
                              "CFDNA": []
                          })

    def test_vep_is_set(self):
        self.assertEquals(self.test_clinseq_pipeline.vep_is_set(), False)

    def test_get_all_unique_capture(self):
        self.assertEquals(self.test_clinseq_pipeline.get_all_unique_captures(), [])
