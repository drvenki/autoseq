import logging

from autoseq.pipeline.clinseq import ClinseqPipeline
from autoseq.tools.alignment import align_library
from autoseq.tools.cnvcalling import AlasccaCNAPlot
from autoseq.tools.intervals import MsiSensor
from autoseq.tools.picard import PicardCollectHsMetrics
from autoseq.tools.picard import PicardCollectInsertSizeMetrics
from autoseq.tools.picard import PicardCollectOxoGMetrics
from autoseq.tools.qc import *
from autoseq.tools.reports import CompileMetadata, CompileAlasccaGenomicJson, WriteAlasccaReport
from autoseq.tools.variantcalling import Freebayes, VcfAddSample, call_somatic_variants
from autoseq.util.library import get_libdict, find_fastqs
from autoseq.util.path import normpath, stripsuffix
from autoseq.pipeline.clinseq import compose_sample_str

__author__ = 'dankle'


class AlasccaPipeline(ClinseqPipeline):
    def __init__(self, sampledata, refdata, outdir, libdir, analysis_id=None, maxcores=1, scratch="/scratch/tmp/tmp/",
                 referral_db_conf="tests/referrals/referral-db-config.json",
                 addresses="tests/referrals/addresses.csv",
                 **kwargs):
        ClinseqPipeline.__init__(self, sampledata, refdata, outdir, libdir, analysis_id,
                                 maxcores, scratch, **kwargs)

        self.referral_db_conf = referral_db_conf
        self.addresses = addresses

        # Check to ensure that the sample data is valid for an ALASCCA analysis:
        self.validate_sample_data_for_alascca()

        # Remove sample capture items for which data is not available:
        self.check_sampledata()

        # Configure all panel analyses:
        self.configure_panel_analyses()

        # Configure ALASCCA report generation:
        self.configure_alascca_specific_analysis()

        # Configure QC of all panel data:
        self.configure_all_panel_qcs()

        # Configure fastq QCs:
        self.configure_fastq_qcs()

        # Configure MultiQC:
        self.configure_multi_qc()

    def validate_sample_data_for_alascca(self):
        """
        Checks validity of the sample data. Raises a ValueError if the sampledata dictionary does
        not fit into the expected ALASCCA analysis pipeline limitations.
        """

        # FIXME: Implement this
        pass

    def get_normal_and_tumor_captures(self):
        """
        Retrieves the unique normal and tumor capture identifiers for this ALASCCA analysis.
        :return: (normal_capture, tumor_capture) tuple, denoting those unique library captures.
        """

        # There must be exactly one tumor and exactly one normal for this to be valid:
        if len(self.get_unique_cancer_captures()) != 1 or \
           len(self.get_unique_normal_captures()) != 1:
            raise ValueError("Invalid pipeline state for configuration of ALASCCA CNA.")

        normal_capture = self.get_unique_normal_captures()[0]
        tumor_capture = self.get_unique_normal_captures()[0]

        return normal_capture, tumor_capture

    def configure_alascca_cna(self, normal_capture, tumor_capture):
        tumor_vs_normal_results = self.normal_cancer_pair_to_results[(normal_capture, tumor_capture)]
        tumor_results = self.capture_to_results[tumor_capture]

        alascca_cna = AlasccaCNAPlot()
        alascca_cna.input_somatic_vcf = tumor_vs_normal_results.somatic_vcf
        alascca_cna.input_germline_vcf = tumor_vs_normal_results.vcf_addsample_output
        alascca_cna.input_cnr = tumor_results.cnr
        alascca_cna.input_cns = tumor_results.cna
        alascca_cna.chrsizes = self.refdata['chrsizes']

        tumor_str = compose_sample_str(tumor_capture)

        alascca_cna.output_cna = "{}/variants/{}-alascca-cna.json".format(
            self.outdir, tumor_str)
        alascca_cna.output_purity = "{}/variants/{}-alascca-purity.json".format(
            self.outdir, tumor_str)
        alascca_cna.output_png = "{}/qc/{}-alascca-cna.png".format(
            self.outdir, tumor_str)
        alascca_cna.jobname = "alascca-cna/{}".format(tumor_str)
        self.add(alascca_cna)
        return alascca_cna.output_cna

    def configure_compile_metadata(self, normal_capture, tumor_capture):
        blood_barcode = normal_capture.sample_id
        tumor_barcode = tumor_capture.sample_id
        metadata_json = "{}/report/{}-{}.metadata.json".format(self.outdir, blood_barcode, tumor_barcode)
        compile_metadata_json = CompileMetadata(self.referral_db_conf, blood_barcode, tumor_barcode,
                                                output_json = metadata_json,
                                                addresses = self.addresses)
        compile_metadata_json.jobname = "compile-metadata/{}-{}".format(tumor_barcode, blood_barcode)
        self.add(compile_metadata_json)
        return compile_metadata_json.output_json

    def configure_compile_genomic_json(self, normal_capture, tumor_capture,
                                       alascca_cna_output, alascca_cna_purity_call):
        tumor_vs_normal_results = self.normal_cancer_pair_to_results[(normal_capture, tumor_capture)]

        compile_genomic_json = CompileAlasccaGenomicJson(
            input_somatic_vcf=tumor_vs_normal_results.somatic_vcf,
            input_cn_calls=alascca_cna_output,
            input_msisensor=tumor_vs_normal_results.msi_output,
            input_purity_qc=alascca_cna_purity_call,
            input_contam_qc=tumor_vs_normal_results.cancer_contam_call,
            input_tcov_qc=tcov_qc,
            input_ncov_qc=ncov_qc,
            output_json=genomic_json)

        compile_genomic_json.jobname = "compile-genomic/{}-{}".format(tumor_barcode, blood_barcode)
        self.add(compile_genomic_json)

    def configure_write_alascca_report(self, metadata_json, genomic_json):
        pdf = "{}/report/AlasccaReport-{}-{}.pdf".format(self.outdir, blood_barcode, tumor_barcode)
        writeAlasccaPdf = WriteAlasccaReport(input_genomic_json=compile_genomic_json.output_json,
                                             input_metadata_json=compile_metadata_json.output_json,
                                             output_pdf=pdf)
        writeAlasccaPdf.jobname = "writeAlasccaPdf/{}-{}".format(tumor_barcode, blood_barcode)
        self.add(writeAlasccaPdf)

    def configure_alascca_report_generation(self, normal_capture, tumor_capture, alascca_cna_output):
        """
        Configure the generation of the ALASCCA report for this pipeline instance.
        """

        metadata_json = self.configure_compile_metadata(normal_capture, tumor_capture)
        genomic_json = self.configure_compile_genomic_json(
            normal_capture, tumor_capture, alascca_cna_output)
        self.configure_write_alascca_report(metadata_json, genomic_json)

    def configure_alascca_specific_analysis(self):
        """
        Configure the Jobs for specific to the ALASCCA pipeline.
        """

        normal_capture, tumor_capture = self.get_normal_and_tumor_captures()
        alascca_cna_ouput = self.configure_alascca_cna(normal_capture, tumor_capture)
        self.configure_alascca_report_generation(normal_capture, tumor_capture, alascca_cna_ouput)

    def run_panel_bam_qc(self, bams):
        """
        Run QC on panel bams
        :param bams: list of bams
        :return: list of generated files
        """

        qc_files = []
        for bam in bams:
            logging.debug("Adding QC jobs for {}".format(bam))
            basefn = stripsuffix(os.path.basename(bam), ".bam")
            isize = PicardCollectInsertSizeMetrics()
            isize.input = bam
            isize.output_metrics = "{}/qc/picard/panel/{}.picard-insertsize.txt".format(self.outdir, basefn)
            isize.jobname = "picard-isize/{}".format(basefn)
            self.add(isize)

            oxog = PicardCollectOxoGMetrics()
            oxog.input = bam
            oxog.reference_sequence = self.refdata['reference_genome']
            oxog.output_metrics = "{}/qc/picard/panel/{}.picard-oxog.txt".format(self.outdir, basefn)
            oxog.jobname = "picard-oxog/{}".format(basefn)
            self.add(oxog)

            hsmetrics = PicardCollectHsMetrics()
            hsmetrics.input = bam
            hsmetrics.reference_sequence = self.refdata['reference_genome']
            hsmetrics.target_regions = self.refdata['targets'][self.targets_name][
                'targets-interval_list-slopped20']
            hsmetrics.bait_regions = self.refdata['targets'][self.targets_name][
                'targets-interval_list-slopped20']
            hsmetrics.bait_name = self.targets_name
            hsmetrics.output_metrics = "{}/qc/picard/panel/{}.picard-hsmetrics.txt".format(self.outdir, basefn)
            hsmetrics.jobname = "picard-hsmetrics/{}".format(basefn)
            self.add(hsmetrics)

            sambamba = SambambaDepth()
            sambamba.targets_bed = self.refdata['targets'][self.targets_name]['targets-bed-slopped20']
            sambamba.input = bam
            sambamba.output = "{}/qc/sambamba/{}.sambamba-depth-targets.txt".format(self.outdir, basefn)
            sambamba.jobname = "sambamba-depth/{}".format(basefn)
            self.add(sambamba)

            alascca_coverage_hist = CoverageHistogram()
            if 'alascca_targets' in self.refdata['targets']:
                alascca_coverage_hist.input_bed = self.refdata['targets']['alascca_targets']['targets-bed-slopped20']
            else:
                alascca_coverage_hist.input_bed = self.refdata['targets'][self.targets_name][
                    'targets-bed-slopped20']
            alascca_coverage_hist.input_bam = bam
            alascca_coverage_hist.output = "{}/qc/{}.coverage-histogram.txt".format(self.outdir, basefn)
            alascca_coverage_hist.jobname = "alascca-coverage-hist/{}".format(basefn)
            self.add(alascca_coverage_hist)

            coverage_qc_call = CoverageCaveat()
            coverage_qc_call.input_histogram = alascca_coverage_hist.output
            coverage_qc_call.output = "{}/qc/{}.coverage-qc-call.json".format(self.outdir, basefn)
            coverage_qc_call.jobname = "coverage-qc-call/{}".format(basefn)
            self.add(coverage_qc_call)

            qc_files += [isize.output_metrics, oxog.output_metrics,
                         hsmetrics.output_metrics, sambamba.output, alascca_coverage_hist.output,
                         coverage_qc_call.output]

        return qc_files
