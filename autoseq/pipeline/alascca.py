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
        self.configure_alascca_report_generation()

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

    def configure_alascca_report_generation(self):
        """
        Configure the Jobs for generating the ALASCCA report from the processed data.
        """

        pass

    # FIXME: BREAK OUT THIS REMAINING PANEL ANALYSIS STUFF AND PUT IT IN OTHER ALASCCA-SPECIFIC PIPELINE GENERATION METHOD(S):
    def analyze_panel(self):
        alascca_cna = AlasccaCNAPlot()
        alascca_cna.input_somatic_vcf = somatic_vcfs['vardict']
        alascca_cna.input_germline_vcf = vcfaddsample.output
        alascca_cna.input_cnr = cnvkit.output_cnr
        alascca_cna.input_cns = cnvkit.output_cns
        alascca_cna.chrsizes = self.refdata['chrsizes']
        alascca_cna.output_cna = "{}/variants/{}-alascca-cna.json".format(self.outdir,
                                                                          self.sampledata['panel']['T'])
        alascca_cna.output_purity = "{}/variants/{}-alascca-purity.json".format(self.outdir,
                                                                                self.sampledata['panel']['T'])
        alascca_cna.output_png = "{}/qc/{}-alascca-cna.png".format(self.outdir, self.sampledata['panel']['T'])
        alascca_cna.jobname = "alascca-cna/{}".format(self.sampledata['panel']['T'])
        self.add(alascca_cna)

        tlib = get_libdict(self.sampledata['panel']['T'])
        nlib = get_libdict(self.sampledata['panel']['N'])
        blood_barcode = nlib['sample_id']
        tumor_barcode = tlib['sample_id']
        metadata_json = "{}/report/{}-{}.metadata.json".format(self.outdir, blood_barcode, tumor_barcode)
        compile_metadata_json = CompileMetadata(self.referral_db_conf, blood_barcode, tumor_barcode,
                                                output_json=metadata_json,
                                                addresses=self.addresses)
        compile_metadata_json.jobname = "compile-metadata/{}-{}".format(tumor_barcode, blood_barcode)
        self.add(compile_metadata_json)

        genomic_json = "{}/report/{}-{}.genomic.json".format(self.outdir, blood_barcode, tumor_barcode)
        # FIXME: This is getting quite nasty (linking outputs and inputs, some conditional and others
        # generated in a different scope). Not sure how Daniel intended this sort of thing
        # to be done...
        contam_qc = None
        if contest_tumor.population_af_vcf:
            contam_qc = process_contest_tumor.output
        tumor_prefix = stripsuffix(os.path.basename(tbam), ".bam")
        normal_prefix = stripsuffix(os.path.basename(nbam), ".bam")
        tcov_qc = "{}/qc/{}.coverage-qc-call.json".format(self.outdir, tumor_prefix)
        ncov_qc = "{}/qc/{}.coverage-qc-call.json".format(self.outdir, normal_prefix)

        compile_genomic_json = CompileAlasccaGenomicJson(input_somatic_vcf=somatic_vcfs['vardict'],
                                                         input_cn_calls=alascca_cna.output_cna,
                                                         input_msisensor=msisensor.output,
                                                         input_purity_qc=alascca_cna.output_purity,
                                                         input_contam_qc=contam_qc,
                                                         input_tcov_qc=tcov_qc,
                                                         input_ncov_qc=ncov_qc,
                                                         output_json=genomic_json)
        compile_genomic_json.jobname = "compile-genomic/{}-{}".format(tumor_barcode, blood_barcode)
        self.add(compile_genomic_json)

        pdf = "{}/report/AlasccaReport-{}-{}.pdf".format(self.outdir, blood_barcode, tumor_barcode)
        writeAlasccaPdf = WriteAlasccaReport(input_genomic_json=compile_genomic_json.output_json,
                                             input_metadata_json=compile_metadata_json.output_json,
                                             output_pdf=pdf)
        writeAlasccaPdf.jobname = "writeAlasccaPdf/{}-{}".format(tumor_barcode, blood_barcode)
        self.add(writeAlasccaPdf)

        return {'tbam': tbam, 'nbam': nbam}

    def call_germline_variants(self, bam, library):
        """
        Call germline variants from a bam
        :param bam:
        :return:
        """
        freebayes = Freebayes()
        freebayes.input_bams = [bam]
        freebayes.somatic_only = False
        freebayes.params = None
        freebayes.reference_sequence = self.refdata['reference_genome']
        freebayes.target_bed = self.refdata['targets'][self.targets_name]['targets-bed-slopped20']
        freebayes.threads = self.maxcores
        freebayes.scratch = self.scratch
        freebayes.output = "{}/variants/{}.freebayes-germline.vcf.gz".format(self.outdir, library)
        freebayes.jobname = "freebayes-germline/{}".format(library)
        self.add(freebayes)
        return freebayes.output

    def run_fastq_qc(self, fastq_files):
        """
        Run QC on fastq files
        :param fastq_files:
        :return:
        """
        qc_files = []
        for fq in fastq_files:
            basefn = stripsuffix(os.path.basename(fq), ".fastq.gz")
            fastqc = FastQC()
            fastqc.input = fq
            fastqc.outdir = "{}/qc/fastqc/".format(self.outdir)
            fastqc.output = "{}/qc/fastqc/{}_fastqc.zip".format(self.outdir, basefn)
            fastqc.jobname = "fastqc/{}".format(basefn)
            qc_files.append(fastqc.output)
            self.add(fastqc)
        return qc_files

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

    def get_pop_af_vcf(self):
        """Get the path to the population allele frequency vcf to use in ContEst.
        Currently only available for panels CB (big design), CS (clinseq v3) & CZ (clinseq v4)"""

        tpanel = get_libdict(self.sampledata['panel']['T'])['capture_kit_name']
        npanel = get_libdict(self.sampledata['panel']['N'])['capture_kit_name']

        # TODO: Move this selection of which vcf to use somewhere else in the pipeline? Maybe possible to specify on command line?
        if tpanel == "big_design" and npanel == "big_design":
            return self.refdata['contest_vcfs']['big'] # "path/to/big_swegene_contest.vcf"
        elif tpanel in ["clinseq_v3_targets", "clinseq_v4"] and npanel in ["clinseq_v3_targets", "clinseq_v4"]:
            return self.refdata['contest_vcfs']['clinseqV3V4'] # "path/to/clinseqV3V4_exac_contest.vcf"
        elif tpanel in ["big_design", "clinseq_v3_targets", "clinseq_v4"] and npanel in ["big_design", "clinseq_v3_targets", "clinseq_v4"]:
            return self.refdata['contest_vcfs']['clinseqV3V4big'] # "path/to/clinseqV3V4big_intersection_exac_contest.vcf"
        elif tpanel == "test-regions" and npanel == "test-regions":
            return self.refdata['contest_vcfs']['test-regions'] # "path/to/test-regions_contest.vcf"
        else:
            raise ValueError("Invalid tpanel/npanel combination: {}/{}".format(tpanel, npanel))
