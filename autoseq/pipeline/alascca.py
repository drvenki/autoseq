import logging

from pypedream.pipeline.pypedreampipeline import PypedreamPipeline

from autoseq.tools.alignment import align_library
from autoseq.tools.cnvcalling import QDNASeq, CNVkit, AlasccaCNAPlot
from autoseq.tools.intervals import MsiSensor
from autoseq.tools.picard import PicardCollectGcBiasMetrics, PicardCollectWgsMetrics, \
    PicardCollectHsMetrics
from autoseq.tools.picard import PicardCollectInsertSizeMetrics
from autoseq.tools.picard import PicardCollectOxoGMetrics
from autoseq.tools.qc import *
from autoseq.tools.reports import CompileMetadata, CompileAlasccaGenomicJson, WriteAlasccaReport
from autoseq.tools.variantcalling import Freebayes, VcfAddSample, call_somatic_variants
from autoseq.util.library import get_libdict
from autoseq.util.path import normpath, stripsuffix

__author__ = 'dankle'


class AlasccaPipeline(PypedreamPipeline):

    def __init__(self, sampledata, refdata, outdir, libdir, analysis_id=None, maxcores=1, scratch="/tmp/",
                 referral_db_conf="tests/referrals/referral-db-config.json",
                 addresses="tests/referrals/addresses.csv",
                 **kwargs):
        PypedreamPipeline.__init__(self, normpath(outdir), **kwargs)
        logging.debug("Unnormalized outdir is {}".format(outdir))
        logging.debug("self.outdir is {}".format(self.outdir))
        self.libdir = libdir
        self.sampledata = sampledata
        self.refdata = refdata
        self.maxcores = maxcores
        self.analysis_id = analysis_id
        self.scratch = scratch
        self.referral_db_conf = referral_db_conf
        self.addresses = addresses
        self.targets_name = get_libdict(self.sampledata['panel']['T'])['capture_kit_name']


        panel_bams = self.analyze_panel()
        wgs_bams = self.analyze_lowpass_wgs()

        ################################################
        # QC
        qc_files = []

        # per-bam qc
        # panel
        all_panel_bams = [bam for bam in panel_bams.values() if bam is not None]
        qc_files += self.run_panel_bam_qc(all_panel_bams)
        # wgs
        all_wgs_bams = [bam for bam in wgs_bams.values() if bam is not None]
        qc_files += self.run_wgs_bam_qc(all_wgs_bams)

        # per-fastq qc
        fqs = self.get_all_fastqs()
        qc_files += self.run_fastq_qc(fqs)

        multiqc = MultiQC()
        multiqc.input_files = qc_files
        multiqc.search_dir = self.outdir
        multiqc.output = "{}/multiqc/{}-multiqc".format(self.outdir, self.sampledata['REPORTID'])
        multiqc.jobname = "multiqc/{}".format(self.sampledata['REPORTID'])
        self.add(multiqc)

    def get_all_fastqs(self):
        fqs = []

        def add_item(item, lst):
            if item != [] and item is not None:
                return lst + item
            else:
                return lst

        fqs = add_item(self.sampledata['PANEL_TUMOR_FQ1'], fqs)
        fqs = add_item(self.sampledata['PANEL_TUMOR_FQ2'], fqs)
        fqs = add_item(self.sampledata['PANEL_NORMAL_FQ1'], fqs)
        fqs = add_item(self.sampledata['PANEL_NORMAL_FQ2'], fqs)
        fqs = add_item(self.sampledata['WGS_TUMOR_FQ1'], fqs)
        fqs = add_item(self.sampledata['WGS_TUMOR_FQ2'], fqs)

        return fqs

    def analyze_lowpass_wgs(self):
        if self.sampledata['WGS_TUMOR_LIB'] is None or self.sampledata['WGS_TUMOR_LIB'] == "NA":
            logging.info("No low-pass WGS data found.")
            return {'tbam': None, 'nbam': None}

        tbam = align_library(self,
                             self.sampledata['WGS_TUMOR_FQ1'],
                             self.sampledata['WGS_TUMOR_FQ2'],
                             self.sampledata['WGS_TUMOR_LIB'],
                             self.refdata['bwaIndex'],
                             self.outdir + "/bams/wgs",
                             maxcores=self.maxcores)

        qdnaseq_t = QDNASeq(tbam,
                            output_segments=os.path.join(self.outdir, "cnv", "{}-qdnaseq.segments.txt".format(
                                self.sampledata['WGS_TUMOR_LIB'])),
                            background=None)
        self.add(qdnaseq_t)

        return {'tbam': tbam}

    def analyze_panel(self):
        if self.sampledata['PANEL_TUMOR_LIB'] is None or self.sampledata['PANEL_TUMOR_LIB'] == "NA":
            logging.info("No panel data found.")
            return {}

        tbam = align_library(self,
                             fq1_files=self.sampledata['PANEL_TUMOR_FQ1'],
                             fq2_files=self.sampledata['PANEL_TUMOR_FQ2'],
                             lib=self.sampledata['PANEL_TUMOR_LIB'],
                             ref=self.refdata['bwaIndex'],
                             outdir=self.outdir + "/bams/panel",
                             maxcores=self.maxcores)

        nbam = align_library(self,
                             fq1_files=self.sampledata['PANEL_NORMAL_FQ1'],
                             fq2_files=self.sampledata['PANEL_NORMAL_FQ2'],
                             lib=self.sampledata['PANEL_NORMAL_LIB'],
                             ref=self.refdata['bwaIndex'],
                             outdir=self.outdir + "/bams/panel",
                             maxcores=self.maxcores)

        vep = False
        if self.refdata['vep_dir']:
            vep = True

        somatic_vcfs = call_somatic_variants(self, tbam, nbam,
                                             tlib=self.sampledata['PANEL_TUMOR_LIB'],
                                             nlib=self.sampledata['PANEL_NORMAL_LIB'],
                                             target_name=self.targets,
                                             refdata=self.refdata,
                                             outdir=self.outdir,
                                             callers=['vardict'],
                                             vep=vep)

        germline_vcf = self.call_germline_variants(nbam, library=self.sampledata['PANEL_NORMAL_LIB'])

        libdict = get_libdict(self.sampledata['PANEL_NORMAL_LIB'])
        rg_sm = "{}-{}-{}".format(libdict['sdid'], libdict['type'], libdict['sample_id'])

        hzconcordance = HeterzygoteConcordance()
        hzconcordance.input_vcf = germline_vcf
        hzconcordance.input_bam = tbam
        hzconcordance.reference_sequence = self.refdata['reference_genome']
        hzconcordance.target_regions = self.refdata['targets'][self.targets_name]['targets-interval_list-slopped20']
        hzconcordance.normalid = rg_sm
        hzconcordance.filter_reads_with_N_cigar = True
        hzconcordance.jobname = "hzconcordance/{}".format(self.sampledata['PANEL_TUMOR_LIB'])
        hzconcordance.output = "{}/bams/{}-{}-hzconcordance.txt".format(self.outdir, self.sampledata['PANEL_TUMOR_LIB'],
                                                                        self.sampledata['PANEL_NORMAL_LIB'])
        self.add(hzconcordance)

        vcfaddsample = VcfAddSample()
        vcfaddsample.input_bam = tbam
        vcfaddsample.input_vcf = germline_vcf
        vcfaddsample.samplename = self.sampledata['PANEL_TUMOR_LIB']
        vcfaddsample.filter_hom = True
        vcfaddsample.scratch = self.scratch
        vcfaddsample.output = "{}/variants/{}-and-{}.germline.vcf.gz".format(self.outdir,
                                                                             self.sampledata['PANEL_TUMOR_LIB'],
                                                                             self.sampledata['PANEL_NORMAL_LIB'])
        vcfaddsample.jobname = "vcf-add-sample/{}".format(self.sampledata['PANEL_TUMOR_LIB'])
        self.add(vcfaddsample)

        msisensor = MsiSensor()
        msisensor.msi_sites = self.refdata['targets'][self.targets]['msisites']
        msisensor.input_normal_bam = nbam
        msisensor.input_tumor_bam = tbam
        msisensor.output = "{}/msisensor.tsv".format(self.outdir)
        msisensor.threads = self.maxcores
        msisensor.scratch = self.scratch
        msisensor.jobname = "msisensor/{}".format(self.sampledata['PANEL_TUMOR_LIB'])
        self.add(msisensor)

        cnvkit = CNVkit(input_bam=tbam,
                        output_cnr="{}/cnv/{}.cnr".format(self.outdir,
                                                          self.sampledata['PANEL_TUMOR_LIB']),
                        output_cns="{}/cnv/{}.cns".format(self.outdir,
                                                          self.sampledata['PANEL_TUMOR_LIB']),
                        scratch=self.scratch
                        )
        # If we have a CNVkit reference
        if self.refdata['targets'][self.targets]['cnvkit-ref']:
            cnvkit.reference = self.refdata['targets'][self.targets]['cnvkit-ref']
        else:
            cnvkit.targets_bed = self.refdata['targets'][self.targets]['targets-bed-slopped20']

        cnvkit.jobname = "cnvkit/{}".format(self.sampledata['PANEL_TUMOR_LIB'])
        self.add(cnvkit)

        alascca_cna = AlasccaCNAPlot()
        alascca_cna.input_somatic_vcf = somatic_vcfs['vardict']
        alascca_cna.input_germline_vcf = vcfaddsample.output
        alascca_cna.input_cnr = cnvkit.output_cnr
        alascca_cna.input_cns = cnvkit.output_cns
        alascca_cna.chrsizes = self.refdata['chrsizes']
        alascca_cna.output_json = "{}/variants/{}-alascca-cna.json".format(self.outdir,
                                                                           self.sampledata['PANEL_TUMOR_LIB'])
        alascca_cna.output_png = "{}/qc/{}-alascca-cna.png".format(self.outdir, self.sampledata['PANEL_TUMOR_LIB'])
        self.add(alascca_cna)

        tlib = get_libdict(self.sampledata['PANEL_TUMOR_LIB'])
        nlib = get_libdict(self.sampledata['PANEL_NORMAL_LIB'])
        blood_barcode = nlib['sample_id']
        tumor_barcode = tlib['sample_id']
        metadata_json = "{}/report/{}-{}.metadata.json".format(self.outdir, blood_barcode, tumor_barcode)
        compile_metadata_json = CompileMetadata(self.referral_db_conf, blood_barcode, tumor_barcode,
                                                output_json=metadata_json,
                                                addresses=self.addresses)
        compile_metadata_json.jobname = "compile-metadata/{}-{}".format(tumor_barcode, blood_barcode)
        self.add(compile_metadata_json)

        genomic_json = "{}/report/{}-{}.genomic.json".format(self.outdir, blood_barcode, tumor_barcode)
        compile_genomic_json = CompileAlasccaGenomicJson(input_somatic_vcf=somatic_vcfs['vardict'],
                                                         input_cn_calls=alascca_cna.output_json,
                                                         input_msisensor=msisensor.output,
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
        freebayes.target_bed = self.refdata['targets'][self.targets]['targets-bed-slopped20']
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

    def run_wgs_bam_qc(self, bams):
        """
        Run QC on wgs bams
        :param bams: list of bams
        :return: list of generated files
        """
        qc_files = []
        logging.debug("bams are {}".format(bams))
        for bam in bams:
            basefn = stripsuffix(os.path.basename(bam), ".bam")
            isize = PicardCollectInsertSizeMetrics()
            isize.input = bam
            isize.jobname = "picard-isize/{}".format(basefn)
            isize.output_metrics = "{}/qc/picard/wgs/{}.picard-insertsize.txt".format(self.outdir, basefn)
            self.add(isize)

            qc_files += [isize.output_metrics]

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
            hsmetrics.target_regions = self.refdata['targets'][self.targets][
                'targets-interval_list-slopped20']
            hsmetrics.bait_regions = self.refdata['targets'][self.targets][
                'targets-interval_list-slopped20']
            hsmetrics.bait_name = self.targets
            hsmetrics.output_metrics = "{}/qc/picard/panel/{}.picard-hsmetrics.txt".format(self.outdir, basefn)
            hsmetrics.jobname = "picard-hsmetrics/{}".format(basefn)
            self.add(hsmetrics)

            sambamba = SambambaDepth()
            sambamba.targets_bed = self.refdata['targets'][self.targets]['targets-bed-slopped20']
            sambamba.input = bam
            sambamba.output = "{}/qc/sambamba/{}.sambamba-depth-targets.txt".format(self.outdir, basefn)
            sambamba.jobname = "sambamba-depth/{}".format(basefn)
            self.add(sambamba)

            alascca_coverage_hist = CoverageHistogram()
            if 'alascca_targets' in self.refdata['targets']:
                alascca_coverage_hist.input_bed = self.refdata['targets']['alascca_targets']['targets-bed-slopped20']
            else:
                alascca_coverage_hist.input_bed = self.refdata['targets'][self.targets][
                    'targets-bed-slopped20']
            alascca_coverage_hist.input_bam = bam
            alascca_coverage_hist.output = "{}/qc/{}.coverage-histogram.txt".format(self.outdir, basefn)
            alascca_coverage_hist.jobname = "alascca-coverage-hist/{}".format(basefn)
            self.add(alascca_coverage_hist)

            qc_files += [isize.output_metrics, oxog.output_metrics,
                         hsmetrics.output_metrics, sambamba.output, alascca_coverage_hist.output]

        return qc_files
