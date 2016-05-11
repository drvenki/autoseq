import json
import logging

from pypedream.pipeline.pypedreampipeline import PypedreamPipeline
from pypedream.tools.unix import Cat

from autoseq.tools.alignment import Bwa, SkewerSE, SkewerPE
from autoseq.tools.cnvcalling import QDNASeq, CNVkit, AlasccaCNAPlot
from autoseq.tools.intervals import MsiSensor
from autoseq.tools.picard import PicardCollectGcBiasMetrics, PicardCalculateHsMetrics, PicardCollectWgsMetrics
from autoseq.tools.picard import PicardCollectInsertSizeMetrics
from autoseq.tools.picard import PicardCollectOxoGMetrics
from autoseq.tools.qc import *
from autoseq.tools.variantcalling import Mutect2, Freebayes, VEP, VcfAddSample, VarDict
from autoseq.util.path import normpath, stripsuffix

__author__ = 'dankle'


class AlasccaPipeline(PypedreamPipeline):
    analysis_id = None
    sampledata = None
    refdata = None
    outdir = None
    maxcores = None
    scratch = "/tmp"

    def __init__(self, sampledata, refdata, outdir, analysis_id=None, maxcores=1, scratch="/tmp/", debug=False,
                 **kwargs):
        PypedreamPipeline.__init__(self, normpath(outdir), **kwargs)
        logging.debug("Unnormalized outdir is {}".format(outdir))
        logging.debug("self.outdir is {}".format(self.outdir))
        self.sampledata = sampledata
        self.refdata = refdata
        self.maxcores = maxcores
        self.analysis_id = analysis_id
        self.scratch = scratch

        panel_bams = self.analyze_panel(debug=debug)
        wgs_bams = self.analyze_lowpass_wgs()

        ################################################
        # QC
        qc_files = []

        # per-bam qc
        # panel
        all_panel_bams = [bam for bam in panel_bams.values() if bam is not None]
        qc_files += self.run_panel_bam_qc(all_panel_bams, debug=debug)
        # wgs
        all_wgs_bams = [bam for bam in wgs_bams.values() if bam is not None]
        qc_files += self.run_wgs_bam_qc(all_wgs_bams, debug=debug)

        # per-fastq qc
        fqs = self.get_all_fastqs()
        qc_files += self.run_fastq_qc(fqs)

        multiqc = MultiQC()
        multiqc.input_files = qc_files
        multiqc.search_dir = self.outdir
        multiqc.output = "{}/multiqc/{}-multiqc".format(self.outdir, self.sampledata['REPORTID'])
        multiqc.jobname = "multiqc-{}".format(self.sampledata['REPORTID'])
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
        fqs = add_item(self.sampledata['WGS_NORMAL_FQ1'], fqs)
        fqs = add_item(self.sampledata['WGS_NORMAL_FQ2'], fqs)

        return fqs

    def analyze_lowpass_wgs(self):
        if self.sampledata['WGS_TUMOR_LIB'] is None or self.sampledata['WGS_TUMOR_LIB'] == "NA":
            logging.info("No low-pass WGS data found.")
            return {'tbam': None, 'nbam': None}

        tbam = self.align_library(self.sampledata['WGS_TUMOR_FQ1'], self.sampledata['WGS_TUMOR_FQ2'],
                                  self.sampledata['WGS_TUMOR_LIB'], self.refdata['bwaIndex'], self.outdir + "/bams/wgs",
                                  maxcores=self.maxcores)
        if self.sampledata['WGS_NORMAL_FQ1']:
            nbam = self.align_library(self.sampledata['WGS_NORMAL_FQ1'], self.sampledata['WGS_NORMAL_FQ2'],
                                      self.sampledata['WGS_NORMAL_LIB'], self.refdata['bwaIndex'],
                                      self.outdir + "/bams/wgs",
                                      maxcores=self.maxcores)
        else:
            nbam = None

        qdnaseq_t = QDNASeq()
        qdnaseq_t.input = tbam
        qdnaseq_t.output_bed = self.outdir + "/cnv/qdnaseq.bed"
        qdnaseq_t.output_segments = self.outdir + "/cnv/qdnaseq.segments.txt"
        qdnaseq_t.genes_gtf = self.refdata['genesGtfGenesOnly']
        qdnaseq_t.background = self.refdata["qdnaseq_background"]
        self.add(qdnaseq_t)

        return {'tbam': tbam, 'nbam': nbam}

    def analyze_panel(self, debug=False):
        if self.sampledata['PANEL_TUMOR_LIB'] is None or self.sampledata['PANEL_TUMOR_LIB'] == "NA":
            logging.info("No panel data found.")
            return {}

        tbam = self.align_library(fq1_files=self.sampledata['PANEL_TUMOR_FQ1'],
                                  fq2_files=self.sampledata['PANEL_TUMOR_FQ2'],
                                  lib=self.sampledata['PANEL_TUMOR_LIB'],
                                  ref=self.refdata['bwaIndex'],
                                  outdir=self.outdir + "/bams/panel",
                                  maxcores=self.maxcores)

        nbam = self.align_library(fq1_files=self.sampledata['PANEL_NORMAL_FQ1'],
                                  fq2_files=self.sampledata['PANEL_NORMAL_FQ2'],
                                  lib=self.sampledata['PANEL_NORMAL_LIB'],
                                  ref=self.refdata['bwaIndex'],
                                  outdir=self.outdir + "/bams/panel",
                                  maxcores=self.maxcores)

        somatic_vcfs = self.call_somatic_variants(tbam, nbam)
        germline_vcf = self.call_germline_variants(nbam, library=self.sampledata['PANEL_NORMAL_LIB'])

        hzconcordance = HeterzygoteConcordance()
        hzconcordance.input_vcf = germline_vcf
        hzconcordance.input_bam = tbam
        hzconcordance.reference_sequence = self.refdata['reference_genome']
        hzconcordance.target_regions = self.refdata['targets'][self.sampledata['TARGETS']][
            'targets-interval_list-slopped20']
        hzconcordance.normalid = self.sampledata['PANEL_NORMAL_LIB']
        hzconcordance.filter_reads_with_N_cigar = True
        hzconcordance.jobname = "hzconcordance-{}".format(self.sampledata['PANEL_TUMOR_LIB'])
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
        vcfaddsample.jobname = "vcf-add-sample-{}".format(self.sampledata['PANEL_TUMOR_LIB'])
        self.add(vcfaddsample)

        msisensor = MsiSensor()
        msisensor.msi_sites = self.refdata['targets'][self.sampledata['TARGETS']]['msisites']
        msisensor.input_normal_bam = nbam
        msisensor.input_tumor_bam = tbam
        msisensor.output = "{}/msisensor.tsv".format(self.outdir)
        msisensor.threads = self.maxcores
        msisensor.scratch = self.scratch
        msisensor.jobname = "msisensor-{}".format(self.sampledata['PANEL_TUMOR_LIB'])
        if not debug:
            self.add(msisensor)

        # If we have a CNVkit reference
        if self.refdata['targets'][self.sampledata['TARGETS']]['cnvkit-ref']:
            cnvkit = CNVkit(input_bam=tbam,
                            reference=self.refdata['targets'][self.sampledata['TARGETS']]['cnvkit-ref'],
                            output_cnr="{}/variants/{}.cnr".format(self.outdir,
                                                                   self.sampledata['PANEL_TUMOR_LIB']),
                            output_cns="{}/variants/{}.cns".format(self.outdir,
                                                                   self.sampledata['PANEL_TUMOR_LIB']),
                            scratch=self.scratch
                            )
            cnvkit.jobname = "cnvkit-{}".format(self.sampledata['PANEL_TUMOR_LIB'])
            self.add(cnvkit)

            alascca_cna = AlasccaCNAPlot()
            alascca_cna.input_somatic_vcf = somatic_vcfs['vardict']
            alascca_cna.input_germline_vcf = germline_vcf
            alascca_cna.input_cnr = cnvkit.output_cnr
            alascca_cna.input_cns = cnvkit.output_cns
            alascca_cna.chrsizes = self.refdata['chrsizes']
            alascca_cna.output_json = "{}/variants/{}-alascca-cna.json".format(self.outdir,
                                                                               self.sampledata['PANEL_TUMOR_LIB'])
            alascca_cna.output_png = "{}/qc/{}-alascca-cna.png".format(self.outdir, self.sampledata['PANEL_TUMOR_LIB'])
            self.add(alascca_cna)

        return {'tbam': tbam, 'nbam': nbam}

    def call_germline_variants(self, bam, library):
        """
        Call germline variants from a bam
        :param bam:
        :return:
        """
        freebayes = Freebayes()
        freebayes.input_bams = [bam]
        freebayes.normalid = library
        freebayes.somatic_only = False
        freebayes.params = None
        freebayes.reference_sequence = self.refdata['reference_genome']
        freebayes.target_bed = self.refdata['targets'][self.sampledata['TARGETS']]['targets-bed-slopped20']
        freebayes.threads = self.maxcores
        freebayes.scratch = self.scratch
        freebayes.output = "{}/variants/{}.freebayes-germline.vcf.gz".format(self.outdir, library)
        freebayes.jobname = "freebayes-germline-{}".format(library)
        self.add(freebayes)
        return freebayes.output

    def call_somatic_variants(self, tbam, nbam):
        """
        Call somatic variants with freebayes and vardict.
        :param tbam: tumor bam
        :param nbam: normal bam
        :return:
        """
        mutect2 = Mutect2()
        mutect2.input_tumor = tbam
        mutect2.input_normal = nbam
        mutect2.reference_sequence = self.refdata['reference_genome']
        mutect2.target_regions = self.refdata['targets'][self.sampledata['TARGETS']]['targets-interval_list-slopped20']
        mutect2.scratch = self.scratch
        mutect2.jobname = "mutect2-{}".format(self.sampledata['PANEL_TUMOR_LIB'])
        mutect2.output = "{}/variants/{}-{}.mutect.vcf.gz".format(self.outdir, self.sampledata['PANEL_TUMOR_LIB'],
                                                                  self.sampledata['PANEL_NORMAL_LIB'])

        freebayes = Freebayes()
        freebayes.input_bams = [tbam, nbam]
        freebayes.tumorid = self.sampledata['PANEL_TUMOR_LIB']
        freebayes.normalid = self.sampledata['PANEL_NORMAL_LIB']
        freebayes.somatic_only = True
        freebayes.reference_sequence = self.refdata['reference_genome']
        freebayes.target_bed = self.refdata['targets'][self.sampledata['TARGETS']]['targets-bed-slopped20']
        freebayes.threads = self.maxcores
        freebayes.scratch = self.scratch
        freebayes.jobname = "freebayes-somatic-{}".format(self.sampledata['PANEL_TUMOR_LIB'])
        freebayes.output = "{}/variants/{}-{}.freebayes-somatic.vcf.gz".format(self.outdir,
                                                                               self.sampledata['PANEL_TUMOR_LIB'],
                                                                               self.sampledata['PANEL_NORMAL_LIB'])
        self.add(freebayes)

        vardict = VarDict(input_tumor=tbam, input_normal=nbam, tumorid=self.sampledata['PANEL_TUMOR_LIB'],
                          normalid=self.sampledata['PANEL_NORMAL_LIB'],
                          reference_sequence=self.refdata['reference_genome'],
                          target_bed=self.refdata['targets'][self.sampledata['TARGETS']]['targets-bed-slopped20'],
                          output="{}/variants/{}-{}.vardict-somatic.vcf.gz".format(self.outdir,
                                                                                   self.sampledata['PANEL_TUMOR_LIB'],
                                                                                   self.sampledata['PANEL_NORMAL_LIB']
                                                                                   )
                          )
        vardict.jobname = "vardict-{}".format(self.sampledata['PANEL_TUMOR_LIB'])
        self.add(vardict)

        vep_vardict = VEP()
        vep_vardict.input_vcf = vardict.output
        vep_vardict.threads = self.maxcores
        vep_vardict.reference_sequence = self.refdata['reference_genome']
        vep_vardict.vep_dir = self.refdata['vep_dir']
        vep_vardict.output_vcf = "{}/variants/{}-{}.vardict-somatic.vep.vcf.gz".format(self.outdir,
                                                                                       self.sampledata[
                                                                                           'PANEL_TUMOR_LIB'],
                                                                                       self.sampledata[
                                                                                           'PANEL_NORMAL_LIB'])
        vep_vardict.jobname = "vep-vardict-{}".format(self.sampledata['PANEL_TUMOR_LIB'])
        self.add(vep_vardict)

        vep_freebayes = VEP()
        vep_freebayes.input_vcf = freebayes.output
        vep_freebayes.threads = self.maxcores
        vep_freebayes.reference_sequence = self.refdata['reference_genome']
        vep_freebayes.vep_dir = self.refdata['vep_dir']
        vep_freebayes.output_vcf = "{}/variants/{}-{}.freebayes-somatic.vep.vcf.gz".format(self.outdir,
                                                                                           self.sampledata[
                                                                                               'PANEL_TUMOR_LIB'],
                                                                                           self.sampledata[
                                                                                               'PANEL_NORMAL_LIB'])
        vep_freebayes.jobname = "vep-freebayes-somatic-{}".format(self.sampledata['PANEL_TUMOR_LIB'])
        self.add(vep_freebayes)

        return {'freebayes': vep_freebayes.output_vcf, 'vardict': vep_vardict.output_vcf}

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
            fastqc.jobname = "fastqc-{}".format(basefn)
            qc_files.append(fastqc.output)
            self.add(fastqc)
        return qc_files

    def run_wgs_bam_qc(self, bams, debug=False):
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
            isize.jobname = "picard-isize-{}".format(basefn)
            isize.output_metrics = "{}/qc/picard/wgs/{}.picard-insertsize.txt".format(self.outdir, basefn)
            self.add(isize)

            if not debug:
                gcbias = PicardCollectGcBiasMetrics()
                gcbias.input = bam
                gcbias.reference_sequence = self.refdata['reference_genome']
                gcbias.output_summary = "{}/qc/picard/wgs/{}.picard-gcbias-summary.txt".format(self.outdir, basefn)
                gcbias.output_metrics = "{}/qc/picard/wgs/{}.picard-gcbias.txt".format(self.outdir, basefn)
                gcbias.jobname = "picard-gcbias-{}".format(basefn)
                gcbias.stop_after = 100
                self.add(gcbias)

            wgsmetrics = PicardCollectWgsMetrics()
            wgsmetrics.input = bam
            wgsmetrics.reference_sequence = self.refdata['reference_genome']
            wgsmetrics.output_metrics = "{}/qc/picard/wgs/{}.picard-wgsmetrics.txt".format(self.outdir, basefn)
            wgsmetrics.jobname = "picard-wgsmetrics-{}".format(basefn)
            self.add(wgsmetrics)

            qc_files += [isize.output_metrics, wgsmetrics.output_metrics]
            if not debug:
                qc_files += [gcbias.output_summary, gcbias.output_metrics]

        return qc_files

    def run_panel_bam_qc(self, bams, debug=False):
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
            isize.jobname = "picard-isize-{}".format(basefn)
            self.add(isize)

            if not debug:
                gcbias = PicardCollectGcBiasMetrics()
                gcbias.input = bam
                gcbias.reference_sequence = self.refdata['reference_genome']
                gcbias.output_summary = "{}/qc/picard/panel/{}.picard-gcbias-summary.txt".format(self.outdir, basefn)
                gcbias.output_metrics = "{}/qc/picard/panel/{}.picard-gcbias.txt".format(self.outdir, basefn)
                gcbias.jobname = "picard-gcbias-{}".format(basefn)
                gcbias.stop_after = 100
                self.add(gcbias)

            oxog = PicardCollectOxoGMetrics()
            oxog.input = bam
            oxog.reference_sequence = self.refdata['reference_genome']
            oxog.output_metrics = "{}/qc/picard/panel/{}.picard-oxog.txt".format(self.outdir, basefn)
            oxog.jobname = "picard-oxog-{}".format(basefn)
            self.add(oxog)

            hsmetrics = PicardCalculateHsMetrics()
            hsmetrics.input = bam
            hsmetrics.reference_sequence = self.refdata['reference_genome']
            hsmetrics.target_regions = self.refdata['targets'][self.sampledata['TARGETS']][
                'targets-interval_list-slopped20']
            hsmetrics.bait_regions = self.refdata['targets'][self.sampledata['TARGETS']][
                'targets-interval_list-slopped20']
            hsmetrics.bait_name = self.sampledata['TARGETS']
            hsmetrics.output_metrics = "{}/qc/picard/panel/{}.picard-hsmetrics.txt".format(self.outdir, basefn)
            hsmetrics.jobname = "picard-hsmetrics-{}".format(basefn)
            self.add(hsmetrics)

            sambamba = SambambaDepth()
            sambamba.targets_bed = self.refdata['targets'][self.sampledata['TARGETS']]['targets-bed-slopped20']
            sambamba.input = bam
            sambamba.output = "{}/qc/sambamba/{}.sambamba-depth-targets.txt".format(self.outdir, basefn)
            sambamba.jobname = "sambamba-depth-{}".format(basefn)
            self.add(sambamba)

            alascca_coverage_hist = CoverageHistogram()
            alascca_coverage_hist.input_bed = self.refdata['targets']['alascca_targets']['targets-bed-slopped20']
            alascca_coverage_hist.input_bam = bam
            alascca_coverage_hist.output = "{}/qc/{}.coverage-histogram.txt".format(self.outdir, basefn)
            alascca_coverage_hist.jobname = "alascca-coverage-hist-{}".format(basefn)
            self.add(alascca_coverage_hist)

            qc_files += [isize.output_metrics, oxog.output_metrics,
                         hsmetrics.output_metrics, sambamba.output, alascca_coverage_hist.output]
            if not debug:
                qc_files += [gcbias.output_summary, gcbias.output_metrics]

        return qc_files

    def align_library(self, fq1_files, fq2_files, lib, ref, outdir, maxcores=1):
        """
        Align fastq files for a PE library
        :param fq1_files:
        :param fq2_files:
        :param lib:
        :param ref:
        :param outdir:
        :param maxcores:
        :return:
        """
        if not fq2_files:
            logging.debug("lib {} is SE".format(lib))
            return self.align_se(fq1_files, lib, ref, outdir, maxcores)
        else:
            logging.debug("lib {} is PE".format(lib))
            return self.align_pe(fq1_files, fq2_files, lib, ref, outdir, maxcores)

    def align_se(self, fq1_files, lib, ref, outdir, maxcores):
        logging.debug("Aligning files: {}".format(fq1_files))
        fq1_abs = [normpath(x) for x in fq1_files]
        fq1_trimmed = []
        for fq1 in fq1_abs:
            skewer = SkewerSE()
            skewer.input = fq1
            skewer.output = outdir + "/skewer/{}".format(os.path.basename(fq1))
            skewer.stats = outdir + "/skewer/skewer-stats-{}.log".format(os.path.basename(fq1))
            skewer.threads = maxcores
            skewer.jobname = "skewer-{}".format(os.path.basename(fq1))
            skewer.scratch = self.scratch
            skewer.is_intermediate = True
            fq1_trimmed.append(skewer.output)
            self.add(skewer)

        cat1 = Cat()
        cat1.input = fq1_trimmed
        cat1.output = outdir + "/skewer/{}_1.fastq.gz".format(lib)
        cat1.jobname = "cat-{}".format(lib)
        cat1.is_intermediate = False
        self.add(cat1)

        bwa = Bwa()
        bwa.input_fastq1 = cat1.output
        bwa.input_reference_sequence = ref
        bwa.readgroup = "\"@RG\\tID:{lib}\\tSM:{lib}\\tLB:{lib}\\tPL:ILLUMINA\"".format(lib=lib)
        bwa.threads = maxcores
        bwa.output = "{}/{}.bam".format(outdir, lib)
        bwa.scratch = self.scratch
        bwa.jobname = "bwa-{}".format(lib)
        bwa.is_intermediate = False
        self.add(bwa)

        return bwa.output

    def align_pe(self, fq1_files, fq2_files, lib, ref, outdir, maxcores=1):
        fq1_abs = [normpath(x) for x in fq1_files]
        fq2_abs = [normpath(x) for x in fq2_files]
        logging.debug("Trimming {} and {}".format(fq1_abs, fq2_abs))
        pairs = [(fq1_abs[k], fq2_abs[k]) for k in range(len(fq1_abs))]

        fq1_trimmed = []
        fq2_trimmed = []

        for fq1, fq2 in pairs:
            skewer = SkewerPE()
            skewer.input1 = fq1
            skewer.input2 = fq2
            skewer.output1 = outdir + "/skewer/libs/{}".format(os.path.basename(fq1))
            skewer.output2 = outdir + "/skewer/libs/{}".format(os.path.basename(fq2))
            skewer.stats = outdir + "/skewer/libs/skewer-stats-{}.log".format(os.path.basename(fq1))
            skewer.threads = maxcores
            skewer.jobname = "skewer-{}".format(os.path.basename(fq1))
            skewer.scratch = self.scratch
            skewer.is_intermediate = True
            fq1_trimmed.append(skewer.output1)
            fq2_trimmed.append(skewer.output2)
            self.add(skewer)

        cat1 = Cat()
        cat1.input = fq1_trimmed
        cat1.output = outdir + "/skewer/{}-concatenated_1.fastq.gz".format(lib)
        cat1.jobname = "cat1-{}".format(lib)
        cat1.is_intermediate = True
        self.add(cat1)

        cat2 = Cat()
        cat2.input = fq2_trimmed
        cat2.jobname = "cat2-{}".format(lib)
        cat2.output = outdir + "/skewer/{}-concatenated_2.fastq.gz".format(lib)
        cat2.is_intermediate = True
        self.add(cat2)

        # cutadapt = Cutadapt()
        # cutadapt.input1 = fq1_abs
        # cutadapt.input2 = fq2_abs
        # cutadapt.output1 = outdir + "/cutadapt/cutadapt_{lib}_1.fastq.gz".format(lib=lib)
        # cutadapt.output2 = outdir + "/cutadapt/cutadapt_{lib}_2.fastq.gz".format(lib=lib)
        # #cutadapt.threads = maxcores
        # cutadapt.jobname = "cutadapt-{}".format(lib)
        # cutadapt.is_intermediate = True
        # self.add(cutadapt)

        bwa = Bwa()
        bwa.input_fastq1 = cat1.output
        bwa.input_fastq2 = cat2.output
        bwa.input_reference_sequence = ref
        bwa.readgroup = "\"@RG\\tID:{lib}\\tSM:{lib}\\tLB:{lib}\\tPL:ILLUMINA\"".format(lib=lib)
        bwa.threads = maxcores
        bwa.output = "{}/{}.bam".format(outdir, lib)
        bwa.jobname = "bwa-{}".format(lib)
        bwa.scratch = self.scratch
        bwa.is_intermediate = False
        self.add(bwa)

        return bwa.output

    def load_sampledata(self, json_file):
        with open(json_file, 'r') as f:
            self.sampledata = json.load(f)

    def load_refdata(self, json_file):
        with open(json_file, 'r') as f:
            self.refdata = json.load(f)
