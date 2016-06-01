import json
import logging

from pypedream.pipeline.pypedreampipeline import PypedreamPipeline
from pypedream.tools.unix import Cat

from autoseq.tools.alignment import Bwa, SkewerSE, SkewerPE
from autoseq.tools.cnvcalling import QDNASeq
from autoseq.tools.intervals import MsiSensor
from autoseq.tools.picard import PicardCollectGcBiasMetrics, PicardCollectWgsMetrics, PicardCollectHsMetrics
from autoseq.tools.picard import PicardCollectInsertSizeMetrics
from autoseq.tools.picard import PicardCollectOxoGMetrics
from autoseq.tools.qc import *
from autoseq.tools.variantcalling import Mutect2, Freebayes, VEP, VcfAddSample, VarDict
from autoseq.util.library import get_libdict
from autoseq.util.path import normpath, stripsuffix

__author__ = 'dankle'


class HopkinsMappingPipeline(PypedreamPipeline):
    analysis_id = None
    sampledata = None
    refdata = None
    outdir = None
    maxcores = None

    def __init__(self, sampledata, refdata, outdir, analysis_id=None, maxcores=1, debug=False, **kwargs):
        PypedreamPipeline.__init__(self, normpath(outdir), **kwargs)
        self.sampledata = sampledata
        self.refdata = refdata
        self.maxcores = maxcores
        self.analysis_id = analysis_id

        bams = []
        fqs = []
        with open(self.sampledata, 'r') as f:
            for ln in f:
                if "COUNT" in ln:
                    continue
                parts = ln.strip().split("\t")
                if len(parts) < 9:
                    continue
                    
                lib = parts[3]
                fq1 = parts[9]
                fq2 = parts[10]
                bam = self.align_pe(fq1=fq1, fq2=fq2, lib=lib, ref=self.refdata['bwaIndex'], outdir=self.outdir)
                bams.append(bam)
                fqs.append(fq1)
                fqs.append(fq2)

        ################################################
        # QC
        qc_files = []

        # per-bam qc
        # panel
        logging.debug("Bam files are {}".format(bams))
        qc_files += self.run_panel_bam_qc(bams, debug=debug)

        # per-fastq qc
        logging.debug("fqs = {}".format(fqs))
        qc_files += self.run_fastq_qc(fqs)

        multiqc = MultiQC()
        multiqc.input_files = qc_files
        multiqc.search_dir = self.outdir
        multiqc.output = "{}/multiqc/{}-multiqc".format(self.outdir, self.sampledata['sdid'])
        multiqc.jobname = "multiqc-{}".format(self.sampledata['sdid'])
        self.add(multiqc)

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

    def run_panel_bam_qc(self, bams, debug=False):
        """
        Run QC on panel bams
        :param bams: list of bams
        :return: list of generated files
        """

        qc_files = []
        for bam in bams:
            lib = stripsuffix(os.path.basename(bam), ".bam")
            targets = 'clinseq_v1_targets'
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

            hsmetrics = PicardCollectHsMetrics()
            hsmetrics.input = bam
            hsmetrics.reference_sequence = self.refdata['reference_genome']
            hsmetrics.target_regions = self.refdata['targets'][targets][
                'targets-interval_list-slopped20']
            hsmetrics.bait_regions = self.refdata['targets'][targets][
                'targets-interval_list-slopped20']
            hsmetrics.bait_name = targets
            hsmetrics.output_metrics = "{}/qc/picard/panel/{}.picard-hsmetrics.txt".format(self.outdir, basefn)
            hsmetrics.jobname = "picard-hsmetrics-{}".format(basefn)
            self.add(hsmetrics)

            sambamba = SambambaDepth()
            sambamba.targets_bed = self.refdata['targets'][targets]['targets-bed-slopped20']
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

    def align_pe(self, fq1, fq2, lib, ref, outdir):

        skewer = SkewerPE()
        skewer.input1 = fq1
        skewer.input2 = fq2
        skewer.output1 = outdir + "/skewer/libs/{}".format(os.path.basename(fq1))
        skewer.output2 = outdir + "/skewer/libs/{}".format(os.path.basename(fq2))
        skewer.stats = outdir + "/skewer/libs/skewer-stats-{}.log".format(os.path.basename(fq1))
        skewer.threads = self.maxcores
        skewer.jobname = "skewer-{}".format(os.path.basename(fq1))
        skewer.is_intermediate = True
        self.add(skewer)

        bwa = Bwa()
        bwa.input_fastq1 = skewer.output1
        bwa.input_fastq2 = skewer.output2
        bwa.input_reference_sequence = ref
        bwa.readgroup = "\"@RG\\tID:{lib}\\tSM:{lib}\\tLB:{lib}\\tPL:ILLUMINA\"".format(lib=lib)
        bwa.threads = self.maxcores
        bwa.output = "{}/{}.bam".format(outdir, lib)
        bwa.jobname = "bwa-{}".format(lib)
        bwa.is_intermediate = False
        self.add(bwa)

        return bwa.output
