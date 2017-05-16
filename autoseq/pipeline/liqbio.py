import logging

from autoseq.pipeline.clinseq import ClinseqPipeline
from autoseq.tools.alignment import align_library
from autoseq.tools.cnvcalling import QDNASeq
from autoseq.tools.picard import PicardCollectWgsMetrics
from autoseq.tools.picard import PicardCollectInsertSizeMetrics
from autoseq.tools.qc import *
from autoseq.util.library import find_fastqs
from autoseq.util.path import stripsuffix

__author__ = 'dankle'


class LiqBioPipeline(ClinseqPipeline):
    def __init__(self, sampledata, refdata, outdir, libdir, maxcores=1, scratch="/scratch/tmp/tmp",
                 **kwargs):
        ClinseqPipeline.__init__(self, sampledata, refdata, outdir, libdir,
                                 maxcores, scratch, **kwargs)

        # Remove clinseq barcodes for which data is not available:
        self.check_sampledata()

        # Configure alignment and merging of fastq data for all clinseq barcodes:
        self.configure_align_and_merge()

        # Configure all panel analyses:
        self.configure_panel_analyses()

        # Configure QC of all panel data:
        self.configure_all_panel_qcs()

        # Configure fastq QCs:
        self.configure_fastq_qcs()

        # Configure the low-pass whole genome analysis:
        self.configure_lowpass_analyses()

        # Configure MultiQC:
        self.configure_multi_qc()

    def align_and_qdnaseq(self, lib):
        bam = align_library(self,
                            fq1_files=find_fastqs(lib, self.libdir)[0],
                            fq2_files=find_fastqs(lib, self.libdir)[1],
                            lib=lib,
                            ref=self.refdata['bwaIndex'],
                            outdir=self.outdir + "/bams/wgs",
                            maxcores=self.maxcores)

        qdnaseq = QDNASeq(bam,
                          output_segments="{}/cnv/{}-qdnaseq.segments.txt".format(self.outdir, lib),
                          background=None
                          )
        self.add(qdnaseq)

        return {'bam': bam}  # , 'qdnaseq-bed': qdnaseq.output_bed, 'qdnaseq-segments': qdnaseq.output_segments}

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
            isize.jobname = "picard-isize-{}".format(basefn)
            isize.output_metrics = "{}/qc/picard/wgs/{}.picard-insertsize.txt".format(self.outdir, basefn)
            self.add(isize)

            wgsmetrics = PicardCollectWgsMetrics()
            wgsmetrics.input = bam
            wgsmetrics.reference_sequence = self.refdata['reference_genome']
            wgsmetrics.output_metrics = "{}/qc/picard/wgs/{}.picard-wgsmetrics.txt".format(self.outdir, basefn)
            wgsmetrics.jobname = "picard-wgsmetrics-{}".format(basefn)
            self.add(wgsmetrics)

            qc_files += [isize.output_metrics, wgsmetrics.output_metrics]

        return qc_files