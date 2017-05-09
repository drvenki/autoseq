import json
import logging

from pypedream.pipeline.pypedreampipeline import PypedreamPipeline

from autoseq.tools.alignment import align_library
from autoseq.tools.cnvcalling import QDNASeq
from autoseq.tools.picard import PicardCollectGcBiasMetrics, PicardCollectWgsMetrics, \
    PicardMergeSamFiles, PicardMarkDuplicates, PicardCollectHsMetrics
from autoseq.tools.picard import PicardCollectInsertSizeMetrics
from autoseq.tools.picard import PicardCollectOxoGMetrics
from autoseq.tools.qc import *
from autoseq.tools.unix import Copy
from autoseq.tools.variantcalling import Freebayes, VEP
from autoseq.util.library import get_libdict, find_fastqs
from autoseq.util.path import normpath, stripsuffix

__author__ = 'dankle'


class LiqBioPipeline(PypedreamPipeline):
    def __init__(self, sampledata, refdata, outdir, libdir, analysis_id=None, maxcores=1, scratch="/scratch/tmp/tmp",
                 **kwargs):
        PypedreamPipeline.__init__(self, normpath(outdir), **kwargs)
        self.sampledata = sampledata
        self.refdata = refdata
        self.maxcores = maxcores
        self.analysis_id = analysis_id
        self.libdir = libdir
        self.qc_files = []
        self.scratch = scratch

        self.check_sampledata()

        panel_files = self.analyze_panel()

        wgs_bams = self.analyze_lowpass_wgs()

        ################################################
        # QC

        #
        # # per-bam qc
        # # panel
        all_panel_bams = [panel_files['tbam'], panel_files['nbam']] + panel_files['pbams']
        all_panel_bams = [bam for bam in all_panel_bams if bam is not None]
        self.qc_files += self.run_panel_bam_qc(all_panel_bams)

        # # wgs
        # all_wgs_bams = [bam for bam in wgs_bams.values() if bam is not None]
        # #qc_files += self.run_wgs_bam_qc(all_wgs_bams)
        #
        # # per-fastq qc
        # fqs = self.get_all_fastqs()
        # logging.debug("fqs = {}".format(fqs))
        # qc_files += self.run_fastq_qc(fqs)
        #

        multiqc = MultiQC()
        multiqc.input_files = self.qc_files
        multiqc.search_dir = self.outdir
        multiqc.output = "{}/multiqc/{}-multiqc".format(self.outdir, self.sampledata['sdid'])
        multiqc.jobname = "multiqc-{}".format(self.sampledata['sdid'])
        self.add(multiqc)

    def check_sampledata(self):
        def check_lib(lib):
            if lib:
                filedir = os.path.join(self.libdir, lib)
                if not os.path.exists(filedir):
                    logging.warn("Dir {} does not exists for {}. Not using library.".format(filedir, lib))
                    return None
                if find_fastqs(lib, self.libdir) == (None, None):
                    logging.warn("No fastq files found for {} in dir {}".format(lib, filedir))
                    return None
            logging.debug("Library {} has data. Using it.".format(lib))
            return lib

        for datatype in ['panel', 'wgs']:
            self.sampledata[datatype]['T'] = check_lib(self.sampledata[datatype]['T'])
            self.sampledata[datatype]['N'] = check_lib(self.sampledata[datatype]['N'])
            plibs_with_data = []
            for plib in self.sampledata[datatype]['CFDNA']:
                plib_checked = check_lib(plib)
                if plib_checked:
                    plibs_with_data.append(plib_checked)

            self.sampledata[datatype]['CFDNA'] = plibs_with_data

    def get_all_fastqs(self):
        fqs = []
        if self.sampledata['panel']['T']:
            fqs.extend(find_fastqs(self.sampledata['panel']['T'], self.libdir)[0])
            fqs.extend(find_fastqs(self.sampledata['panel']['T'], self.libdir)[1])
        if self.sampledata['panel']['N']:
            fqs.extend(find_fastqs(self.sampledata['panel']['N'], self.libdir)[0])
            fqs.extend(find_fastqs(self.sampledata['panel']['N'], self.libdir)[1])
        for plib in self.sampledata['panel']['CFDNA']:
            fqs.extend(find_fastqs(plib, self.libdir)[0])
            fqs.extend(find_fastqs(plib, self.libdir)[1])

        return [fq for fq in fqs if fq is not None]

    def analyze_lowpass_wgs(self):
        tbam = None
        nbam = None
        pbams = []

        tlib = self.sampledata['wgs']['T']
        nlib = self.sampledata['wgs']['N']
        plibs = self.sampledata['wgs']['CFDNA']

        if nlib:
            nfiles = self.align_and_qdnaseq(nlib)
            nbam = nfiles['bam']

        if tlib:
            tfiles = self.align_and_qdnaseq(tlib)
            tbam = tfiles['bam']

        for plib in plibs:
            pfiles = self.align_and_qdnaseq(plib)
            pbam = pfiles['bam']
            pbams.append(pbam)

        return {'tbam': tbam, 'nbam': nbam, 'pbams': pbams}

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

    def call_germline_variants(self, bam, library):
        """
        Call germline variants from a bam and run VEP on it
        :param bam:
        :return:
        """
        targets = get_libdict(library)['capture_kit_name']
        freebayes = Freebayes()
        freebayes.input_bams = [bam]
        freebayes.somatic_only = False
        freebayes.params = None
        freebayes.reference_sequence = self.refdata['reference_genome']
        freebayes.target_bed = self.refdata['targets'][targets]['targets-bed-slopped20']
        freebayes.threads = self.maxcores
        freebayes.scratch = self.scratch
        freebayes.output = "{}/variants/{}.freebayes-germline.vcf.gz".format(self.outdir, library)
        freebayes.jobname = "freebayes-germline-{}".format(library)
        self.add(freebayes)

        if self.refdata['vep_dir']:
            vep_freebayes = VEP()
            vep_freebayes.input_vcf = freebayes.output
            vep_freebayes.threads = self.maxcores
            vep_freebayes.reference_sequence = self.refdata['reference_genome']
            vep_freebayes.vep_dir = self.refdata['vep_dir']
            vep_freebayes.output_vcf = "{}/variants/{}.freebayes-germline.vep.vcf.gz".format(self.outdir, library)
            vep_freebayes.jobname = "vep-freebayes-germline-{}".format(library)
            self.add(vep_freebayes)

            return vep_freebayes.output_vcf
        else:
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
            fastqc.jobname = "fastqc-{}".format(basefn)
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

    def run_panel_bam_qc(self, bams):
        """
        Run QC on panel bams
        :param bams: list of bams
        :return: list of generated files
        """

        qc_files = []
        for bam in bams:
            lib = stripsuffix(os.path.basename(bam), ".bam")
            lib = stripsuffix(lib, '-nodups')
            targets = get_libdict(lib)['capture_kit_name']
            logging.debug("Adding QC jobs for {}".format(bam))
            basefn = stripsuffix(os.path.basename(bam), ".bam")
            isize = PicardCollectInsertSizeMetrics()
            isize.input = bam
            isize.output_metrics = "{}/qc/picard/panel/{}.picard-insertsize.txt".format(self.outdir, basefn)
            isize.jobname = "picard-isize-{}".format(basefn)
            self.add(isize)

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

            qc_files += [isize.output_metrics, oxog.output_metrics,
                         hsmetrics.output_metrics, sambamba.output]

        return qc_files
