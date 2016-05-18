import json
import logging

from pypedream.pipeline.pypedreampipeline import PypedreamPipeline
from pypedream.tools.unix import Cat

from autoseq.tools.alignment import Bwa, SkewerSE, SkewerPE
from autoseq.tools.cnvcalling import QDNASeq
from autoseq.tools.intervals import MsiSensor
from autoseq.tools.picard import PicardCollectGcBiasMetrics, PicardCalculateHsMetrics, PicardCollectWgsMetrics
from autoseq.tools.picard import PicardCollectInsertSizeMetrics
from autoseq.tools.picard import PicardCollectOxoGMetrics
from autoseq.tools.qc import *
from autoseq.tools.variantcalling import Mutect2, Freebayes, VEP, VcfAddSample, VarDict
from autoseq.util.library import get_libdict
from autoseq.util.path import normpath, stripsuffix

__author__ = 'dankle'


class LiqBioPipeline(PypedreamPipeline):
    analysis_id = None
    sampledata = None
    refdata = None
    outdir = None
    maxcores = None
    scratch = None

    def __init__(self, sampledata, refdata, outdir, libdir, analysis_id=None, maxcores=1, scratch="/tmp", debug=False,
                 **kwargs):
        PypedreamPipeline.__init__(self, normpath(outdir), **kwargs)
        self.sampledata = sampledata
        self.refdata = refdata
        self.maxcores = maxcores
        self.analysis_id = analysis_id
        self.libdir = libdir
        self.scratch = scratch

        self.check_sampledata()

        panel_files = self.analyze_panel(debug=debug)

        wgs_bams = self.analyze_lowpass_wgs()

        ################################################
        # QC
        # qc_files = []
        #
        # # per-bam qc
        # # panel
        # all_panel_bams = [panel_files['tbam'], panel_files['nbam']] + panel_files['pbams']
        # all_panel_bams = [bam for bam in all_panel_bams if bam is not None]
        # logging.debug("Bam files are {}".format(all_panel_bams))
        # #qc_files += self.run_panel_bam_qc(all_panel_bams, debug=debug)
        # # wgs
        # all_wgs_bams = [bam for bam in wgs_bams.values() if bam is not None]
        # #qc_files += self.run_wgs_bam_qc(all_wgs_bams, debug=debug)
        #
        # # per-fastq qc
        # fqs = self.get_all_fastqs()
        # logging.debug("fqs = {}".format(fqs))
        # qc_files += self.run_fastq_qc(fqs)
        #
        # multiqc = MultiQC()
        # multiqc.input_files = qc_files
        # multiqc.search_dir = self.outdir
        # multiqc.output = "{}/multiqc/{}-multiqc".format(self.outdir, self.sampledata['sdid'])
        # multiqc.jobname = "multiqc-{}".format(self.sampledata['sdid'])
        # self.add(multiqc)

    def check_sampledata(self):
        def check_lib(lib):
            if lib:
                dir = os.path.join(self.libdir, lib)
                if not os.path.exists(dir):
                    logging.warn("Dir {} does not exists for {}. Not using library.".format(dir, lib))
                    return None
                if self.find_fastqs(lib) == (None, None):
                    logging.warn("No fastq files found for {} in dir {}".format(lib, dir))
                    return None
            logging.debug("Library {} has data. Using it.".format(lib))
            return lib

        for datatype in ['panel', 'wgs']:
            self.sampledata[datatype]['T'] = check_lib(self.sampledata[datatype]['T'])
            self.sampledata[datatype]['N'] = check_lib(self.sampledata[datatype]['N'])
            plibs_with_data = []
            for plib in self.sampledata[datatype]['P']:
                plib_checked = check_lib(plib)
                if plib_checked:
                    plibs_with_data.append(plib_checked)

            self.sampledata[datatype]['P'] = plibs_with_data

    def get_all_fastqs(self):
        fqs = []
        if self.sampledata['panel']['T']:
            fqs.extend(self.find_fastqs(self.sampledata['panel']['T'])[0])
            fqs.extend(self.find_fastqs(self.sampledata['panel']['T'])[1])
        if self.sampledata['panel']['N']:
            fqs.extend(self.find_fastqs(self.sampledata['panel']['N'])[0])
            fqs.extend(self.find_fastqs(self.sampledata['panel']['N'])[1])
        for plib in self.sampledata['panel']['P']:
            fqs.extend(self.find_fastqs(plib)[0])
            fqs.extend(self.find_fastqs(plib)[1])

        return [fq for fq in fqs if fq is not None]

    def find_fastqs(self, lib):
        """Find fastq files for a given library id, return a tuple of lists for _1 and _2 files."""
        d = os.path.join(self.libdir, lib)
        logging.debug("Looking for fastq files for library {} in {}".format(
            lib, d
        ))
        fq1s = ["{}/{}".format(d, f) for f in os.listdir(d) if
                f.endswith("_1.fastq.gz") or f.endswith("_1.fq.gz")]

        fq2s = ["{}/{}".format(d, f) for f in os.listdir(d) if
                f.endswith("_2.fastq.gz") or f.endswith("_2.fq.gz")]

        logging.debug("Found {}".format((fq1s, fq2s)))
        return fq1s, fq2s

    def analyze_lowpass_wgs(self):
        tbam = None
        nbam = None
        pbams = []

        tlib = self.sampledata['wgs']['T']
        nlib = self.sampledata['wgs']['N']
        plibs = self.sampledata['wgs']['P']

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
        bam = self.align_library(fq1_files=self.find_fastqs(lib)[0],
                                 fq2_files=self.find_fastqs(lib)[1],
                                 lib=lib,
                                 ref=self.refdata['bwaIndex'],
                                 outdir=self.outdir + "/bams/wgs",
                                 maxcores=self.maxcores)

        qdnaseq = QDNASeq()
        qdnaseq.input = bam
        qdnaseq.output_bed = "{}/cnv/{}-qdnaseq.bed".format(self.outdir, lib)
        qdnaseq.output_segments = "{}/cnv/{}-qdnaseq.segments.txt".format(self.outdir, lib)
        qdnaseq.genes_gtf = self.refdata['genesGtfGenesOnly']
        qdnaseq.background = self.refdata["qdnaseq_background"]
        self.add(qdnaseq)

        return {'bam': bam, 'qdnaseq-bed': qdnaseq.output_bed, 'qdnaseq-segments': qdnaseq.output_segments}

    def analyze_panel(self, debug=False):
        tbam = None
        nbam = None
        pbams = []
        somatic_vcfs = []
        germline_vcf = None
        tlib = self.sampledata['panel']['T']
        nlib = self.sampledata['panel']['N']
        plibs = self.sampledata['panel']['P']

        # align germline normal
        if nlib:
            nbam = self.align_library(fq1_files=self.find_fastqs(nlib)[0],
                                      fq2_files=self.find_fastqs(nlib)[1],
                                      lib=nlib,
                                      ref=self.refdata['bwaIndex'],
                                      outdir=self.outdir + "/bams/panel",
                                      maxcores=self.maxcores)

            germline_vcf = self.call_germline_variants(nbam, library=nlib)

        # process tumor and plasma samples
        libs = [x for x in plibs + [tlib] if x is not None]
        for lib in libs:
            bam = self.align_library(fq1_files=self.find_fastqs(lib)[0],
                                     fq2_files=self.find_fastqs(lib)[1],
                                     lib=lib,
                                     ref=self.refdata['bwaIndex'],
                                     outdir=self.outdir + "/bams/panel",
                                     maxcores=self.maxcores)
            if nlib:
                #  If we have a normal, call variants, verify identity and run msisensor

                somatic_vcfs.append(self.call_somatic_variants(bam, nbam))

                targets = get_libdict(lib)['capture_kit_name']
                hzconcordance = HeterzygoteConcordance()
                hzconcordance.input_vcf = germline_vcf
                hzconcordance.input_bam = bam
                hzconcordance.reference_sequence = self.refdata['reference_genome']
                hzconcordance.target_regions = self.refdata['targets'][targets]['targets-interval_list-slopped20']
                hzconcordance.normalid = nlib
                hzconcordance.filter_reads_with_N_cigar = True
                hzconcordance.jobname = "hzconcordance-{}".format(lib)
                hzconcordance.output = "{}/bams/{}-{}-hzconcordance.txt".format(self.outdir, lib, nlib)
                self.add(hzconcordance)

                vcfaddsample = VcfAddSample()
                vcfaddsample.input_bam = bam
                vcfaddsample.input_vcf = germline_vcf
                vcfaddsample.samplename = lib
                vcfaddsample.filter_hom = True
                vcfaddsample.output = "{}/variants/{}-and-{}.germline.vcf.gz".format(self.outdir,
                                                                                     lib,
                                                                                     nlib)
                vcfaddsample.jobname = "vcf-add-sample-{}".format(lib)
                self.add(vcfaddsample)

                msisensor = MsiSensor()
                msisensor.msi_sites = self.refdata['targets'][targets]['msisites']
                msisensor.input_normal_bam = nbam
                msisensor.input_tumor_bam = bam
                msisensor.output = "{}/msisensor.tsv".format(self.outdir)
                msisensor.threads = self.maxcores
                msisensor.jobname = "msisensor-{}".format(lib)
                if not debug:
                    self.add(msisensor)

        return {'tbam': tbam, 'nbam': nbam, 'pbams': pbams,
                'somatic_vcfs': somatic_vcfs}

    def call_germline_variants(self, bam, library):
        """
        Call germline variants from a bam and run VEP on it
        :param bam:
        :return:
        """
        targets = get_libdict(library)['capture_kit_name']
        freebayes = Freebayes()
        freebayes.input_bams = [bam]
        freebayes.normalid = library
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

    def call_somatic_variants(self, tbam, nbam):
        """
        Call somatic variants with freebayes and vardict.
        :param tbam: tumor bam
        :param nbam: normal bam
        :return:
        """
        tlib = stripsuffix(os.path.basename(tbam), ".bam")
        nlib = stripsuffix(os.path.basename(nbam), ".bam")

        targets = get_libdict(tlib)['capture_kit_name']

        mutect2 = Mutect2()
        mutect2.input_tumor = tbam
        mutect2.input_normal = nbam
        mutect2.reference_sequence = self.refdata['reference_genome']
        mutect2.target_regions = self.refdata['targets'][targets]['targets-interval_list-slopped20']
        mutect2.scratch = self.scratch
        mutect2.jobname = "mutect2-{}".format(tlib)
        mutect2.output = "{outdir}/variants/{tlib}/{tlib}-{nlib}.mutect2.vcf.gz".format(outdir=self.outdir,
                                                                                        tlib=tlib,
                                                                                        nlib=nlib)
        self.add(mutect2)

        freebayes = Freebayes()
        freebayes.input_bams = [tbam, nbam]
        freebayes.tumorid = tlib
        freebayes.normalid = nlib
        freebayes.somatic_only = True
        freebayes.reference_sequence = self.refdata['reference_genome']
        freebayes.target_bed = self.refdata['targets'][targets]['targets-bed-slopped20']
        freebayes.threads = self.maxcores
        freebayes.scratch = self.scratch
        freebayes.jobname = "freebayes-somatic-{}".format(tlib)
        freebayes.output = "{outdir}/variants/{tlib}/{tlib}-{nlib}.freebayes-somatic.vcf.gz".format(outdir=self.outdir,
                                                                                                    tlib=tlib,
                                                                                                    nlib=nlib)

        self.add(freebayes)

        vardict = VarDict(input_tumor=tbam, input_normal=nbam, tumorid=tlib, normalid=nlib,
                          reference_sequence=self.refdata['reference_genome'],
                          target_bed=self.refdata['targets'][targets]['targets-bed-slopped20'],
                          output="{outdir}/variants/{tlib}/{tlib}-{nlib}.vardict-somatic.vcf.gz".format(
                              outdir=self.outdir,
                              tlib=tlib,
                              nlib=nlib)
                          )
        vardict.jobname = "vardict-{}".format(tlib)
        self.add(vardict)

        if self.refdata['vep_dir']:
            vep_mutect2 = VEP()
            vep_mutect2.input_vcf = vardict.output
            vep_mutect2.threads = self.maxcores
            vep_mutect2.reference_sequence = self.refdata['reference_genome']
            vep_mutect2.vep_dir = self.refdata['vep_dir']
            vep_mutect2.output_vcf = "{outdir}/variants/{tlib}/{tlib}-{nlib}.mutect2.vep.vcf.gz".format(
                outdir=self.outdir,
                tlib=tlib,
                nlib=nlib)
            vep_mutect2.jobname = "vep-mutect2-{}".format(tlib)
            self.add(vep_mutect2)

            vep_vardict = VEP()
            vep_vardict.input_vcf = vardict.output
            vep_vardict.threads = self.maxcores
            vep_vardict.reference_sequence = self.refdata['reference_genome']
            vep_vardict.vep_dir = self.refdata['vep_dir']
            vep_vardict.output_vcf = "{outdir}/variants/{tlib}/{tlib}-{nlib}.vardict-somatic.vep.vcf.gz".format(
                outdir=self.outdir,
                tlib=tlib,
                nlib=nlib)
            vep_vardict.jobname = "vep-vardict-{}".format(tlib)
            self.add(vep_vardict)

            vep_freebayes = VEP()
            vep_freebayes.input_vcf = freebayes.output
            vep_freebayes.threads = self.maxcores
            vep_freebayes.reference_sequence = self.refdata['reference_genome']
            vep_freebayes.vep_dir = self.refdata['vep_dir']
            vep_freebayes.output_vcf = "{outdir}/variants/{tlib}/{tlib}-{nlib}.freebayes-somatic.vep.vcf.gz".format(
                outdir=self.outdir,
                tlib=tlib,
                nlib=nlib)
            vep_freebayes.jobname = "vep-freebayes-somatic-{}".format(tlib)
            self.add(vep_freebayes)

            return [vep_freebayes.output_vcf, vep_vardict.output_vcf, vep_mutect2.output_vcf]
        else:
            return [freebayes.output, vardict.output, mutect2.output]

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
            lib = stripsuffix(os.path.basename(bam), ".bam")
            targets = get_libdict(lib)['capture_kit_name']
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
            skewer.scratch = self.scratch
            skewer.jobname = "skewer-{}".format(os.path.basename(fq1))
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
            skewer.scratch = self.scratch
            skewer.jobname = "skewer-{}".format(os.path.basename(fq1))
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
