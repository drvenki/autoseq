import collections
import json
import logging

from pypedream.pipeline.pypedreampipeline import PypedreamPipeline

from autoseq.tools.alignment import align_library
from autoseq.tools.cnvcalling import QDNASeq
from autoseq.tools.intervals import MsiSensor
from autoseq.tools.picard import PicardCollectGcBiasMetrics, PicardCollectWgsMetrics, \
    PicardMergeSamFiles, PicardMarkDuplicates, PicardCollectHsMetrics
from autoseq.tools.picard import PicardCollectInsertSizeMetrics
from autoseq.tools.picard import PicardCollectOxoGMetrics
from autoseq.tools.qc import *
from autoseq.tools.unix import Copy
from autoseq.tools.variantcalling import Mutect2, Freebayes, VEP, VcfAddSample, VarDict, call_somatic_variants
from autoseq.util.library import get_libdict, get_capture_kit_name_from_id, find_fastqs
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

    def analyze_panel(self):
        tbam = None
        nbam = None
        pbams = []
        germline_vcf = None
        somatic_variants = {}
        tlib = self.sampledata['panel']['T']
        nlib = self.sampledata['panel']['N']
        plibs = self.sampledata['panel']['CFDNA']

        # align germline normal
        if nlib:
            nbam = align_library(self,
                                 fq1_files=find_fastqs(nlib, self.libdir)[0],
                                 fq2_files=find_fastqs(nlib, self.libdir)[1],
                                 lib=nlib,
                                 ref=self.refdata['bwaIndex'],
                                 outdir=self.outdir + "/bams/panel",
                                 maxcores=self.maxcores)

            germline_vcf = self.call_germline_variants(nbam, library=nlib)

        # process tumor and plasma samples
        libs = [x for x in plibs + [tlib] if x is not None]
        sample_bams = collections.defaultdict(list)
        sample_dicts = collections.defaultdict(list)
        for lib in libs:
            libdict = get_libdict(lib)
            sample = "{}-{}-{}".format(libdict['sdid'], libdict['type'], libdict['sample_id'])
            bam = align_library(self,
                                fq1_files=find_fastqs(lib, self.libdir)[0],
                                fq2_files=find_fastqs(lib, self.libdir)[1],
                                lib=lib,
                                ref=self.refdata['bwaIndex'],
                                outdir=self.outdir + "/bams/panel",
                                maxcores=self.maxcores,
                                remove_duplicates=False)
            sample_bams[sample].append(bam)
            sample_dicts[sample].append(libdict)

        for sample, bams in sample_bams.items():
            # don't allow samples with different captures to be merged.
            capture_kits_used_for_sample = [sampledict['capture_id'][0:2] for sampledict in sample_dicts[sample]]
            prep_kits_used_for_sample = [sampledict['prep_id'][0:2] for sampledict in sample_dicts[sample]]
            if len(set(capture_kits_used_for_sample)) > 1:
                raise ValueError("Multiple capture kits used for libraries for sample {} ({})".format(
                    sample, sample_dicts[sample]))

            targets_short = capture_kits_used_for_sample[0]  # short name, ex CB
            targets_long = get_capture_kit_name_from_id(capture_kits_used_for_sample[0])  # long name, ex big_design
            prep_kit_short = prep_kits_used_for_sample[0]

            merge_bams = PicardMergeSamFiles(input_bams=bams,
                                             output_bam="{}/bams/panel/{}-{}-{}.bam".format(
                                                 self.outdir, sample, prep_kit_short, targets_short))

            merge_bams.is_intermediate = True
            merge_bams.jobname = "picard-mergebams-{}".format(sample)
            self.add(merge_bams)

            markdups = PicardMarkDuplicates(merge_bams.output_bam,
                                            output_bam="{}/bams/panel/{}-{}-{}-nodups.bam".format(
                                                self.outdir, sample, prep_kit_short, targets_short),
                                            output_metrics="{}/qc/picard/panel/{}-{}-{}-markdups-metrics.txt".format(
                                                self.outdir, sample, prep_kit_short, targets_short))
            markdups.is_intermediate = False
            self.qc_files.append(markdups.output_metrics)
            self.add(markdups)

            if set([smp['type'] for smp in sample_dicts[sample]]) == set('CFDNA'):
                # if it's a plasma sample
                pbams.append(markdups.output_bam)
            elif set([smp['type'] for smp in sample_dicts[sample]]) == set('T'):
                # if it's a tumor sample
                tbam = markdups.output_bam

            vep = False
            if self.refdata['vep_dir']:
                vep = True

            if nlib:
                somatic_variants = call_somatic_variants(self, tbam=markdups.output_bam, nbam=nbam, tlib=sample,
                                                         nlib=nlib, target_name=targets_long, refdata=self.refdata,
                                                         outdir=self.outdir,
                                                         callers=['vardict'],
                                                         vep=vep)

                vcfaddsample = VcfAddSample()
                vcfaddsample.input_bam = markdups.output_bam
                vcfaddsample.input_vcf = germline_vcf
                vcfaddsample.samplename = sample
                vcfaddsample.filter_hom = True
                vcfaddsample.output = "{}/variants/{}-and-{}.germline-variants-with-somatic-afs.vcf.gz".format(
                    self.outdir,
                    nlib,
                    sample)
                vcfaddsample.jobname = "vcf-add-sample-{}".format(sample)
                self.add(vcfaddsample)

                msisensor = MsiSensor()
                msisensor.msi_sites = self.refdata['targets'][targets_long]['msisites']
                msisensor.input_normal_bam = nbam
                msisensor.input_tumor_bam = markdups.output_bam
                msisensor.output = "{}/msisensor.tsv".format(self.outdir)
                msisensor.threads = self.maxcores
                msisensor.jobname = "msisensor-{}".format(sample)
                self.add(msisensor)

                libdict = get_libdict(nlib)
                rg_sm = "{}-{}-{}".format(libdict['sdid'], libdict['type'], libdict['sample_id'])

                hzconcordance = HeterzygoteConcordance()
                hzconcordance.input_vcf = germline_vcf
                hzconcordance.input_bam = markdups.output_bam
                hzconcordance.reference_sequence = self.refdata['reference_genome']
                hzconcordance.target_regions = self.refdata['targets'][targets_long]['targets-interval_list-slopped20']
                hzconcordance.normalid = rg_sm
                hzconcordance.filter_reads_with_N_cigar = True
                hzconcordance.jobname = "hzconcordance-{}".format(lib)
                hzconcordance.output = "{}/bams/{}-{}-hzconcordance.txt".format(self.outdir, lib, nlib)
                self.add(hzconcordance)

        return {'tbam': tbam, 'nbam': nbam, 'pbams': pbams,
                'somatic_variants': somatic_variants}

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
