from pypedream.pipeline.pypedreampipeline import PypedreamPipeline
from autoseq.util.path import normpath, stripsuffix
from autoseq.tools.alignment import align_library
from autoseq.util.library import find_fastqs
from autoseq.tools.picard import PicardCollectInsertSizeMetrics, PicardCollectOxoGMetrics, \
    PicardMergeSamFiles, PicardMarkDuplicates, PicardCollectHsMetrics
from autoseq.tools.variantcalling import Freebayes, VEP, VcfAddSample, call_somatic_variants
from autoseq.tools.intervals import MsiSensor
from autoseq.tools.qc import *
import collections, logging


class CancerPanelResults(object):
    """"
    A collection of results produced by comparing a clinseq cancer sample library capture
    against a corresponding normal sample library capture.
    """

    def __init__(self, somatic_vcf, msi_output, hzconcordance_output, normal_contest_output,
                 cancer_contest_output, cancer_contam_call):
        self.somatic_vcf = somatic_vcf
        self.msi_output = msi_output
        self.hzconcordance_output = hzconcordance_output
        self.normal_contest_output = normal_contest_output
        self.cancer_contest_output = cancer_contest_output
        self.cancer_contam_call = cancer_contam_call


def compose_sample_str(sample_library_capture_tup):
    return "{}-{}-{}-{}".format(sample_library_capture_tup[0],
                                sample_library_capture_tup[1],
                                sample_library_capture_tup[2],
                                sample_library_capture_tup[3])


class ClinseqPipeline(PypedreamPipeline):
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

        # Overall dictionary of unique captures (tuples) to corresponding merged bams:
        self.capture_to_merged_bam = collections.defaultdict(list)

        # Dictionary of unique normal captures (tuples) to results obtained by
        # analysis based on that normal capture. Results values tuples containing:
        # - The resulting germlineVCF,
        # - A dictionary of with unique cancer captures (tuples) as keys and
        # CancerPanelResults objects as values.
        self.normal_capture_to_results = {}

    def get_unique_normal_captures(self):
        """
        Obtain tuples for all unique normal sample library captures in this pipeline instance.

        :return: List of tuples.
        """

        all_capture_tuples = self.capture_to_merged_bam.keys()
        return filter(lambda curr_tup: curr_tup[0] == "N", all_capture_tuples)

    def get_vep(self):
        vep = False
        if self.refdata['vep_dir']:
            vep = True

        return vep

    def get_unique_cancer_captures(self):
        """
        Obtain tuples for all unique cancer sample library captures in this pipeline instance.

        :return: List of tuples.
        """

        all_capture_tuples = self.capture_to_merged_bam.keys()
        return filter(lambda curr_tup: curr_tup[0] != "N", all_capture_tuples)

    def get_prep_kit_name(self, prep_kit_code):
        """
        Convert a two-letter library kit code to the corresponding library kit name.

        :param prep_kit_code: Two-letter library prep code. 
        :return: The library prep kit name.
        """

        # FIXME: Move this information to a config JSON file.
        prep_kit_lookup = {"BN": "BIOO_NEXTFLEX",
                           "KH": "KAPA_HYPERPREP",
                           "TD": "THRUPLEX_DNASEQ",
                           "TP": "THRUPLEX_PLASMASEQ",
                           "TF": "THRUPLEX_FD",
                           "TS": "TRUSEQ_RNA",
                           "NN": "NEBNEXT_RNA",
                           "VI": "VILO_RNA"}

        return prep_kit_lookup[prep_kit_code]

    def get_capture_name(self, capture_kit_code):
        """
        Convert a two-letter capture kit code to the corresponding capture kit name.

        :param capture_kit_code: The two-letter capture kit code.
        :return: The capture-kit name.
        """

        # FIXME: Move this information to a config JSON file.
        capture_kit_loopkup = {"CS": "clinseq_v3_targets",
                               "CZ": "clinseq_v4",
                               "EX": "EXOMEV3",
                               "EO": "EXOMEV1",
                               "RF": "fusion_v1",
                               "CC": "core_design",
                               "CD": "discovery_coho",
                               "CB": "big_design",
                               "AL": "alascca_targets",
                               "TT": "test-regions",
                               "CP": "progression",
                               "CM": "monitor"
                               }

        if capture_kit_code == 'WG':
            return 'lowpass_wgs'

        else:
            return capture_kit_loopkup[capture_kit_code]

    def get_all_clinseq_barcodes(self):
        """
        :return: All clinseq barcodes included in this clinseq analysis pipeline's panel data.
        """
        all_panel_clinseq_barcodes = \
            self.sampledata['panel']['T'] + \
            self.sampledata['panel']['N'] + \
            self.sampledata['panel']['CFDNA']
        return filter(lambda bc: bc != None, all_panel_clinseq_barcodes)

    def get_unique_capture_to_clinseq_barcodes(self):
        """
        Retrieves all clinseq barcodes for this clinseq analysis, and organises them according
        to unique library captures.

        :return: A dictionary with tuples indicating unique library captures as keys,
        and barcode lists as values. 
        """
        capture_to_barcodes = collections.defaultdict(list)
        for clinseq_barcode in self.get_all_clinseq_barcodes():
            capture_tuple = parse_capture_tuple(clinseq_barcode)
            capture_to_barcodes[capture_tuple].append(clinseq_barcode)

        return capture_to_barcodes

    def merge_and_rm_dup(self, sample_type, sample_id, prep_kit_id, capture_kit_id, input_bams):
        """
        Configures Picard merging and duplicate marking, for the specified group input bams,
        which should all correspond to the specified sample library capture.

        :param sample_type: Clinseq sample type
        :param sample_id: Sample ID
        :param prep_kit_id: Two-letter prep kit ID
        :param capture_kit_id: Two-letter capture kit ID
        :input_bams: The bam filenames for which to do merging and duplicate marking

        :return: The filename of the output merged bam
        """

        # Strings indicating the sample and capture, for use in output file names below:
        sample_str = "{}-{}".format(sample_type, sample_id)
        capture_str = "{}-{}-{}".format(sample_str, prep_kit_id, capture_kit_id)

        # Configure merging:
        merged_bam_filename = \
            "{}/bams/panel/{}.bam".format(self.outdir, capture_str)
        merge_bams = PicardMergeSamFiles(input_bams, merged_bam_filename)
        merge_bams.is_intermediate = True
        merge_bams.jobname = "picard-mergebams-{}".format(sample_str)
        self.add(merge_bams)

        # Configure duplicate marking:
        mark_dups_bam_filename = \
            "{}/bams/panel/{}-nodups.bam".format(self.outdir, capture_str)
        mark_dups_metrics_filename = \
            "{}/qc/picard/panel/{}-markdups-metrics.txt".format(self.outdir, capture_str)
        markdups = PicardMarkDuplicates(\
            merge_bams.output_bam, mark_dups_bam_filename, mark_dups_metrics_filename)
        markdups.is_intermediate = False
        self.add(markdups)

        self.qc_files.append(markdups.output_metrics)
        return markdups.output_bam

    def configure_align_and_merge(self):
        capture_to_barcodes = self.get_unique_capture_to_clinseq_barcodes()
        for capture_tuple in capture_to_barcodes.keys():
            curr_bamfiles = []
            for clinseq_barcode in capture_to_barcodes[capture_tuple]:
                curr_bamfiles.append(\
                    align_library(self,
                                  fq1_files=find_fastqs(clinseq_barcode, self.libdir)[0],
                                  fq2_files=find_fastqs(clinseq_barcode, self.libdir)[1],
                                  lib=clinseq_barcode,
                                  ref=self.refdata['bwaIndex'],
                                  outdir=self.outdir + "/bams/panel",
                                  maxcores=self.maxcores))

            self.capture_to_merged_bam = self.merge_and_rm_dup(\
                capture_tuple[0], capture_tuple[1], capture_tuple[2], capture_tuple[3], curr_bamfiles)

    def call_germline_variants(self, sample_id, prep_kit_id, capture_kit_id, bam):
        """
        Call germline variants for a normal sample library capture, and run VEP if configured.

        :param sample_id: Sample ID 
        :param prep_kit_id: Library prep two-letter code 
        :param capture_kit_id: Panel capture two-letter code
        :param bam: Corresponding library capture bam filename
        :return: Germline VCF filename
        """

        targets = self.get_capture_name(capture_kit_id)
        capture_str = "{}-{}-{}".format(sample_id, prep_kit_id, capture_kit_id)

        freebayes = Freebayes()
        freebayes.input_bams = [bam]
        freebayes.somatic_only = False
        freebayes.params = None
        freebayes.reference_sequence = self.refdata['reference_genome']
        freebayes.target_bed = self.refdata['targets'][targets]['targets-bed-slopped20']
        freebayes.threads = self.maxcores
        freebayes.scratch = self.scratch
        freebayes.output = "{}/variants/{}.freebayes-germline.vcf.gz".format(self.outdir, capture_str)
        freebayes.jobname = "freebayes-germline-{}".format(capture_str)
        self.add(freebayes)

        if self.refdata['vep_dir']:
            vep_freebayes = VEP()
            vep_freebayes.input_vcf = freebayes.output
            vep_freebayes.threads = self.maxcores
            vep_freebayes.reference_sequence = self.refdata['reference_genome']
            vep_freebayes.vep_dir = self.refdata['vep_dir']
            vep_freebayes.output_vcf = "{}/variants/{}.freebayes-germline.vep.vcf.gz".format(self.outdir, capture_str)
            vep_freebayes.jobname = "vep-freebayes-germline-{}".format(capture_str)
            self.add(vep_freebayes)

            return vep_freebayes.output_vcf
        else:
            return freebayes.output

    def configure_panel_analysis_with_normal(self, normal_capture_tuple):
        """
        Configure panel analyses focused on a specific unique normal library capture.
        """

        if normal_capture_tuple[0] is not "N":
            raise ValueError("Invalid input tuple: " + normal_capture_tuple)

        normal_sample_id = normal_capture_tuple[1]
        normal_library_prep_id = normal_capture_tuple[2]
        normal_capture_kit_id = normal_capture_tuple[3]

        normal_bam = self.capture_to_merged_bam[normal_capture_tuple]

        # Configure germline variant calling:
        germline_vcf = self.call_germline_variants(\
            normal_sample_id, normal_library_prep_id, normal_capture_kit_id, normal_bam)

        # For each unique cancer library capture, configure a comparative analysis against
        # this normal capture:
        cancer_capture_to_results = {}
        for cancer_capture in self.get_unique_cancer_captures():
            cancer_capture_to_results[cancer_capture] = \
                self.configure_panel_analysis_cancer_vs_normal(\
                    normal_capture_tuple, cancer_capture, germline_vcf)

        self.normal_capture_to_results = \
            (germline_vcf, cancer_capture_to_results)

    def configure_panel_analyses(self):
        """
        Configure generic analyses of all panel data for this clinseq pipeline.

        Populates self.capture_to_merged_bam and normal_capture_to_results. 
        """

        # Configure alignment and merging for each unique sample library capture:
        self.configure_align_and_merge()

        # Configure a separate group of analyses for each unique normal library capture:
        for normal_capture in self.get_unique_normal_captures():
            self.configure_panel_analysis_with_normal(normal_capture)

    def configure_panel_analysis_cancer_vs_normal(self, normal_capture_tuple, cancer_capture_tuple, normal_vcf):
        normal_bam = self.capture_to_merged_bam[normal_capture_tuple]
        normal_sample_str = compose_sample_str(normal_capture_tuple)

        cancer_bam = self.capture_to_merged_bam[cancer_capture_tuple]
        cancer_sample_str = compose_sample_str(cancer_capture_tuple)
        cancer_target_name = self.get_capture_name(cancer_capture_tuple[3])

        # Configure somatic variant calling:
        somatic_variants = call_somatic_variants(self, tbam=cancer_bam, nbam=normal_bam, tlib=cancer_sample_str,
                                                 nlib=normal_sample_str, target_name=cancer_target_name,
                                                 refdata=self.refdata, outdir=self.outdir,
                                                 callers=['vardict'], vep=self.get_vep())

        # Configure VCF add sample:
        vcfaddsample = VcfAddSample()
        vcfaddsample.input_bam = cancer_bam
        vcfaddsample.input_vcf = normal_vcf
        vcfaddsample.samplename = cancer_sample_str
        vcfaddsample.filter_hom = True
        vcfaddsample.output = "{}/variants/{}-and-{}.germline-variants-with-somatic-afs.vcf.gz".format( \
            self.outdir, normal_sample_str, cancer_sample_str)
        vcfaddsample.jobname = "vcf-add-sample-{}".format(cancer_sample_str)
        self.add(vcfaddsample)

        # Configure MSI sensor:
        msisensor = MsiSensor()
        msisensor.msi_sites = self.refdata['targets'][cancer_target_name]['msisites']
        msisensor.input_normal_bam = normal_bam
        msisensor.input_tumor_bam = cancer_bam
        msisensor.output = "{}/msisensor.tsv".format(self.outdir)
        msisensor.threads = self.maxcores
        msisensor.jobname = "msisensor-{}".format(cancer_sample_str)
        self.add(msisensor)

        hzconcordance = HeterzygoteConcordance()
        hzconcordance.input_vcf = normal_vcf
        hzconcordance.input_bam = cancer_bam
        hzconcordance.reference_sequence = self.refdata['reference_genome']
        hzconcordance.target_regions = self.refdata['targets'][cancer_target_name]['targets-interval_list-slopped20']
        hzconcordance.normalid = normal_sample_str
        hzconcordance.filter_reads_with_N_cigar = True
        hzconcordance.jobname = "hzconcordance-{}".format(cancer_sample_str)
        hzconcordance.output = "{}/bams/{}-{}-hzconcordance.txt".format(self.outdir, cancer_sample_str,
                                                                        normal_sample_str)
        self.add(hzconcordance)

        # FIXME: ADD IN CONTAMINATION CALL CALCULATION USING CONTEST AND SUBSEQUENT
        normal_contest_output = None
        cancer_contest_output = None
        cancer_contam_call = None

        return CancerPanelResults(somatic_variants, msisensor.output, hzconcordance.output,
                                  normal_contest_output, cancer_contest_output, cancer_contam_call)

    def configure_all_panel_qcs(self):
        for capture_tup in self.capture_to_merged_bam.keys():
            self.qc_files += \
                self.configure_panel_bam_qc(self.capture_to_merged_bam[capture_tup],
                                      capture_tup[3])

    def configure_multi_qc(self):
        multiqc = MultiQC()
        multiqc.input_files = self.qc_files
        multiqc.search_dir = self.outdir
        multiqc.output = "{}/multiqc/{}-multiqc".format(self.outdir, self.sampledata['sdid'])
        multiqc.jobname = "multiqc-{}".format(self.sampledata['sdid'])
        self.add(multiqc)

    def configure_panel_bam_qc(self, bam, capture_kit_code):
        """
        Run QC on panel bams
        :param bams: list of bams
        :return: list of generated files
        """

        targets = self.get_capture_name(capture_kit_code)
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

        return [isize.output_metrics, oxog.output_metrics, hsmetrics.output_metrics, sambamba.output]


def parse_capture_tuple(clinseq_barcode):
    """
    Convenience function for use in the context of joint panel analysis.

    Extracts the sample type, sample ID, library prep ID, and capture kit ID,
    from the specified clinseq barcode.

    :param clinseq_barcode: List of one or more clinseq barcodes 
    :return: (sample type, sample ID, capture kit ID) tuple
    """
    return (parse_sample_type(clinseq_barcode),
            parse_sample_id(clinseq_barcode),
            parse_prep_kit_id(clinseq_barcode),
            parse_capture_kit_id(clinseq_barcode))


def parse_sample_type(clinseq_barcode):
    return clinseq_barcode.split("-")[3]


def parse_sample_id(clinseq_barcode):
    return clinseq_barcode.split("-")[4]


def parse_prep_kit_id(clinseq_barcode):
    return clinseq_barcode.split("-")[5][:2]


def parse_capture_kit_id(clinseq_barcode):
    return clinseq_barcode.split("-")[6][:2]