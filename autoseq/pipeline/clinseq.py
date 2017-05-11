from pypedream.pipeline.pypedreampipeline import PypedreamPipeline
from autoseq.util.path import normpath
from autoseq.tools.alignment import align_library
from autoseq.util.library import find_fastqs
from autoseq.tools.picard import PicardMergeSamFiles, PicardMarkDuplicates
import collections


class CancerPanelResults(object):
    """"
    A collection of results produced by comparing a clinseq cancer sample library capture
    against a corresponding normal sample library capture.
    """

    def __init__(self):
        self.somatic_vcf = None
        self.hzconcordance_output = None
        self.msi_output = None
        self.contam_call = None


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

        self.panel_germline_vcf = None
        self.capture_to_merged_bam = collections.defaultdict(list)
        self.cancer_capture_to_results = collections.defaultdict(list)

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


    def configure_panel_analysis(self):
        """
        Configure generic analysis of all panel data for this clinseq pipeline.
        """

        # Configure alignment and merging for each unique sample library capture:
        self.configure_align_and_merge()

        # Configure germline calling:
        self.germline_vcf = self.call_germline_variants(
            clinseq_barcode_to_bamfile[normal_clinseq_barcode], normal_clinseq_barcode)
        




def analyze_panel(pipeline):
    """
Configures core analysis of all panel-captured libraries for a given pipeline.

The panel-capture libraries can include:
- Zero or one normal library captures
- Zero or one tumor library captures
- Zero or more cfDNA library captures

Configures the following core analyses:
- Alignment for all normal, cfDNA, and tumor library captures
- Merging and duplicate removal for each unique tumor and cfDNA (sample, capture kit) pairing.
- If germline is present:
-- Germline calling
-- Additional core non_normal vs normal analyses

    :param pipeline: Analysis pipeline specifying the library captures to analyse and other relevant parameters
    :return: A dictionary with (sample type, sample barcode, capture kit code) tuples as keys, and CancerPanelResults
    objects as values.
    """

    # If there is a normal library capture, then do additional core panel analyses
    # supported by that item:
    if normal_clinseq_barcode != None:
        # Configure germline calling:
        

        # For each unique cfDNA and tumor capture:
        for capture_tup in capture_tup_to_merged_bam.keys():
            # Configure additional core analysis with the resulting bam file together
            # with the normal bam file and the germline VCF:
            analyze_panel_cancer_vs_normal(pipeline, capture_tup[0], capture_tup[1], capture_tup[2],
                                           capture_tup_to_merged_bam[capture_tup],
                                           clinseq_barcode_to_bamfile[normal_clinseq_barcode],
                                           germline_vcf,
                                           normal_clinseq_barcode)

    # Create and return a merged cancer and normal final bam file dictionary:
    # XXX CONTINUE HERE; FIGURE OUT THE REQUIRED STRUCTURE FOR THIS FINAL OUTPUT DICTIONARY, AND IMPLEMENT IT
    # HERE BY AGGREGATING THE REQUIRED DATA. Can work with + adapt these data structures and functions:
    # XXX IS IT / WILL IT ALWAYS BE SUFFICIENT TO HAVE ONLY THE SAME TYPE, SIMPLE ID, AND CAPTURE KIT ID FOR
    # EACH OF THE CANCER MERGED BAM FILES? If so then just include the capture_tup tuples as keys in the dictionary
    # as planned.
    # get_capture_tup_to_clinseq_barcodes(normal_clinseq_barcode)
    # capture_tup_to_merged_bam[]


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
