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

    def configure_panel_analysis(self):
        """
        Configure generic analysis of all panel data for this clinseq pipeline.
        """

        # Configure alignment and merging for each unique sample library capture:
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

        # String indicating the sample, for use in output file names below:
        sample_str = "{}-{}".format(sample_type, sample_id)
        merged_bam_filename = "{}/bams/panel/{}-{}-{}.bam".format(self.outdir,
                                                                  sample_str,
                                                                  prep_kit_id,
                                                                  capture_kit_id)

        merge_bams = PicardMergeSamFiles(input_bams, merged_bam_filename)
        merge_bams.is_intermediate = True
        merge_bams.jobname = "picard-mergebams-{}".format(sample_str)
        self.add(merge_bams)

        markdups = PicardMarkDuplicates(merge_bams.output_bam,
                                        output_bam="{}/bams/panel/{}-{}-{}-nodups.bam".format(
                                            self.outdir, sample_str, prep_kit_id, capture_kit_id),
                                        output_metrics="{}/qc/picard/panel/{}-{}-{}-markdups-metrics.txt".format(
                                            self.outdir, sample_str, prep_kit_id, capture_kit_id))

        markdups.is_intermediate = False
        self.qc_files.append(markdups.output_metrics)
        self.add(markdups)
        return markdups.output_bam


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

    all_clinseq_barcodes = get_all_clinseq_barcodes(pipeline)

    # Configure alignment jobs for all library capture items:
    clinseq_barcode_to_bamfile = {}
    for clinseq_barcode in filter(lambda clinseq_barcode: clinseq_barcode != None,
                                  non_normal_clinseq_barcodes + [normal_clinseq_barcode]):
        clinseq_barcode_to_bamfile[clinseq_barcode] = \
            align_library(pipeline,
                          fq1_files=find_fastqs(clinseq_barcode, pipeline.libdir)[0],
                          fq2_files=find_fastqs(clinseq_barcode, pipeline.libdir)[1],
                          lib=clinseq_barcode,
                          ref=pipeline.refdata['bwaIndex'],
                          outdir=pipeline.outdir + "/bams/panel",
                          maxcores=pipeline.maxcores)

    # Configure merging and duplicate removal for each unique tumor and cfDNA
    # (sample type, sample ID, capture kit) pairing:
    capture_tup_to_clinseq_barcodes = get_capture_tup_to_clinseq_barcodes(non_normal_clinseq_barcodes)
    capture_tup_to_merged_bam = dict.fromkeys(capture_tup_to_clinseq_barcodes.keys(),[])
    for capture_tup in capture_tup_to_merged_bam.keys():
        capture_tup_to_merged_bam[capture_tup] = \
            merge_and_rm_dup(pipeline,
                             capture_tup_to_clinseq_barcodes[capture_tup],
                             clinseq_barcode_to_bamfile)

    # If there is a normal library capture, then do additional core panel analyses
    # supported by that item:
    if normal_clinseq_barcode != None:
        # Configure germline calling:
        germline_vcf = pipeline.call_germline_variants(\
            clinseq_barcode_to_bamfile[normal_clinseq_barcode], normal_clinseq_barcode)

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
