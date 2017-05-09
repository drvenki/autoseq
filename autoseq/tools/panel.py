import collections

from autoseq.util.library import find_fastqs, parse_sdid, parse_sample_type, parse_sample_id, parse_capture_kit_id, parse_prep_kit_id
from autoseq.tools.alignment import align_library
from autoseq.tools.picard import PicardMergeSamFiles, PicardMarkDuplicates
from autoseq.tools.variantcalling import VcfAddSample, call_somatic_variants
from autoseq.tools.intervals import MsiSensor
from autoseq.tools.qc import *


def get_non_normal_clinseq_barcodes(pipeline):
    return filter(lambda bc: bc != None,
                  [pipeline.sampledata['panel']['T']]
                  + pipeline.sampledata['panel']['CFDNA'])


def get_capture_tup_to_clinseq_barcodes(clinseq_barcodes):
    capture_tup_to_barcodes = collections.defaultdict(list)
    for clinseq_barcode in clinseq_barcodes:
        curr_tup = (parse_sample_type(clinseq_barcode), parse_sample_id(clinseq_barcode),
                    parse_capture_kit_id(clinseq_barcode))
        capture_tup_to_barcodes[curr_tup].append(clinseq_barcode)

    return capture_tup_to_barcodes


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


def merge_and_rm_dup(pipeline, clinseq_barcodes, barcode_to_bams):
    """
    Configures Picard merging and duplicate marking, for the specified group of clinseq barcodes,
    with bam file names to be retrieved from the specified dictionary. The specified clinseq_barcodes
    must indicate the same sample type, sample ID, and capture kit ID, and be of length >= 1.

    :param clinseq_barcodes: List of clinseq barcodes for which to configure the analysis 
    :param barcode_to_bams: Dictionary of clinseq_barcode -> bam_filename key, value pairs
    :return: The filename of the output merged bam
    """
    
    capture_tups = [parse_capture_tuple(clinseq_barcode) for clinseq_barcode in clinseq_barcodes]
    if len(set(capture_tups)) != 1:
        raise ValueError("Invalid input clinseq_barcodes: " + clinseq_barcodes)

    (sample_type, sample_id, capture_kit_id, prep_kit_id) = \
        (capture_tups[0][0], capture_tups[0][1], capture_tups[0][2], capture_tups[0][3])

    sample_str = "{}-{}".format(sample_type, sample_id)

    input_bams = [barcode_to_bams[clinseq_barcode] for clinseq_barcode in clinseq_barcodes]

    merged_bam_filename = "{}/bams/panel/{}-{}-{}.bam".format(pipeline.outdir,
                                                              sample_str,
                                                              prep_kit_id,
                                                              capture_kit_id)

    merge_bams = PicardMergeSamFiles(input_bams, merged_bam_filename)
    merge_bams.is_intermediate = True
    merge_bams.jobname = "picard-mergebams-{}".format(sample_str)
    pipeline.add(merge_bams)

    markdups = PicardMarkDuplicates(merge_bams.output_bam,
                                    output_bam="{}/bams/panel/{}-{}-{}-nodups.bam".format(
                                        pipeline.outdir, sample_str, prep_kit_id, capture_kit_id),
                                    output_metrics="{}/qc/picard/panel/{}-{}-{}-markdups-metrics.txt".format(
                                        pipeline.outdir, sample_str, prep_kit_id, capture_kit_id))

    markdups.is_intermediate = False
    pipeline.qc_files.append(markdups.output_metrics)
    pipeline.add(markdups)
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
    :return: A dictionary with (sample type, sample barcode, capture kit code) tuples as keys, and bam filename lists
    as keys.
    """

    non_normal_clinseq_barcodes = get_non_normal_clinseq_barcodes(pipeline)
    normal_clinseq_barcode = pipeline.sampledata['panel']['N']

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
            merge_and_rm_dup(capture_tup_to_clinseq_barcodes[capture_tup],
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
    # get_capture_tup_to_clinseq_barcodes(normal_clinseq_barcode)
    # capture_tup_to_merged_bam[]


def analyze_panel_cancer_vs_normal(pipeline, sample_type, sample_id, capture_kit_id,
                                   cancer_bam, normal_bam, germline_vcf, normal_clinseq_barcode):
    """
    Configures several core analyses for a cancer vs normal comparison with library panel capture data.
    
    :param pipeline: 
    :param sample_type: 
    :param sample_id: 
    :param capture_kit_id: XXX
    :param cancer_bam: A bam file of the cancer reads, potentially merged
    :param normal_bam: 
    :param germline_vcf: 
    :param normal_clinseq_barcode: 
    :return: 
    """

    vep = False
    if pipeline.refdata['vep_dir']:
        vep = True

    sample_str = "{}-{}".format(sample_type, sample_id)

    # XXXFIXME: Need to figure out how to specify target_name below (i.e. the full name)

    # Configure somatic variant calling:
    somatic_variants = call_somatic_variants(pipeline, tbam=cancer_bam, nbam=normal_bam, tlib=sample_str,
                                             nlib=normal_clinseq_barcode, target_name=XXXFIXTHIS,
                                             refdata=pipeline.refdata, outdir=pipeline.outdir,
                                             callers=['vardict'], vep=vep)

    # Configure VCF add sample:
    vcfaddsample = VcfAddSample()
    vcfaddsample.input_bam = cancer_bam
    vcfaddsample.input_vcf = germline_vcf
    vcfaddsample.samplename = sample_str
    vcfaddsample.filter_hom = True
    vcfaddsample.output = "{}/variants/{}-and-{}.germline-variants-with-somatic-afs.vcf.gz".format(\
        pipeline.outdir, normal_clinseq_barcode, sample_str)
    vcfaddsample.jobname = "vcf-add-sample-{}".format(sample_str)
    pipeline.add(vcfaddsample)

    # Configure MSI sensor:
    msisensor = MsiSensor()
    msisensor.msi_sites = pipeline.refdata['targets'][targets_long]['msisites']
    msisensor.input_normal_bam = normal_bam
    msisensor.input_tumor_bam = cancer_bam
    msisensor.output = "{}/msisensor.tsv".format(pipeline.outdir)
    msisensor.threads = pipeline.maxcores
    msisensor.jobname = "msisensor-{}".format(sample_str)
    
    pipeline.add(msisensor)

    hzconcordance = HeterzygoteConcordance()
    hzconcordance.input_vcf = germline_vcf
    hzconcordance.input_bam = cancer_bam
    hzconcordance.reference_sequence = pipeline.refdata['reference_genome']
    hzconcordance.target_regions = pipeline.refdata['targets'][targets_long]['targets-interval_list-slopped20']
    hzconcordance.normalid = "{}-{}-{}".format(parse_sdid(normal_clinseq_barcode),
                                               parse_sample_type(normal_clinseq_barcode),
                                               parse_sample_id(normal_clinseq_barcode))
    hzconcordance.filter_reads_with_N_cigar = True
    hzconcordance.jobname = "hzconcordance-{}".format(sample_str)
    hzconcordance.output = "{}/bams/{}-{}-hzconcordance.txt".format(pipeline.outdir, sample_str, normal_clinseq_barcode)
    pipeline.add(hzconcordance)