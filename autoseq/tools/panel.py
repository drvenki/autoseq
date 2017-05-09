import collections

from autoseq.util.library import find_fastqs, parse_capture_kit_id, parse_sample_type
from autoseq.tools.alignment import align_library
from autoseq.tools.picard import PicardMergeSamFiles, PicardMarkDuplicates
from autoseq.tools.variantcalling import VcfAddSample, call_somatic_variants
from autoseq.tools.intervals import MsiSensor
from autoseq.tools.qc import *


def get_non_normal_clinseq_barcodes(pipeline):
    return filter(lambda bc: bc != None,
                  [pipeline.sampledata['panel']['T']]
                  + pipeline.sampledata['panel']['CFDNA'])


def analyze_panel(pipeline):
    """
Configures core analysis of all panel-captured libraries. The panel-capture libraries can include:
- Zero or one normal library captures
- Zero or one tumor library captures
- Zero or more cfDNA library captures

Performs the following core analyses:
- Configures alignment for all normal, cfDNA, and tumor library captures
- If germline is present, then it configures some additional downstream analyses:
-- Germline calling
-- Merging and duplicate removal for all bam files for each of the tumor and cfDNA library captures,
followed by a core non_normal vs normal analysis.

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

    # If there is a normal library capture, then do additional core panel analyses
    # supported by that item:
    if normal_clinseq_barcode != None:
        # Configure germline calling:
        germline_vcf = pipeline.call_germline_variants(\
            clinseq_barcode_to_bamfile[normal_clinseq_barcode], normal_clinseq_barcode)

        # For each unique cfDNA and tumor *sample* (i.e. uniquifying on the sample IDs):
        # XXX CONTINUE HERE: GENERATE DICTIONARY OF UNIQUE SAMPLE IDS, AND MODIFY/IMPLEMENT
        # THE CODE BELOW.
        
        for clinseq_barcode in non_normal_clinseq_barcodes:
            # Configure merging and duplicate removal for all bams for this library capture:
            merge_and_rm_dup(clinseq_barcode, XXX)
    
            # Configure core analysis with the resulting bam file together with the normal
            # bam file and the germline VCF:
            XXX
        

    # process tumor and plasma samples
    libs = [x for x in plibs + [tlib] if x is not None]
    sample_bams = collections.defaultdict(list)
    sample_dicts = collections.defaultdict(list)
    for lib in libs:
        libdict = get_libdict(lib)
        sample = "{}-{}-{}-{}".format(libdict['sdid'], libdict['type'], libdict['sample_id'],
                                      parse_capture_kit_id(lib))
        bam = align_library(pipeline,
                            fq1_files=find_fastqs(lib, pipeline.libdir)[0],
                            fq2_files=find_fastqs(lib, pipeline.libdir)[1],
                            lib=lib,
                            ref=pipeline.refdata['bwaIndex'],
                            outdir=pipeline.outdir + "/bams/panel",
                            maxcores=pipeline.maxcores,
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
                                             pipeline.outdir, sample, prep_kit_short, targets_short))

        merge_bams.is_intermediate = True
        merge_bams.jobname = "picard-mergebams-{}".format(sample)
        pipeline.add(merge_bams)

        markdups = PicardMarkDuplicates(merge_bams.output_bam,
                                        output_bam="{}/bams/panel/{}-{}-{}-nodups.bam".format(
                                            pipeline.outdir, sample, prep_kit_short, targets_short),
                                        output_metrics="{}/qc/picard/panel/{}-{}-{}-markdups-metrics.txt".format(
                                            pipeline.outdir, sample, prep_kit_short, targets_short))
        markdups.is_intermediate = False
        pipeline.qc_files.append(markdups.output_metrics)
        pipeline.add(markdups)

        if set([smp['type'] for smp in sample_dicts[sample]]) == set('CFDNA'):
            # if it's a plasma sample
            pbams.append(markdups.output_bam)
        elif set([smp['type'] for smp in sample_dicts[sample]]) == set('T'):
            # if it's a tumor sample
            tbam = markdups.output_bam

        vep = False
        if pipeline.refdata['vep_dir']:
            vep = True

        if nlib:
            somatic_variants = call_somatic_variants(pipeline, tbam=markdups.output_bam, nbam=nbam, tlib=sample,
                                                     nlib=nlib, target_name=targets_long, refdata=pipeline.refdata,
                                                     outdir=pipeline.outdir,
                                                     callers=['vardict'],
                                                     vep=vep)

            vcfaddsample = VcfAddSample()
            vcfaddsample.input_bam = markdups.output_bam
            vcfaddsample.input_vcf = germline_vcf
            vcfaddsample.samplename = sample
            vcfaddsample.filter_hom = True
            vcfaddsample.output = "{}/variants/{}-and-{}.germline-variants-with-somatic-afs.vcf.gz".format(
                pipeline.outdir,
                nlib,
                sample)
            vcfaddsample.jobname = "vcf-add-sample-{}".format(sample)
            pipeline.add(vcfaddsample)

            msisensor = MsiSensor()
            msisensor.msi_sites = pipeline.refdata['targets'][targets_long]['msisites']
            msisensor.input_normal_bam = nbam
            msisensor.input_tumor_bam = markdups.output_bam
            msisensor.output = "{}/msisensor.tsv".format(pipeline.outdir)
            msisensor.threads = pipeline.maxcores
            msisensor.jobname = "msisensor-{}".format(sample)
            pipeline.add(msisensor)

            libdict = get_libdict(nlib)
            rg_sm = "{}-{}-{}".format(libdict['sdid'], libdict['type'], libdict['sample_id'])

            hzconcordance = HeterzygoteConcordance()
            hzconcordance.input_vcf = germline_vcf
            hzconcordance.input_bam = markdups.output_bam
            hzconcordance.reference_sequence = pipeline.refdata['reference_genome']
            hzconcordance.target_regions = pipeline.refdata['targets'][targets_long]['targets-interval_list-slopped20']
            hzconcordance.normalid = rg_sm
            hzconcordance.filter_reads_with_N_cigar = True
            hzconcordance.jobname = "hzconcordance-{}".format(lib)
            hzconcordance.output = "{}/bams/{}-{}-hzconcordance.txt".format(pipeline.outdir, lib, nlib)
            pipeline.add(hzconcordance)

    return {'tbam': tbam, 'nbam': nbam, 'pbams': pbams,
            'somatic_variants': somatic_variants}
