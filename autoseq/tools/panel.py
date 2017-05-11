import collections

from autoseq.util.library import find_fastqs, parse_sdid, parse_sample_type, parse_sample_id, parse_capture_kit_id, parse_prep_kit_id
from autoseq.tools.variantcalling import VcfAddSample, call_somatic_variants
from autoseq.tools.intervals import MsiSensor
from autoseq.tools.qc import *


def contrast_cancer_vs_normal_panel(pipeline, sample_type, sample_id, capture_kit_id,
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