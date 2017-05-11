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

