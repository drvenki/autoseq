import logging
import os
import uuid

from pypedream.job import Job, repeat, required, optional, conditional, stripsuffix


class QDNASeq(Job):
    def __init__(self, input_bam, output_segments, background=None):
        Job.__init__(self)
        self.input = input_bam
        self.output = output_segments
        self.background = background
        self.jobname = "qdnaseq"

    def command(self):
        qdnaseq_cmd = "qdnaseq.R " + \
                      required("--bam ", self.input) + \
                      required("--output ", self.output) + \
                      optional("--background ", self.background)

        return qdnaseq_cmd


class QDNASeq2Bed(Job):
    def __init__(self, input_segments, output_bed, genes_gtf):
        Job.__init__(self)
        self.input_segments = input_segments
        self.output_bed = output_bed
        self.genes_gtf = genes_gtf

    def command(self):
        qdnaseq2bed_cmd = "qdnaseq2bed.py -n segments " + \
                          required("-i ", self.input_segments) + \
                          "| sort -k1,1 -k2,2n " + \
                          "| bedtools median -c 5 -o mean " + \
                          required("-a ", self.genes_gtf) + " -b - " + \
                          "| cnvgtf2bed.py -i /dev/stdin -n gene_id " + \
                          required("> ", self.output_bed)
        return qdnaseq2bed_cmd


class AlasccaCNAPlot(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_cnr = None
        self.input_cns = None
        self.input_germline_vcf = None
        self.input_somatic_vcf = None
        self.chrsizes = None
        self.output_png = None
        self.output_json = None
        self.jobname = "alascca-cna"

    def command(self):
        return "alasccaCNA.R " + \
               required("--cnr ", self.input_cnr) + \
               required("--cns ", self.input_cns) + \
               required("--germlinevcf ", self.input_germline_vcf) + \
               required("--somaticvcf ", self.input_somatic_vcf) + \
               required("--chrsizes ", self.chrsizes) + \
               required("--png ", self.output_png) + \
               required("--json ", self.output_json)


class CNVkit(Job):
    """Runs CNVkit. Either reference or targets_bed must be supplied"""

    def __init__(self, input_bam, output_cns, output_cnr, reference=None, targets_bed=None, scratch="/tmp"):
        self.input_bam = input_bam
        self.reference = reference
        self.output_cnr = output_cnr
        self.output_cns = output_cns
        self.targets_bed = targets_bed
        self.scratch = scratch

    def command(self):
        if not self.reference and not self.targets_bed:
            raise ValueError("Either reference or targets_bed must be supplied")
        if self.reference and self.targets_bed:
            raise ValueError("Supply either reference OR targets_bed")

        tmpdir = "{}/cnvkit-{}".format(self.scratch, uuid.uuid4())
        sample_prefix = stripsuffix(os.path.basename(self.input_bam), ".bam")
        cnvkit_cmd = "cnvkit.py batch " + required("", self.input_bam) + \
                     optional("-r ", self.reference) + \
                     conditional(self.targets_bed, "-n") + \
                     optional("-t ", self.targets_bed) + \
                     required("-d ", tmpdir)
        copy_cns_cmd = "cp {}/{}.cns ".format(tmpdir, sample_prefix) + required(" ", self.output_cns)
        copy_cnr_cmd = "cp {}/{}.cnr ".format(tmpdir, sample_prefix) + required(" ", self.output_cnr)
        rm_cmd = "rm -r {}".format(tmpdir)
        return " && ".join([cnvkit_cmd, copy_cns_cmd, copy_cnr_cmd, rm_cmd])
