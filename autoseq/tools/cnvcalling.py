import logging
import os
import uuid

from pypedream.job import Job, repeat, required, optional, conditional, stripsuffix


class QDNASeq(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.output_segments = None
        self.output_bed = None
        self.background = None
        self.genes_gtf = None
        self.jobname = "qdnaseq"

    def command(self):
        qdnaseq_cmd = "qdnaseq.R " + \
                      required("--bam ", self.input) + \
                      required("--output ", self.output_segments) + \
                      required("--background ", self.background)

        qdnaseq2bed_cmd = "qdnaseq2bed.py -n call " + \
                          required("-i ", self.output_segments) + \
                          "| sort -k1,1 -k2,2n " + \
                          "| bedtools median -c 5 -o mean " + \
                          required("-a ", self.genes_gtf) + " -b - " + \
                          "| cnvgtf2bed.py -i /dev/stdin -n gene_id " + \
                          required("> ", self.output_bed)

        return " && ".join([qdnaseq_cmd, qdnaseq2bed_cmd])

        # python tools/scripts/qdnaseq2bed.py -i $QDNASEQ |
        # sort -k1,1 -k2,2n |
        # bedtools map -c 5 -o mean -a $GENESSORTED -b - |
        # python tools/scripts/cnvgtf2bed.py -i /dev/stdin -n gene_id > ~/bla.bed


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
    def __init__(self, input_bam, reference, output_cns, output_cnr, scratch="/tmp"):
        self.input_bam = input_bam
        self.reference = reference
        self.output_cnr = output_cnr
        self.output_cns = output_cns
        self.scratch = scratch

    def command(self):
        tmpdir = "{}/cnvkit-{}".format(self.scratch, uuid.uuid4())
        sample_prefix = stripsuffix(os.path.basename(self.input_bam), ".bam")
        cnvkit_cmd = "cnvkit.py batch " + required("", self.input_bam) + \
            required("-r ", self.reference) + \
                     required("-d ", tmpdir)
        copy_cns_cmd = "cp {}/{}.cns ".format(tmpdir, sample_prefix) + required(" ", self.output_cns)
        copy_cnr_cmd = "cp {}/{}.cnr ".format(tmpdir, sample_prefix) + required(" ", self.output_cnr)
        rm_cmd = "rm -r {}".format(tmpdir)
        return " && ".join([cnvkit_cmd, copy_cns_cmd, copy_cnr_cmd, rm_cmd])


