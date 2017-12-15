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
        self.output_cna = None
        self.output_purity = None
        self.jobname = "alascca-cna"

    def command(self):
        return "alasccaCNA.R " + \
               required("--cnr ", self.input_cnr) + \
               required("--cns ", self.input_cns) + \
               required("--germlinevcf ", self.input_germline_vcf) + \
               required("--somaticvcf ", self.input_somatic_vcf) + \
               required("--chrsizes ", self.chrsizes) + \
               required("--png ", self.output_png) + \
               required("--json.cna ", self.output_cna) + \
               required("--json.purity ", self.output_purity)


class LiqbioCNAPlot(Job):
    def __init__(self):
        Job.__init__(self)
        self.tumor_cnr = None
        self.tumor_cns = None
        self.normal_cnr = None
        self.normal_cns = None
        self.het_snps_vcf = None
        self.purecn_csv = None
        self.purecn_genes_csv = None
        self.purecn_loh_csv = None
        self.purecn_variants_csv = None
        self.svcaller_T_DEL = None
        self.svcaller_T_DUP = None
        self.svcaller_T_INV = None
        self.svcaller_T_TRA = None
        self.svcaller_N_DEL = None
        self.svcaller_N_DUP = None
        self.svcaller_N_INV = None
        self.svcaller_N_TRA = None
        self.germline_mut_vcf = None
        self.somatic_mut_vcf = None
        self.plot_png = None
        self.output_cna = None
        self.output_purity = None

        self.jobname = "liqbio-cna"

    def command(self):
        return "liqbioCNA.R" + \
        required("--tumor_cnr", self.tumor_cnr) + \
        required("--tumor_cns", self.tumor_cns) + \
        required("--normal_cnr", self.normal_cnr) + \
        required("--normal_cns", self.normal_cns) + \
        required("--het_snps_vcf", self.het_snps_vcf) + \
        required("--purecn_csv", self.purecn_csv) + \
        required("--purecn_genes_csv", self.purecn_genes_csv) + \
        required("--purecn_loh_csv", self.purecn_loh_csv) + \
        required("--purecn_variants_csv", self.purecn_variants_csv) + \
        required("--svcaller_T_DEL", self.svcaller_T_DEL) + \
        required("--svcaller_T_DUP", self.svcaller_T_DUP) + \
        required("--svcaller_T_INV", self.svcaller_T_INV) + \
        required("--svcaller_T_TRA", self.svcaller_T_TRA) + \
        required("--svcaller_N_DEL", self.svcaller_N_DEL) + \
        required("--svcaller_N_DUP", self.svcaller_N_DUP) + \
        required("--svcaller_N_INV", self.svcaller_N_INV) + \
        required("--svcaller_N_TRA", self.svcaller_N_TRA) + \
        required("--germline_mut_vcf", self.germline_mut_vcf) + \
        required("--somatic_mut_vcf", self.somatic_mut_vcf) + \
        required("--plot_png", self.plot_png) + \
        required("--cna_json", self.output_cna) + \
        required("--purity_json", self.output_purity)


class CNVkit(Job):
    """Runs CNVkit. Either reference or targets_bed must be supplied"""

    def __init__(self, input_bam, output_cns, output_cnr, reference=None,
                 targets_bed=None, scratch="/tmp", fasta=None):
        self.input_bam = input_bam
        self.reference = reference
        self.fasta = fasta
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
                     conditional(self.targets_bed, "--fasta " + str(self.fasta) + " --split") + \
                     conditional(self.targets_bed, "-n") + \
                     optional("-t ", self.targets_bed) + \
                     required("-d ", tmpdir)
        copy_cns_cmd = "cp {}/{}.cns ".format(tmpdir, sample_prefix) + required(" ", self.output_cns)
        copy_cnr_cmd = "cp {}/{}.cnr ".format(tmpdir, sample_prefix) + required(" ", self.output_cnr)
        rm_cmd = "rm -r {}".format(tmpdir)
        return " && ".join([cnvkit_cmd, copy_cns_cmd, copy_cnr_cmd, rm_cmd])


class CNVkitFix(Job):
    """Fixes the output of CNV-kit, using a specified table of reference data that should be
    specific to the sample type, capture kit, and library prep kit.
    """

    def __init__(self, input_cnr, input_cns, input_ref, output_cnr, output_cns):
        self.input_cnr = input_cnr
        self.input_cns = input_cns
        self.input_ref = input_ref
        self.output_cns = output_cns
        self.output_cnr = output_cnr

    def command(self):
        return "fix_cnvkit.py --input-cnr {input_cnr} --input-cns {input_cns} " + \
               "--input-reference {input_ref} --output-cnr {output_cnr} " + \
               "--output-cns {output_cns}".format(
                   input_cnr=self.input_cnr,
                   input_cns=self.input_cns,
                   input_ref=self.input_ref,
                   output_cns=self.output_cns,
                   output_cnr=self.output_cnr)


class Cns2Seg(Job):
    """Converts CNVkit segment format (.cns files) to DNAcopy segment format (.seg files). """

    def __init__(self, input_cns, output_seg):
        self.input_cns = input_cns
        self.output_seg = output_seg

    def command(self):
        cmd = "cnvkit.py export seg " + required("-o ", self.output_seg) + required("", self.input_cns)
        return cmd
