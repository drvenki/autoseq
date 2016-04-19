import logging
import uuid

from pypedream.job import Job, repeat, required, optional, conditional


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


class AlasccaGenomePlot(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_segments = None
        self.output_jpg = None
        self.chrsizes = None
        self.jobname = "alascca-genomeplot"

    def command(self):
        return "alascca_genomeplot.R " + \
               required("--input ", self.input_segments) + \
               required("--chrsizes ", self.chrsizes) + \
               required("--output ", self.output_jpg)


