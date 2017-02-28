from pypedream.job import required, Job, conditional, optional

__author__ = 'rebber'

class ContEst(Job):
    """Runs ContEst to estimate contamination level in bam file "input_eval_bam"."""

    def __init__(self):
        Job.__init__(self)
        self.reference_genome = None
        self.input_eval_bam = None
        self.input_genotype_bam = None
        self.population_af_vcf = None
        self.output = None
        self.jobname = "contest"

    def command(self):
        min_genotype_ratio = "0.95"

        return "java -Xmx15g -jar GenomeAnalysisTK.jar -T ContEst " + \
            required("-R ", self.reference_genome) + \
            required("-I:eval ", self.input_eval_bam) + \
            required("-I:genotype ", self.input_genotype_bam) + \
            required("--popfile ", self.population_af_vcf) + \
            required("--min_genotype_ratio ", min_genotype_ratio) + \
            required("-o ", self.output)


class ContEstToContamCaveat(Job):
    """Runs script to convert ContEst output to JSON file with contamination QC
    estimatimate."""

    def __init__(self):
        Job.__init__(self)
        self.input_contest_results = None
        self.output = None
        self.jobname = "contest"

    def command(self):
        return "contest_to_contam_caveat.py " + \
            required(" ", self.input_contest_results) + \
            required("> ", self.output)