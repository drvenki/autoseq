from pypedream.job import required, Job, conditional, optional

__author__ = 'Thomas Whitington'


class MakeAllelicFractionTrack(Job):
    """
    Generate an IGV track representing variant allelic fraction information (bedGraph file).
    """

    def __init__(self):
        Job.__init__(self)
        self.input_vcf = None
        self.output_bedgraph = None
        self.jobname = "make_allelic_fraction_track"

    def command(self):        
        return "generate_allelic_fraction_bedGraph.py " + \
               required("--output ", self.output_bedgraph) + \
               required(" ", self.input_vcf)
