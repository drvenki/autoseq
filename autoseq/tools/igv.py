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


class MakeCNVkitTracks(Job):
    """
    Generate a IGV tracks representing the profile and segment information from a CNV-kit run.
    """

    def __init__(self):
        Job.__init__(self)
        self.input_cns = None
        self.input_cnr = None
        self.output_profile_bedgraph = None
        self.output_segments_bedgraph = None
        self.jobname = "make_cnvkit_tracks"

    def command(self):
        awk_cmd1 = "awk '$1 != \"chromosome\" \{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5\}' {} > {}".format(
            self.input_cnr, self.output_profile_bedgraph)
        awk_cmd2 = "awk '$1 != \"chromosome\" \{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5\}' {} > {}".format(
            self.input_cns, self.output_segments_bedgraph)
        return "{} && {}".format(awk_cmd1, awk_cmd2)
