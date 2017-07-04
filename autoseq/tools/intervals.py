import uuid
from pypedream.job import required, Job, conditional


class SlopIntervalList(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.output = None
        self.jobname = "slop-interval-list"

    def command(self):
        return "slopIntervalList.py " + \
               required(" < ", self.input) + \
               required(" > ", self.output)


class IntervalListToBed(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.output = None
        self.jobname = "interval-list-to-bed"

    def command(self):
        return "picard_interval_list_to_bed6_converter.py " + \
               required(" ", self.input) + \
               required(" ", self.output)


class MsiSensorScan(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_fasta = None
        self.output = None
        self.homopolymers_only = True
        self.jobname = "msisensor-scan"

    def command(self):
        return "msisensor scan " + \
               required("-d ", self.input_fasta) + \
               required("-o ", self.output) + \
               conditional(self.homopolymers_only, " -p 1 ")


class IntersectMsiSites(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_msi_sites = None
        self.target_bed = None
        self.output_msi_sites = True
        self.jobname = "msi-intersect"

    def command(self):
        return "intersect-msi-sites.sh " + \
               required(" ", self.input_msi_sites) + \
               required(" ", self.target_bed) + \
               required(" ", self.output_msi_sites)


class MsiSensor(Job):
    def __init__(self):
        Job.__init__(self)
        self.msi_sites = None
        self.input_normal_bam = None
        self.input_tumor_bam = None
        self.output = None
        self.jobname = "msisensor"

    def command(self):
        output_prefix = "{scratch}/msisensor-{uuid}".format(scratch=self.scratch,
                                                            uuid=uuid.uuid4())
        output_table = "{}".format(output_prefix)
        output_dis = "{}_dis".format(output_prefix)
        output_germline = "{}_germline".format(output_prefix)
        output_somatic = "{}_somatic".format(output_prefix)

        return "msisensor msi " + \
               required("-d ", self.msi_sites) + \
               required("-n ", self.input_normal_bam) + \
               required("-t ", self.input_tumor_bam) + \
               required("-o ", output_prefix) + \
               required("-b ", self.threads) + \
               " && cp {} {}".format(output_prefix, self.output) + \
               " && rm {} {} {} {}".format(output_table, output_dis,
                                           output_germline, output_somatic)
