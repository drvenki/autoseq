from pypedream.job import Job, required, optional


class PicardCollectInsertSizeMetrics(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.output_metrics = None
        self.jobname = "picard-isize"

    def command(self):
        return "picard CollectInsertSizeMetrics H=/dev/null" + \
               required("I=", self.input) + \
               required("O=", self.output_metrics)


class PicardCollectGcBiasMetrics(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.reference_sequence = None
        self.output_metrics = None
        self.output_summary = None
        self.stop_after = None
        self.jobname = "picard-gcbias"

    def command(self):
        return "picard CollectGcBiasMetrics CHART=/dev/null" + \
               required("I=", self.input) + \
               required("O=", self.output_metrics) + \
               required("S=", self.output_summary) + \
               required("R=", self.reference_sequence) + \
               optional("STOP_AFTER=", self.stop_after)


class PicardCollectOxoGMetrics(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.reference_sequence = None
        self.output_metrics = None
        self.jobname = "picard-oxog"

    def command(self):
        return "picard -Xmx2g CollectOxoGMetrics " + \
               required("I=", self.input) + \
               required("R=", self.reference_sequence) + \
               required("O=", self.output_metrics)


class PicardCalculateHsMetrics(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.reference_sequence = None
        self.target_regions = None
        self.bait_regions = None
        self.bait_name = None
        self.output_metrics = None
        self.jobname = "picard-hsmetrics"

    def command(self):
        return "picard CalculateHsMetrics " + \
               required("I=", self.input) + \
               required("R=", self.reference_sequence) + \
               required("O=", self.output_metrics) + \
               required("TI=", self.target_regions) + \
               required("BI=", self.bait_regions) + \
               optional("BAIT_SET_NAME=", self.bait_name)


class PicardCollectWgsMetrics(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.reference_sequence = None
        self.minimum_mapping_quality = None
        self.minimum_base_quality = None
        self.coverage_cap = None
        self.output_metrics = None
        self.jobname = "picard-wgsmetrics"

    def command(self):
        return "picard CollectWgsMetrics " + \
               required("I=", self.input) + \
               required("R=", self.reference_sequence) + \
               required("O=", self.output_metrics) + \
               optional("MINIMUM_MAPPING_QUALITY=", self.minimum_mapping_quality) + \
               optional("MINIMUM_BASE_QUALITY=", self.minimum_base_quality) + \
               optional("COVERAGE_CAP=", self.coverage_cap)


class PicardCreateSequenceDictionary(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.output_dict = None
        self.jobname = "picard-createdict"

    def command(self):
        return "picard CreateSequenceDictionary " + \
               required("REFERENCE=", self.input) + \
               required("OUTPUT=", self.output_dict)


class PicardBedToIntervalList(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.reference_dict = None
        self.output = None
        self.jobname = "picard-bedtointervallist"

    def command(self):
        return "picard BedToIntervalList " + \
               required("INPUT=", self.input) + \
               required("SEQUENCE_DICTIONARY=", self.reference_dict) + \
               required("OUTPUT=", self.output)

