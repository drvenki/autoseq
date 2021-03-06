from pypedream.job import Job, required, optional, repeat, conditional


class PicardCollectInsertSizeMetrics(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.output_metrics = None
        self.jobname = "picard-isize"

    def command(self):
        return "picard -XX:ParallelGCThreads=1 CollectInsertSizeMetrics H=/dev/null" + \
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
        return "picard -XX:ParallelGCThreads=1 -Xmx5g CollectGcBiasMetrics CHART=/dev/null" + \
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
        return "picard -XX:ParallelGCThreads=1 -Xmx2g CollectOxoGMetrics " + \
               required("I=", self.input) + \
               required("R=", self.reference_sequence) + \
               required("O=", self.output_metrics)


class PicardCollectHsMetrics(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.reference_sequence = None
        self.target_regions = None
        self.bait_regions = None
        self.bait_name = None
        self.output_metrics = None
        self.accumulation_level = ['LIBRARY']
        self.jobname = "picard-hsmetrics"

    def command(self):
        return "picard -XX:ParallelGCThreads=1 CollectHsMetrics " + \
               required("I=", self.input) + \
               required("R=", self.reference_sequence) + \
               required("O=", self.output_metrics) + \
               required("TI=", self.target_regions) + \
               required("BI=", self.bait_regions) + \
               optional("BAIT_SET_NAME=", self.bait_name) + \
               repeat('METRIC_ACCUMULATION_LEVEL=', self.accumulation_level)


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
        return "picard -XX:ParallelGCThreads=1 CollectWgsMetrics " + \
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
        return "picard -XX:ParallelGCThreads=1 CreateSequenceDictionary " + \
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
        return "picard -XX:ParallelGCThreads=1 BedToIntervalList " + \
               required("INPUT=", self.input) + \
               required("SEQUENCE_DICTIONARY=", self.reference_dict) + \
               required("OUTPUT=", self.output)


class PicardMergeSamFiles(Job):
    def __init__(self, input_bams, output_bam, assume_sorted=True, merge_dicts=True):
        Job.__init__(self)
        self.input_bams = input_bams
        self.output_bam = output_bam
        self.assume_sorted = assume_sorted
        self.merge_dicts = merge_dicts
        self.jobname = "picard-mergesamfiles"

    def command(self):
        return "picard -XX:ParallelGCThreads=1 MergeSamFiles " + \
               repeat("INPUT=", self.input_bams) + \
               required("ASSUME_SORTED=", str(self.assume_sorted).lower()) + \
               required("MERGE_SEQUENCE_DICTIONARIES=", str(self.merge_dicts).lower()) + \
               required("OUTPUT=", self.output_bam) + \
               " && samtools index " + required("", self.output_bam)


class PicardMarkDuplicates(Job):
    def __init__(self, input_bam, output_bam, output_metrics, remove_duplicates=False):
        Job.__init__(self)
        self.input_bam = input_bam
        self.output_bam = output_bam
        self.output_metrics = output_metrics
        self.remove_duplicates = remove_duplicates
        self.jobname = "picard-markdups"

    def command(self):
        return "picard -Xmx5g -XX:ParallelGCThreads=1 " + \
            required("-Djava.io.tmpdir=", self.scratch) + \
                " MarkDuplicates " + \
                required("INPUT=", self.input_bam) + \
                required("METRICS_FILE=", self.output_metrics) + \
                required("OUTPUT=", self.output_bam) + \
                conditional(self.remove_duplicates, "REMOVE_DUPLICATES=true") + \
                " && samtools index " + required("", self.output_bam)
