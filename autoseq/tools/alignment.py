from pypedream.job import *

__author__ = 'dankle'


class Bwa(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_fastq1 = None  # input ports must start with "input"
        self.input_fastq2 = None
        self.input_reference_sequence = None
        self.mark_secondary = None
        self.remove_duplicates = True
        self.readgroup = None
        self.output = None  # output ports must start with "output", can be "output_metrics", "output", etc
        self.duplication_metrics = None
        self.jobname = "bwa"

    def command(self):
        bwalog = self.output + ".bwa.log"
        samblasterlog = self.output + ".samblaster.log"
        tmpprefix = "{}/{}".format(self.scratch, uuid.uuid4())
        return "bwa mem -M -v 1 " + \
               required("-R ", self.readgroup) + \
               optional("-t ", self.threads) + \
               required(" ", self.input_reference_sequence) + \
               required(" ", self.input_fastq1) + \
               optional("", self.input_fastq2) + \
               required("2>", bwalog) + \
               "| samblaster -M --addMateTags " + \
               conditional(self.remove_duplicates, "--removeDups") + \
               optional("--metricsFile ", self.duplication_metrics) + \
               required("2>", samblasterlog) + \
               "| samtools view -Sb -u - " + \
               "| samtools sort " + \
               required("-T ", tmpprefix) + \
               optional("-@ ", self.threads) + \
               required("-o ", self.output) + \
               " - " + \
               " && samtools index " + self.output + \
               " && cat {} {} ".format(bwalog, samblasterlog) + \
               " && rm {} {} ".format(bwalog, samblasterlog)


class SkewerPE(Job):
    def __init__(self):
        Job.__init__(self)
        self.input1 = ""
        self.input2 = ""
        self.output1 = ""
        self.output2 = ""
        self.stats = ""
        self.jobname = "skewer"

    def command(self):
        tmpdir = "{}/skewer-{}".format(self.scratch, uuid.uuid4())
        prefix = "{}/skewer".format(tmpdir)
        out_fq1 = prefix + "-trimmed-pair1.fastq.gz"
        out_fq2 = prefix + "-trimmed-pair2.fastq.gz"
        out_stats = prefix + "-trimmed.log"

        mkdir_cmd = "mkdir -p {}".format(tmpdir)
        skewer_cmd = "skewer -z " + \
                     optional("-t ", self.threads) + " --quiet " + \
                     required("-o ", prefix) + \
                     required(" ", self.input1) + \
                     required(" ", self.input2)
        copy_fq1_cmd = "cp " + out_fq1 + " " + self.output1
        copy_fq2_cmd = "cp " + out_fq2 + " " + self.output2
        copy_stats_cmd = "cp " + out_stats + " " + self.stats
        rm_cmd = "rm -r {}".format(tmpdir)
        return " && ".join([mkdir_cmd, skewer_cmd, copy_fq1_cmd, copy_fq2_cmd, copy_stats_cmd, rm_cmd])


class SkewerSE(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.output = None
        self.stats = None
        self.jobname = "skewer"

    def command(self):
        tmpdir = "{}/skewer-{}".format(self.scratch, uuid.uuid4())
        prefix = "{}/skewer".format(tmpdir)
        out_fq1 = prefix + "-trimmed.fastq.gz"
        out_stats = prefix + "-trimmed.log"

        mkdir_cmd = "mkdir -p {}".format(tmpdir)
        skewer_cmd = "skewer -z " + \
                     optional("-t ", self.threads) + " --quiet " + \
                     required("-o ", prefix) + \
                     required(" ", self.input)
        copy_fq1_cmd = "cp " + out_fq1 + " " + self.output
        copy_stats_cmd = "cp " + out_stats + " " + self.stats
        rm_cmd = "rm -r {}".format(tmpdir)
        return " && ".join([mkdir_cmd, skewer_cmd, copy_fq1_cmd, copy_stats_cmd, rm_cmd])


class CatAndSkewer(Job):
    def __init__(self):
        Job.__init__(self)
        self.input1 = ""
        self.input2 = ""
        self.output1 = ""
        self.output2 = ""
        self.stats = ""
        self.downsample = -1
        self.jobname = "skewer"

    def command(self):
        prefix = "{scratch}/skewer-{uuid}/skewer".format(scratch=self.scratch, uuid=uuid.uuid4())
        out_fq1 = prefix + "-pair1.fastq.gz"
        out_fq2 = prefix + "-pair2.fastq.gz"
        tmp_fq1 = prefix + "-input_1.fastq.gz"
        tmp_fq2 = prefix + "-input_2.fastq.gz"
        out_stats = prefix + ".log"
        head_cmd = ""
        if self.downsample > 0:
            head_cmd = "|gzip -cd|head -n " + 4 * self.downsample + "|gzip"
        mkdir_cmd = "mkdir -p {}".format(os.path.dirname(prefix))
        cat_fq1_cmd = "cat " + repeat("", self.input1) + head_cmd + " > " + tmp_fq1
        cat_fq2_cmd = "cat " + repeat("", self.input2) + head_cmd + " > " + tmp_fq2
        skewer_cmd = "skewer -z -t " + str(self.threads) + " --quiet -o " + prefix + " " + tmp_fq1 + " " + tmp_fq2
        copy_fq1_cmd = "cp " + out_fq1 + " " + self.output1
        copy_fq2_cmd = "cp " + out_fq2 + " " + self.output2
        copy_stats_cmd = "cp " + out_stats + " " + self.stats
        rm_cmd = "rm " + tmp_fq1 + " " + tmp_fq2
        return " && ".join(
            [mkdir_cmd, cat_fq1_cmd, cat_fq2_cmd, skewer_cmd, copy_fq1_cmd, copy_fq2_cmd, copy_stats_cmd, rm_cmd])


class Cutadapt(Job):
    def __init__(self):
        Job.__init__(self)
        self.input1 = ""
        self.input2 = ""
        self.output1 = ""
        self.output2 = ""
        self.scratch = "/tmp/"
        self.downsample = -1
        self.jobname = "cutadapt"

    def command(self):
        adapters = " -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT "
        prefix = "{scratch}/cutadapt-{uuid}/cutadapt".format(scratch=self.scratch, uuid=uuid.uuid4())
        out_fq1 = prefix + "-pair1.fastq.gz"
        out_fq2 = prefix + "-pair2.fastq.gz"
        tmp_fq1 = prefix + "-input_1.fastq.gz"
        tmp_fq2 = prefix + "-input_2.fastq.gz"
        head_cmd = ""
        if self.downsample > 0:
            head_cmd = "|gzip -cd|head -n " + 4 * self.downsample + "|gzip"
        mkdir_cmd = "mkdir -p {}".format(os.path.dirname(prefix))
        cat_fq1_cmd = "cat " + repeat("", self.input1) + head_cmd + " > " + tmp_fq1
        cat_fq2_cmd = "cat " + repeat("", self.input2) + head_cmd + " > " + tmp_fq2
        cutadapt_cmd = "cutadapt " + adapters + " -o " + out_fq1 + " -p " + out_fq2 + " " + tmp_fq1 + " " + tmp_fq2
        copy_fq1_cmd = "cp " + out_fq1 + " " + self.output1
        copy_fq2_cmd = "cp " + out_fq2 + " " + self.output2
        rm_cmd = "rm " + tmp_fq1 + " " + tmp_fq2
        return " && ".join([mkdir_cmd, cat_fq1_cmd, cat_fq2_cmd, cutadapt_cmd, copy_fq1_cmd, copy_fq2_cmd, rm_cmd])
