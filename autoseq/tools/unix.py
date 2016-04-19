from pypedream.job import Job, required
# case class bwaIndex(ref:File) extends ExternalCommonArgs with  SingleCoreJob with OneDayJob {
#   @Input(doc="Input reference") val _ref: File = ref
#   @Output(doc="Input reference") val _out: File = ref + ".bwt"
#   def commandLine = bwaPath + " index -a bwtsw " + ref
#   this.jobName = "bwaIndex"
#   this.analysisName = this.jobName
#   this.isIntermediate = false
# }


class Copy(Job):
    def __init__(self):
	Job.__init__(self)
	self.input = None
	self.output = None
	self.jobname = "copy"

    def command(self):
	return "cp " + \
	       required(" ", self.input) + \
	       required(" ", self.output)


class Gunzip(Job):
    def __init__(self):
	Job.__init__(self)
	self.input = None
	self.output = None
	self.jobname = "gunzip"

    def command(self):
	return "gzip -cd " + \
	       required(" ", self.input) + \
	       required(" > ", self.output)


class Curl(Job):
    def __init__(self):
	Job.__init__(self)
	self.remote = None
	self.output = None
	self.jobname = "curl"

    def command(self):
	return "curl " + \
	       required(" ", self.remote) + \
	       required(" > ", self.output)
