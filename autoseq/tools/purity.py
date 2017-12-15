from pypedream.job import Job, required, optional, conditional


# FIXME: Could be adapted so it's possible to run directly from bam files, instead of using pre-calculated segments

class PureCN(Job):
    def __init__(self, input_seg=None, input_vcf=None, tumorid=None, outdir=None, output=None, gcgene_file=None,
                 minpurity=0.05, hzdev=0.1, seg_sdev=None, genome="hg19", postopt=False):
        Job.__init__(self)
        self.input_seg = input_seg
        self.input_vcf = input_vcf
        self.tumorid = tumorid  # should be the same as the ID in the seg file
        self.outdir = outdir
        self.output = output
        self.gcgene_file = gcgene_file
        self.minpurity = minpurity
        self.hzdev = hzdev
        self.seg_sdev = seg_sdev
        self.genome = genome
        self.postopt = postopt

    def command(self):

        # activating conda env
        activate_cmd = "source activate purecn-env"

        # Determining output here purely as a means of determining when the job has
        # completed; PureCN itself determines the output file paths based on the
        # specified "--out" and "--sampleid" arguments.
        self.output = "{}/{}_genes.csv".format(
            self.outdir, self.tumorid)

        # running PureCN
        running_cmd = "Rscript PureCN.R " + required("--out ", self.outdir) + \
                       required("--sampleid ", self.tumorid) + \
                       required("--segfile ", self.input_seg) + \
                       required("--vcf ", self.input_vcf) + \
                       required("--gcgene ", self.gcgene_file) + \
                       required("--genome ", self.genome) + \
                       optional("--minpurity ", self.minpurity) + \
                       optional("--hzdev ", self.hzdev) + \
                       optional("--segfilesdev ", self.seg_sdev) + \
                       conditional(self.postopt, "--postoptimize")

        # deactivating the conda env
        deactivate_cmd = "source deactivate"

        return " && ".join([activate_cmd, running_cmd, deactivate_cmd])
