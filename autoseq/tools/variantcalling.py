import logging
import sys
import uuid

from pypedream.job import Job, repeat, required, optional, conditional
from autoseq.util.clinseq_barcode import *
from autoseq.util.vcfutils import vt_split_and_leftaln, fix_ambiguous_cl, remove_dup_cl


class Freebayes(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_bams = None
        self.reference_sequence = None
        self.target_bed = None
        self.somatic_only = False
        self.params = "--pooled-discrete --pooled-continuous --genotype-qualities --report-genotype-likelihood-max --allele-balance-priors-off"
        self.min_coverage = 20
        self.min_alt_frac = 0.01
        self.use_harmonic_indel_quals = False
        self.output = ""
        self.jobname = "freebayes-somatic"

    def command(self):
        regions_file = "{scratch}/{uuid}.regions".format(scratch=self.scratch, uuid=uuid.uuid4())
        bed_to_regions_cmd = "cat {} | bed_to_regions.py > {}".format(self.target_bed, regions_file)

        call_somatic_cmd = " | {} -c 'from autoseq.util.bcbio import call_somatic; import sys; print call_somatic(sys.stdin.read())' ".format(
            sys.executable)

        freebayes_cmd = "freebayes-parallel {} {} ".format(regions_file, self.threads) + \
                        required("-f ", self.reference_sequence) + " --use-mapping-quality " + \
                        optional("--min-alternate-fraction ", self.min_alt_frac) + \
                        optional("--min-coverage ", self.min_coverage) + \
                        conditional(self.use_harmonic_indel_quals, "--harmonic-indel-quality") + \
                        optional("", self.params) + \
                        repeat(" ", self.input_bams) + \
                        """| bcftools filter -i 'ALT="<*>" || QUAL > 5' """ + \
                        "| filter_erroneus_alt.py -V /dev/stdin " + \
                        conditional(self.somatic_only, call_somatic_cmd) + \
                        " | " + vt_split_and_leftaln(self.reference_sequence) + \
                        " | vcfuniq | bcftools view --apply-filters .,PASS " + \
                        " | bgzip > {output} && tabix -p vcf {output}".format(output=self.output)
        # reason for 'vcfuniq': freebayes sometimes report duplicate variants that need to be uniqified.
        rm_regions_cmd = "rm {}".format(regions_file)
        return " && ".join([bed_to_regions_cmd, freebayes_cmd, rm_regions_cmd])


class VarDict(Job):
    def __init__(self, input_tumor=None, input_normal=None, tumorid=None, normalid=None, reference_sequence=None,
                 reference_dict=None, target_bed=None, output=None, min_alt_frac=0.1):
        Job.__init__(self)
        self.input_tumor = input_tumor
        self.input_normal = input_normal
        self.tumorid = tumorid
        self.normalid = normalid
        self.reference_sequence = reference_sequence
        self.reference_dict = reference_dict
        self.target_bed = target_bed
        self.output = output
        self.min_alt_frac = min_alt_frac

    def command(self):
        required("", self.input_tumor)
        required("", self.input_normal)

        freq_filter = (" bcftools filter -e 'STATUS !~ \".*Somatic\"' 2> /dev/null "
                       "| %s -c 'from autoseq.util.bcbio import depth_freq_filter_input_stream; import sys; print depth_freq_filter_input_stream(sys.stdin, %s, \"%s\")' " %
                       (sys.executable, 0, 'bwa'))

        somatic_filter = (" sed 's/\\.*Somatic\\\"/Somatic/' "  # changes \".*Somatic\" to Somatic
                          "| sed 's/REJECT,Description=\".*\">/REJECT,Description=\"Not Somatic via VarDict\">/' "
                          "| %s -c 'from autoseq.util.bcbio import call_somatic; import sys; print call_somatic(sys.stdin.read())' " % sys.executable)

        cmd = "vardict-java " + required("-G ", self.reference_sequence) + \
              optional("-f ", self.min_alt_frac) + \
              required("-N ", self.tumorid) + \
              " -b \"{}|{}\" ".format(self.input_tumor, self.input_normal) + \
              " -c 1 -S 2 -E 3 -g 4 -Q 10 " + required("", self.target_bed) + \
              " | testsomatic.R " + \
              " | var2vcf_paired.pl -P 0.9 -m 4.25 -M " + required("-f ", self.min_alt_frac) + \
              " -N \"{}|{}\" ".format(self.tumorid, self.normalid) + \
              " | " + freq_filter + " | " + somatic_filter + " | " + fix_ambiguous_cl() + " | " + remove_dup_cl() + \
              " | vcfstreamsort -w 1000 " + \
              " | " + vt_split_and_leftaln(self.reference_sequence) + \
              " | bcftools view --apply-filters .,PASS " + \
              " | vcfsorter.pl {} /dev/stdin ".format(self.reference_dict) + \
              " | bgzip > {output} && tabix -p vcf {output}".format(output=self.output)
        return cmd


class VEP(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_vcf = None
        self.output_vcf = None
        self.reference_sequence = None
        self.vep_dir = None
        self.jobname = "vep"

    def command(self):
        bgzip = ""
        fork = ""
        if self.threads > 1 and self.vep_dir:  # vep does not accept "--fork 1", so need to check.
            fork = " --fork {} ".format(self.threads)
        if self.output_vcf.endswith('gz'):
            bgzip = " | bgzip "

        cmdstr = "variant_effect_predictor.pl --vcf --output_file STDOUT " + \
                 optional("--dir ", self.vep_dir) + \
                 required("--fasta ", self.reference_sequence) + \
                 required("-i ", self.input_vcf) + \
                 " --check_alleles --check_existing  --total_length --allele_number " + \
                 " --no_escape --no_stats " + \
                 conditional(self.vep_dir, " --everything ") + \
                 conditional(self.vep_dir, " --offline ") + \
                 conditional(not self.vep_dir, " --database ") + \
                 conditional(not self.vep_dir, " --port 3337 ") + \
                 conditional(not self.vep_dir, " --hgvs ") + \
                 fork + bgzip + " > " + required("", self.output_vcf) + \
                 " && tabix -p vcf {}".format(self.output_vcf)

        return cmdstr


class VcfAddSample(Job):
    """
    Add DP, RO and AO tags for a new sample to a VCF, filter low-qual variants on the fly
    """

    def __init__(self):
        Job.__init__(self)
        self.input_vcf = None
        self.input_bam = None
        self.samplename = None
        self.filter_hom = True
        self.output = None
        self.jobname = "vcf-add-sample"

    def command(self):
        filt_vcf = "{scratch}/{uuid}.vcf.gz".format(scratch=self.scratch, uuid=uuid.uuid4())
        bgzip = ""
        tabix = ""
        if self.output.endswith('gz'):
            bgzip = "| bgzip"
            tabix = " && tabix -p vcf {}".format(self.output)

        filt_vcf_cmd = "vcf_filter.py --no-filtered " + required("", self.input_vcf) + " sq --site-quality 5 " + \
                       "|bgzip" + " > " + filt_vcf
        vcf_add_sample_cmd = "vcf_add_sample.py " + \
                             conditional(self.filter_hom, "--filter_hom") + \
                             required("--samplename ", self.samplename) + \
                             filt_vcf + " " + \
                             required("", self.input_bam) + \
                             bgzip + " > " + self.output + tabix
        rm_filt_cmd = "rm " + filt_vcf
        return " && ".join([filt_vcf_cmd, vcf_add_sample_cmd, rm_filt_cmd])


class VcfFilter(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.filter = None
        self.output = None
        self.jobname = "vcffilter"

    def command(self):
        return "zcat" + \
               required(" ", self.input) + \
               "| vcffilter " + \
               required("-f ", self.filter) + \
               "| bgzip " + required(" > ", self.output) + \
               " && tabix -p vcf {output}".format(output=self.output)


class CurlSplitAndLeftAlign(Job):
    def __init__(self):
        Job.__init__(self)
        self.remote = None
        self.input_reference_sequence = None
        self.input_reference_sequence_fai = None
        self.output = None
        self.jobname = "curl-split-leftaln"

    def command(self):
        required("", self.input_reference_sequence_fai)
        return "curl -L " + \
               required(" ", self.remote) + \
               "| gzip -d |" + vt_split_and_leftaln(self.input_reference_sequence, allow_ref_mismatches=True) + \
               "| bgzip " + required(" > ", self.output) + \
               " && tabix -p vcf {output}".format(output=self.output)


class InstallVep(Job):
    def __init__(self):
        Job.__init__(self)
        self.output_dir = None
        self.jobname = "fetch-vep-cache"

    def command(self):
        return "vep_install.pl --SPECIES homo_sapiens_vep --AUTO c --ASSEMBLY GRCh37 --NO_HTSLIB " + \
               required("--CACHEDIR ", self.output_dir) + \
               " && vep_convert_cache.pl " + required("--dir ", self.output_dir) + \
               " --species homo_sapiens --version 83_GRCh37"


def call_somatic_variants(pipeline, tbam, nbam, tumor_capture, normal_capture, target_name, refdata, outdir,
                          callers=['vardict', 'freebayes'], vep=True, min_alt_frac=0.1):
    """
    Call somatic variants.
    :param pipeline:
    :param tbam:
    :param nbam:
    :param tlib:
    :param nlib:
    :param target_name:
    :param refdata:
    :param outdir:
    :param callers: list of callers to use (combination of vardict, freebayes)
    # :return: dict of generated files
    :param vep: boolean whether to run vep on generated vcfs or not
    :param min_alt_frac:
    :return:
    """
    tlib = compose_lib_capture_str(tumor_capture)
    nlib = compose_lib_capture_str(normal_capture)
    normal_sample_str = compose_sample_str(normal_capture)
    tumor_sample_str = compose_sample_str(tumor_capture)

    d = {}
    if 'freebayes' in callers:
        freebayes = Freebayes()
        freebayes.input_bams = [tbam, nbam]
        freebayes.tumorid = tlib
        freebayes.normalid = nlib
        freebayes.somatic_only = True
        freebayes.reference_sequence = refdata['reference_genome']
        freebayes.target_bed = refdata['targets'][target_name]['targets-bed-slopped20']
        freebayes.threads = pipeline.maxcores
        freebayes.min_alt_frac = min_alt_frac
        freebayes.scratch = pipeline.scratch
        freebayes.jobname = "freebayes-somatic/{}".format(tlib)
        freebayes.output = "{}/variants/{}-{}.freebayes-somatic.vcf.gz".format(outdir, tlib, nlib)
        pipeline.add(freebayes)
        d['freebayes'] = freebayes.output

        if vep:
            vep_freebayes = VEP()
            vep_freebayes.input_vcf = freebayes.output
            vep_freebayes.threads = pipeline.maxcores
            vep_freebayes.reference_sequence = refdata['reference_genome']
            vep_freebayes.vep_dir = refdata['vep_dir']
            vep_freebayes.output_vcf = "{}/variants/{}-{}.freebayes-somatic.vep.vcf.gz".format(outdir, tlib, nlib)
            vep_freebayes.jobname = "vep-freebayes-somatic/{}".format(tlib)
            pipeline.add(vep_freebayes)
            d['freebayes'] = vep_freebayes.output_vcf

    if 'vardict' in callers:
        vardict = VarDict(input_tumor=tbam, input_normal=nbam, tumorid=tumor_sample_str,
                          normalid=normal_sample_str,
                          reference_sequence=refdata['reference_genome'],
                          reference_dict=refdata['reference_dict'],
                          target_bed=refdata['targets'][target_name]['targets-bed-slopped20'],
                          output="{}/variants/{}-{}.vardict-somatic.vcf.gz".format(outdir, tlib, nlib),
                          min_alt_frac=min_alt_frac
                          )

        vardict.jobname = "vardict/{}".format(tlib)
        pipeline.add(vardict)

        vep_vardict = VEP()
        vep_vardict.input_vcf = vardict.output
        vep_vardict.threads = pipeline.maxcores
        vep_vardict.reference_sequence = refdata['reference_genome']
        vep_vardict.vep_dir = refdata['vep_dir']
        vep_vardict.output_vcf = "{}/variants/{}-{}.vardict-somatic.vep.vcf.gz".format(outdir, tlib, nlib)
        vep_vardict.jobname = "vep-vardict/{}".format(tlib)
        pipeline.add(vep_vardict)
        d['vardict'] = vep_vardict.output_vcf

    return d
