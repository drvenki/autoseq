import unittest
from autoseq.tools.variantcalling import *


class TestVariantCalling(unittest.TestCase):
    def test_freebayes(self):
        freebayes = Freebayes()
        freebayes.input_bams = ["input.bam"]
        freebayes.reference_sequence = "dummy.fasta"
        freebayes.target_bed = "dummy_targets.bed"
        freebayes.output = "output.txt"
        cmd = freebayes.command()
        self.assertIn('input.bam', cmd)
        self.assertIn('dummy.fasta', cmd)
        self.assertIn('dummy_targets.bed', cmd)
        self.assertIn('output.txt', cmd)

    def test_vardict(self):
        vardict = VarDict()
        vardict.input_tumor = "input_tumor.bam"
        vardict.input_normal = "input_normal.bam"
        vardict.tumorid = "tumor_id"
        vardict.normalid = "normal_id"
        vardict.reference_sequence = "dummy.fasta"
        vardict.reference_dict = {}
        vardict.target_bed = "dummy_targets.bed"
        vardict.output = "output.txt"
        cmd = vardict.command()
        self.assertIn('input_tumor.bam', cmd)
        self.assertIn('input_normal.bam', cmd)
        self.assertIn('tumor_id', cmd)
        self.assertIn('normal_id', cmd)
        self.assertIn('dummy.fasta', cmd)
        self.assertIn('dummy_targets.bed', cmd)
        self.assertIn('output.txt', cmd)

    def test_vep(self):
        vep = VEP()
        vep.input_vcf = "input.vcf"
        vep.output_vcf = "output.vcf"
        vep.reference_sequence = "dummy.fasta"
        vep.vep_dir = "dummy_dir"
        cmd = vep.command()
        self.assertIn('input.vcf', cmd)
        self.assertIn('dummy.fasta', cmd)
        self.assertIn('dummy_dir', cmd)
        self.assertIn('output.vcf', cmd)

    def test_vcf_add_sample(self):
        vcf_add_sample = VcfAddSample()
        vcf_add_sample.input_vcf = "input.vcf"
        vcf_add_sample.input_bam = "input.bam"
        vcf_add_sample.samplename = "dummy_name"
        vcf_add_sample.output = "output.vcf"
        cmd = vcf_add_sample.command()
        self.assertIn('input.vcf', cmd)
        self.assertIn('input.bam', cmd)
        self.assertIn('dummy_name', cmd)
        self.assertIn('output.vcf', cmd)

    # XXX CONTINUE HERE: ADD AN EXTRA TEST TO TEST THE GZ FUNCTIONALITY.
    # ALSO, ADD TESTS FOR BRANCHES IN THE VEP JOB CODE.

    def test_vcf_filter(self):
        vcf_filter = VcfFilter()
        vcf_filter.input = "input.vcf"
        vcf_filter.filter = "test_filter"
        vcf_filter.output = "output.vcf"
        cmd = vcf_filter.command()
        self.assertIn('input.vcf', cmd)
        self.assertIn('test_filter', cmd)
        self.assertIn('output.vcf', cmd)

    def test_curl_split_and_left_align(self):
        curl_split_and_left_align = CurlSplitAndLeftAlign()
        curl_split_and_left_align.remote = "dummy_remote"
        curl_split_and_left_align.input_reference_sequence = "dummy_reference.fasta"
        curl_split_and_left_align.input_reference_sequence_fai = "dummy_reference.fasta.fai"
        curl_split_and_left_align.output = "output.vcf"
        cmd = curl_split_and_left_align.command()
        self.assertIn('dummy_remote', cmd)
        self.assertIn('dummy_reference.fasta', cmd)
        self.assertIn('output.vcf', cmd)

    def test_install_vep(self):
        install_vep = InstallVep()
        install_vep.output_dir = "dummy_output_dir"
        cmd = install_vep.command()
        self.assertIn('dummy_output_dir', cmd)

    # XXX CONTINUE HERE: WRITE TESTS FOR THIS FUNCTION. POSSIBLY ALSO REFACTOR IT TO MAKE THE INTERFACE
    # CLEANER, AND IMPROVE THE DOCUMENTATION TO IT.
    def test_call_somatic_variants(self):
        pass
