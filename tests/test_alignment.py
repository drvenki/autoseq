import unittest

from autoseq.tools.alignment import Bwa


class TestAlignment(unittest.TestCase):

    def test_bwa_pe(self):
        """test bwa paired end"""
        bwa = Bwa()
        bwa.input_fastq1 = "foo_1.fq"
        bwa.input_fastq2 = "foo_2.fq"
        bwa.input_reference_sequence = "ref.fasta"
        bwa.readgroup = "__readgroup__"
        bwa.output = "out.bam"
        cmd = bwa.command()
        self.assertIn('foo_1.fq', cmd)
        self.assertIn('foo_2.fq', cmd)

    def test_bwa_se(self):
        """ test bwa single end"""
        bwa = Bwa()
        bwa.input_fastq1 = "foo_1.fq"
        bwa.input_fastq2 = None
        bwa.input_reference_sequence = "ref.fasta"
        bwa.readgroup = "__readgroup__"
        bwa.output = "out.bam"
        cmd = bwa.command()
        self.assertIn('foo_1.fq', cmd)
        self.assertNotIn('foo_2.fq', cmd)

    def test_bwa_removes_temporary_logs(self):
        bwa = Bwa()
        bwa.input_fastq1 = "foo_1.fq"
        bwa.input_fastq2 = None
        bwa.input_reference_sequence = "ref.fasta"
        bwa.readgroup = "__readgroup__"
        bwa.output = "out.bam"
        cmd = bwa.command()

        self.assertTrue(cmd.endswith('rm out.bam.bwa.log out.bam.samblaster.log'),
                        msg="rm temp logs must be the final part of the command")
