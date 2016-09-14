import unittest

from mock import patch, MagicMock

from autoseq.tools.alignment import Bwa, Skewer, align_library


class TestAlignment(unittest.TestCase):
    def test_bwa_pe(self):
        """
        test bwa paired end, should contain string "foo_1.fq  foo_2.fq" with a double-space between the file names
        """
        bwa = Bwa()
        bwa.input_fastq1 = "foo_1.fq"
        bwa.input_fastq2 = "foo_2.fq"
        bwa.input_reference_sequence = "ref.fasta"
        bwa.readgroup = "__readgroup__"
        bwa.output = "out.bam"
        cmd = bwa.command()
        self.assertIn('foo_1.fq  foo_2.fq', cmd)

    def test_bwa_se(self):
        """
        test bwa single end, should not include any reference to foo_2.fq
        """
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
        """
        test that the command includes a removal of the temporary log files
        """
        bwa = Bwa()
        bwa.input_fastq1 = "foo_1.fq"
        bwa.input_fastq2 = None
        bwa.input_reference_sequence = "ref.fasta"
        bwa.readgroup = "__readgroup__"
        bwa.output = "out.bam"
        cmd = bwa.command()
        self.assertTrue(cmd.endswith('rm out.bam.bwa.log out.bam.samblaster.log'),
                        msg="rm temp logs must be the final part of the command")

    @patch('uuid.uuid4')
    def test_skewer_deletes_tmp(self, mock_uuid):
        """
        test that skewer includes commands to remove the temp dir it uses
        """
        mock_uuid.return_value = "foo"

        skewer = Skewer()
        skewer.input1 = "in_1.fq.gz"
        skewer.input2 = "in_2.fq.gz"
        skewer.output1 = "out_1.fq.gz"
        skewer.output2 = "out_2.fq.gz"
        skewer.stats = "stats.txt"
        skewer.scratch = "/scratch"
        cmd = skewer.command()

        self.assertTrue(cmd.endswith('rm -r /scratch/skewer-foo'))

    @patch('uuid.uuid4')
    def test_skewer_pe_uses_both_fqs(self, mock_uuid):
        """
        test that both input1 and input2 are used in the command line
        and that both output1 and output2 are copied back
        """
        mock_uuid.return_value = "foo"

        skewer = Skewer()
        skewer.input1 = "in_1.fq.gz"
        skewer.input2 = "in_2.fq.gz"
        skewer.output1 = "/path/to/out_1.fq.gz"
        skewer.output2 = "/path/to/out_2.fq.gz"
        skewer.stats = "stats.txt"
        skewer.scratch = "/scratch"
        cmd = skewer.command()
        self.assertIn("in_1.fq.gz  in_2.fq.gz", cmd)
        self.assertIn("cp /scratch/skewer-foo/skewer-trimmed-pair1.fastq.gz /path/to/out_1.fq.gz", cmd)
        self.assertIn("cp /scratch/skewer-foo/skewer-trimmed-pair2.fastq.gz /path/to/out_2.fq.gz", cmd)

    @patch('uuid.uuid4')
    def test_skewer_se_only_uses_input1(self, mock_uuid):
        """
        test that skewer doesn't try to copy back fastq2 in single-end mode
        """
        mock_uuid.return_value = "foo"
        skewer = Skewer()
        skewer.input1 = "in_1.fq.gz"
        skewer.input2 = None
        skewer.output1 = "/path/to/out_1.fq.gz"
        skewer.output2 = "/path/to/out_2.fq.gz"
        skewer.stats = "stats.txt"
        skewer.scratch = "/scratch"
        cmd = skewer.command()

        self.assertNotIn("/path/to/out_2.fq.gz", cmd)

    def test_skewer_raises_valuerror_if_output_not_gzipped(self):
        """
        test that skewer raises a ValueError if the user tries to use uncompressed output files
        """
        skewer = Skewer()
        skewer.input1 = "in_1.fq.gz"
        skewer.input2 = None
        skewer.output1 = "/path/to/out_1.fq"
        skewer.output2 = "/path/to/out_2.fq.gz"
        skewer.stats = "stats.txt"
        skewer.scratch = "/scratch"
        with self.assertRaises(ValueError):
            skewer.command()

        skewer.output1 = "/path/to/out_1.fq.gz"
        skewer.output2 = "/path/to/out_2.fq"
        with self.assertRaises(ValueError):
            skewer.command()
