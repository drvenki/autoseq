import json
import os
import unittest

from autoseq.util.library import find_fastqs

from autoseq.util.path import normpath


class TestLibrary(unittest.TestCase):
    def test_find_fastq_gz(self):
        files = find_fastqs(library='NA12877-N-03098121-TD1-TT1', libdir='tests/libraries')
        files_basenames = [os.path.basename(f) for f in files[0]] + [os.path.basename(f) for f in files[1]]
        self.assertIn('foo_1.fastq.gz', files_basenames)
        self.assertIn('foo_2.fastq.gz', files_basenames)

    def test_find_fq_gz(self):
        files = find_fastqs(library='NA12877-N-03098121-TD1-TT1', libdir='tests/libraries')
        files_basenames = [os.path.basename(f) for f in files[0]] + [os.path.basename(f) for f in files[1]]
        self.assertIn('bar_1.fq.gz', files_basenames)
        self.assertIn('bar_2.fq.gz', files_basenames)

    def test_find_RN_DDD(self):
        files = find_fastqs(library='NA12877-N-03098121-TD1-TT1', libdir='tests/libraries')
        files_basenames = [os.path.basename(f) for f in files[0]] + [os.path.basename(f) for f in files[1]]
        self.assertIn('baz_R1_001.fastq.gz', files_basenames)
        self.assertIn('baz_R2_001.fastq.gz', files_basenames)
        self.assertIn('baz_R1_999.fastq.gz', files_basenames)
        self.assertIn('baz_R2_999.fastq.gz', files_basenames)
