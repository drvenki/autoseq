import os
import unittest

from autoseq.util.library import find_fastqs, get_libdict, Library


class TestLibrary(unittest.TestCase):
    library = 'NA12877-N-03098121-TD1-TT1'
    libdir = 'tests/libraries'

    def test_find_fastq_gz(self):
        """
        test that files on the format *_1.fastq.gz / *_2.fastq.gz are found
        """
        files = find_fastqs(library=self.library, libdir=self.libdir)
        files_basenames = [os.path.basename(f) for f in files[0]] + [os.path.basename(f) for f in files[1]]
        self.assertIn('foo_1.fastq.gz', files_basenames)
        self.assertIn('foo_2.fastq.gz', files_basenames)

    def test_find_fq_gz(self):
        """
        test that files on the format *_1.fq.gz / *_2.fq.gz are found
        """
        files = find_fastqs(library=self.library, libdir=self.libdir)
        files_basenames = [os.path.basename(f) for f in files[0]] + [os.path.basename(f) for f in files[1]]
        self.assertIn('bar_1.fq.gz', files_basenames)
        self.assertIn('bar_2.fq.gz', files_basenames)

    def test_find_RN_DDD(self):
        """
        test that files on the format *R1_nnn.fastq.gz/*R2_nnn.fastq.gz are found
        """
        files = find_fastqs(library=self.library, libdir=self.libdir)
        files_basenames = [os.path.basename(f) for f in files[0]] + [os.path.basename(f) for f in files[1]]
        self.assertIn('baz_R1_001.fastq.gz', files_basenames)
        self.assertIn('baz_R2_001.fastq.gz', files_basenames)
        self.assertIn('baz_R1_999.fastq.gz', files_basenames)
        self.assertIn('baz_R2_999.fastq.gz', files_basenames)

    def test_fastq_path(self):
        """
        test that the components of the path to the files are in place
        """
        files = find_fastqs(library='NA12877-N-03098121-TD1-TT1', libdir='tests/libraries')
        self.assertIn(self.libdir, files[0][0])
        self.assertIn(self.library, files[0][0])

    def test_get_libdict_no_project(self):
        libdict = get_libdict(self.library)

        self.assertEqual(libdict['capture_kit_name'], 'test-regions')
        self.assertEqual(libdict['prep_kit_name'], 'THRUPLEX_DNASEQ')
        self.assertEqual(libdict['sdid'], "NA12877")

    def test_get_libdict_with_project(self):
        library_with_project = 'AL-NA12878-T-sample1-TD1-TT1'
        library_with_wrong_project = 'RX-NA12878-T-sample1-TD1-TT1'
        libdict = get_libdict(library_with_project)
        self.assertEqual(libdict['project_name'], 'ALASCCA')

        with self.assertRaises(ValueError):
            get_libdict(library_with_wrong_project)

    def test_library_formats(self):
        libs = ['NA12877-N-03098121-TD1-TT1',
                'AL-NA12877-N-03098121-TD1-TT1',
                'P-NA12877-N-03098121-TD1-TT1',
                'AL-P-NA12877-N-03098121-TD1-TT1']

        for libnm in libs:
            lib = get_libdict(libnm)
            self.assertIn(lib['sdid'], ['P-NA12877', 'NA12877'])
            self.assertEqual(lib['type'], 'N')
            self.assertEqual(lib['capture_kit_name'], 'test-regions')
