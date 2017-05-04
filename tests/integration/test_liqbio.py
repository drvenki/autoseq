import json
import os
import subprocess
import unittest

from genomicassertions.readassertions import ReadAssertions
from genomicassertions.variantassertions import VariantAssertions

from autoseq.tests import liqbio_test_outdir
from autoseq.util.path import normpath


class TestLiqbio(unittest.TestCase, VariantAssertions, ReadAssertions):
    returncode = None
    tmpdir = None
    outdir = None
    somatic_vcf = None
    outdir = liqbio_test_outdir

    @classmethod
    def setUpClass(cls):
        cls.jobdb = os.path.join(cls.outdir, "jobdb.json")
        subprocess.check_call("autoseq " +
                              " --ref /tmp/test-genome/autoseq-genome.json " +
                              " --outdir {} ".format(cls.outdir) +
                              " --libdir /tmp/libraries/ " +
                              " --scratch /scratch/tmp/autoseq-integration-tests/liqbio --jobdb {} --cores 2 ".format(cls.jobdb) +
                              " liqbio "
                              " tests/liqbio-test-sample.json", shell=True)

    def test_jobdb(self):
        jobdb = json.load(open(self.jobdb, 'r'))
        self.assertEqual(set([job['status'] for job in jobdb['jobs']]),
                         set(['COMPLETED']))

    def test_vardict_somatic(self):
        vcf = os.path.join(self.outdir, "variants",
                           "NA12877-T-03098849-NA12877-N-03098121-TD1-TT1.vardict-somatic.vcf.gz")

        # TP53 insertion: MU2185182, chr17:g.7578475->G
        # TP53 deletion: MU25947, chr17:g.7577558G>-
        # TP53 DNV: MU52971976, chr17:g.7574003GG>AA
        # PIK3CA hotspot E545K, MU5219, chr3:g.178936091G>A
        # PTEN hotspot R130Q, MU29098, chr10:g.89692905G>A
        # PTEN hotspot R233*, MU589331, chr10:g.89717672C>T
        # AR intron variant, MU50988553, chrX:g.66788924G>A

        self.assertVcfHasSample(vcf, 'NA12877-N-03098121-TD1-TT1')  # N lib id is set
        self.assertVcfHasSample(vcf, 'NA12877-T-03098849')  # T lib id is the merged library

        # deletion is called
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 17, 7577557, 'AG', 'A')

        # insertion is called
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 17, 7578474, 'C', 'CG')

        # dnv (GG>AA) is called as a single event
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 17, 7574003, 'GG', 'AA')

        # PTEN hotspots are called
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 10, 89692905, 'G', 'A')
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 10, 89717672, 'C', 'T')

        # PIK3CA hotspot is called
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 3, 178936091, 'G', 'A')

    def test_panel_bam_readgroups(self):
        bam = os.path.join(self.outdir, "bams", "panel", "NA12877-T-03098849-TD-TT-nodups.bam")
        self.assertBamHeaderElementEquals(bam, 'RG', [{'ID': 'NA12877-T-03098849-TD1-TT1',
                                                       'SM': 'NA12877-T-03098849',
                                                       'LB': 'NA12877-T-03098849-TD1',
                                                       'PL': 'ILLUMINA'}])

        bam = os.path.join(self.outdir, "bams", "panel", "NA12877-CFDNA-03098850-TD-TT-nodups.bam")
        self.assertBamHeaderElementEquals(bam, 'RG', [{'LB': 'NA12877-CFDNA-03098850-TD1',
                                                       'ID': 'NA12877-CFDNA-03098850-TD1-TT1',
                                                       'SM': 'NA12877-CFDNA-03098850',
                                                       'PL': 'ILLUMINA'},
                                                      {'LB': 'NA12877-CFDNA-03098850-TD1',
                                                       'ID': 'NA12877-CFDNA-03098850-TD1-TT2',
                                                       'SM': 'NA12877-CFDNA-03098850',
                                                       'PL': 'ILLUMINA'}]
                                          )

    def test_qdnaseq(self):
        qdnaseq_file_names = ["NA12877-T-03098849-TD1-WGS-qdnaseq.segments.txt",
                              "NA12877-CFDNA-03098850-TD1-WGS-qdnaseq.segments.txt",
                              "NA12877-N-03098121-TD1-WGS-qdnaseq.segments.txt"]
        for qdnaseqf in qdnaseq_file_names:
            absf = os.path.join(self.outdir, "cnv", qdnaseqf)

            with open(absf) as segments:
                ln = segments.readline().strip()
                header = ln.split("\t")
                correct_header = ["chromosome", "start", "end", "bases", "gc", "mappability", "blacklist",
                                  "residual", "use", "readcount", "copynumber", "segmented"]
                self.assertListEqual(header, correct_header)

    def test_wgs_bam_coverage(self):
        bam = os.path.join(self.outdir, "bams/wgs/",
                           "NA12877-T-03098849-TD1-WGS.bam")
        self.assertBamHasCoverageAt(bam, 1, '3', 3617655)  # 1x coverage in that position
        self.assertBamHasCoverageAt(bam, 0, '3', 3618655)  # no coverage in that position

    def test_germline_vcf(self):
        vcf = os.path.join(self.outdir, "variants",
                           "NA12877-N-03098121-TD1-TT1.freebayes-germline.vcf.gz")

        self.assertVcfHasSample(vcf, 'NA12877-N-03098121')  # sample is the name of the normal lib id
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, '3', 178925677, 'G', 'A')  # SNP
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, '17', 7579643, 'CCCCCAGCCCTCCAGGT', 'C')  # deletion

    def test_germline_vcf_with_added_tumor_afs(self):
        vcf = os.path.join(self.outdir, "variants",
                           "NA12877-N-03098121-TD1-TT1-and-NA12877-T-03098849.germline-variants-with-somatic-afs.vcf.gz")

        self.assertVcfHasSample(vcf, 'NA12877-N-03098121')
        self.assertVcfHasSample(vcf, 'NA12877-T-03098849')

        vcf = os.path.join(self.outdir, "variants",
                           "NA12877-N-03098121-TD1-TT1-and-NA12877-CFDNA-03098850.germline-variants-with-somatic-afs.vcf.gz")

        self.assertVcfHasSample(vcf, 'NA12877-N-03098121')
        self.assertVcfHasSample(vcf, 'NA12877-CFDNA-03098850')
