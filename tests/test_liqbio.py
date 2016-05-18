import tempfile
import unittest

import sys

import time
from genomicassertions.variantassertions import VariantAssertions

from autoseq.cli.cli import load_ref, get_runner
from autoseq.pipeline.liqbio import LiqBioPipeline
from autoseq.util.path import normpath


class TestWorkflow(unittest.TestCase, VariantAssertions):
    returncode = None
    tmpdir = None
    somatic_vcf = None

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()

        ref = load_ref(normpath("~/test-genome/autoseq-genome.json"))
        sampledata = {
            "sdid": "NA12877",
            "panel": {
                "T": "NA12877-T-03098849-TD1-TT1",
                "N": "NA12877-N-03098121-TD1-TT1",
                "P": ["NA12877-P-03098850-TD1-TT1"]
            },
            "wgs": {
                "T": "NA12877-T-03098849-TD1-WGS",
                "N": "NA12877-N-03098121-TD1-WGS",
                "P": ["NA12877-P-03098850-TD1-WGS"]
            }
        }

        libdir = "tests/libraries"

        outdir = cls.tmpdir
        maxcores = 1
        runner = get_runner("shellrunner", maxcores)

        p = LiqBioPipeline(sampledata, ref, outdir, libdir, analysis_id="test", maxcores=maxcores, runner=runner)
        print p.outdir
        #p.start()

        while p.is_alive():
            time.sleep(5)

    def test_test(self):
        self.assertEqual(1, 1)

        # sampledata=sampledata, refdata=ctx.obj['refdata'],
        #                outdir=ctx.obj['outdir'], libdir=libdir, maxcores=ctx.obj['cores'],
        #                debug=ctx.obj['debug'], runner=ctx.obj['runner'],
        #                jobdb=ctx.obj['jobdb'], dot_file=ctx.obj['dot_file'],
        #                scratch=ctx.obj['scratch'])
        #
        #
        # cls.returncode = subprocess.check_call(["./tests/workflow.sh", cls.tmpdir])
        # cls.somatic_vcf = cls.tmpdir + "/panel/virtual-tumor.somatic.vcf.gz"

    # cls.returncode = 0
