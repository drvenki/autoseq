import json
import logging
import os

import click
import time

import sys

from autoseq.util.path import mkdir
from autoseq.pipeline.alasccapipeline import AlasccaPipeline


@click.command()
@click.argument('sample', type=click.File('r'))
@click.pass_context
def alascca(ctx, sample):
    logging.info("Running Alascca pipeline")
    logging.info("sample is {}".format(sample))

    logging.debug("Reading sample config from {}".format(sample))
    sampledata = json.load(sample)

    if ctx.obj['jobdb']:
        mkdir(os.path.dirname(ctx.obj['jobdb']))

    # sampledata, refdata, outdir, maxcores=1, scriptdir="/tmp/.pypedream", debug=False
    p = AlasccaPipeline(sampledata=ctx.obj['sampledata_to_use'], refdata=ctx.obj['refdata'],
                        outdir=ctx.obj['outdir_to_use'], maxcores=ctx.obj['cores'],
                        debug=ctx.obj['debug'], runner=ctx.obj['runner'],
                        jobdb=ctx.obj['jobdb'], dot_file=ctx.obj['dot_file'])

    logging.debug("outdir = {}".format(p.outdir))

    # start main analysis
    logging.info("Starting AlasccaPipeline.")
    p.start()

    logging.info("Waiting for AlasccaPipeline to finish.")
    while p.is_alive():
        logging.debug("Waiting for AutoseqPipeline")
        time.sleep(5)

    # return_code from run_pipeline() will be != 0 if the pipeline fails
    sys.exit(p.exitcode)
