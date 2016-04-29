import json
import logging
import os

import click
import time

import sys

from autoseq.pipeline.liqbio import LiqBioPipeline
from autoseq.util.path import mkdir
from autoseq.pipeline.alasccapipeline import AlasccaPipeline


@click.command()
@click.option('--libdir', default=None, help="directory to search for libraries")
@click.argument('sample', type=click.File('r'))
@click.pass_context
def liqbio(ctx, libdir, sample):
    logging.info("Running Liquid Biopsy pipeline")
    logging.info("Sample is {}".format(sample))

    logging.debug("Reading sample config from {}".format(sample))
    sampledata = json.load(sample)

    if ctx.obj['jobdb']:
        mkdir(os.path.dirname(ctx.obj['jobdb']))

    ctx.obj['pipeline'] = LiqBioPipeline(sampledata=sampledata, refdata=ctx.obj['refdata'],
                                         outdir=ctx.obj['outdir'], libdir=libdir, maxcores=ctx.obj['cores'],
                                         debug=ctx.obj['debug'], runner=ctx.obj['runner'],
                                         jobdb=ctx.obj['jobdb'], dot_file=ctx.obj['dot_file'])

    # start main analysis
    ctx.obj['pipeline'].start()
    #
    logging.info("Waiting for pipeline to finish.")
    while ctx.obj['pipeline'].is_alive():
        logging.debug("Waiting for LiqBioPipeline")
        time.sleep(5)

    # # return_code from run_pipeline() will be != 0 if the pipeline fails
    sys.exit(ctx.obj['pipeline'].exitcode)
