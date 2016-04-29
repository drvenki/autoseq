import json
import logging
import os
import sys
import time

import click

from autoseq.pipeline.hopkins import HopkinsMappingPipeline
from autoseq.util.path import mkdir


@click.command()
@click.argument('samples', type=str)
@click.pass_context
def hopkins(ctx, samples):
    logging.info("Running Hopkins mapping pipeline")

    if ctx.obj['jobdb']:
        mkdir(os.path.dirname(ctx.obj['jobdb']))

    ctx.obj['pipeline'] = HopkinsMappingPipeline(sampledata=samples, refdata=ctx.obj['refdata'],
                                                 outdir=ctx.obj['outdir'], maxcores=ctx.obj['cores'],
                                                 debug=ctx.obj['debug'], runner=ctx.obj['runner'],
                                                 jobdb=ctx.obj['jobdb'], dot_file=ctx.obj['dot_file'])

    # start main analysis
    ctx.obj['pipeline'].start()
    #
    logging.info("Waiting for pipeline to finish.")
    while ctx.obj['pipeline'].is_alive():
        logging.debug("Waiting for HopkinsMappingPipeline")
        time.sleep(5)

    # # return_code from run_pipeline() will be != 0 if the pipeline fails
    sys.exit(ctx.obj['pipeline'].exitcode)
