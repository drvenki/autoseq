import json
import logging
import signal

import click
from pypedream import runners

from .alascca import alascca as alascca_cmd
from .liqbio import liqbio as liqbio_cmd
from .liqbio import liqbio_prepare as liqbio_prepare_cmd

__author__ = 'dankle'


@click.group()
@click.option('--ref', default='/nfs/ALASCCA/autoseq-genome/autoseq-genome.json',
              help='json with reference files to use',
              type=click.File('r'))
@click.option('--outdir', default='/tmp/pyautoseq-test', help='output directory', type=click.Path())
@click.option('--runner_name', default='shellrunner', help='Runner to use.')
@click.option('--loglevel', default='INFO', help='level of logging')
@click.option('--jobdb', default=None, help="sqlite3 database to write job info and stats")
@click.option('--dot_file', default=None, help="write graph to dot file with this name")
@click.option('--cores', default=1, help="write graph to dot file with this name")
@click.option('--debug', default=False, is_flag=True)
@click.option('--tmpdir', default=None, help="node-local tmpdir to use")
@click.pass_context
def cli(ctx, ref, outdir, runner_name, loglevel, jobdb, dot_file, cores, tmpdir, debug):
    setup_logging(loglevel)
    logging.debug("Reading reference data from {}".format(ref))
    ctx.obj = {}
    ctx.obj['refdata'] = json.load(ref)
    ctx.obj['outdir'] = outdir
    ctx.obj['pipeline'] = None
    ctx.obj['runner'] = get_runner(runner_name, cores)
    ctx.obj['jobdb'] = jobdb
    ctx.obj['dot_file'] = dot_file
    ctx.obj['cores'] = cores
    ctx.obj['tmpdir'] = tmpdir
    ctx.obj['debug'] = debug

    def capture_sigint(sig, frame):
        """
        Capture ctrl-c (or SIGINT sent in other ways).
        1. update remote log
        :param sig:
        :param frame:
        :return:
        """
        try:
            ctx.obj['pipeline'].stop()
            logging.info("Stopping jobs...")
        except AttributeError:
            logging.debug("No pipeline to stop.")

    signal.signal(signal.SIGINT, capture_sigint)
    signal.signal(signal.SIGTERM, capture_sigint)


def get_runner(runner_name, maxcores):
    try:
        module = __import__("pypedream.runners." + runner_name, fromlist="runners")
        runner_class = getattr(module, runner_name.title())

        if runner_name == 'localqrunner':
            runner = runner_class(maxcores)
        else:
            runner = runner_class()
        return runner
    except ImportError:
        print "Couldn't find runner " + runner_name + ". Available Runners:"
        import inspect
        for name, obj in inspect.getmembers(runners):
            if name != "runner" and "runner" in name:
                print "- " + name
        raise ImportError


def setup_logging(loglevel="INFO"):
    """
    Set up logging
    :param loglevel: loglevel to use, one of ERROR, WARNING, DEBUG, INFO (default INFO)
    :return:
    """
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level,
                        format='%(levelname)s %(asctime)s %(funcName)s - %(message)s')
    logging.info("Started log with loglevel %(loglevel)s" % {"loglevel": loglevel})

cli.add_command(alascca_cmd)
cli.add_command(liqbio_cmd)
cli.add_command(liqbio_prepare_cmd)
