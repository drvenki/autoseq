import json
import logging
import os
import signal

import click
from pypedream import runners

from .alascca import alascca as alascca_cmd
from .liqbio import liqbio as liqbio_cmd
from .liqbio import liqbio_prepare as liqbio_prepare_cmd
from .hopkins import hopkins as hopkins_cmd

__author__ = 'dankle'


@click.group()
@click.option('--ref', default='/nfs/ALASCCA/autoseq-genome/autoseq-genome.json',
              help='json with reference files to use',
              type=str)
@click.option('--outdir', default='/tmp/autoseq-test', help='output directory', type=click.Path())
@click.option('--libdir', default="/tmp", help="directory to search for libraries")
@click.option('--runner_name', default='shellrunner', help='Runner to use.')
@click.option('--loglevel', default='INFO', help='level of logging')
@click.option('--jobdb', default=None, help="sqlite3 database to write job info and stats")
@click.option('--dot_file', default=None, help="write graph to dot file with this name")
@click.option('--cores', default=1, help="max number of cores to allow jobs to use")
@click.option('--scratch', default="/tmp", help="scratch dir to use")
@click.pass_context
def cli(ctx, ref, outdir, libdir, runner_name, loglevel, jobdb, dot_file, cores, scratch):
    setup_logging(loglevel)
    logging.debug("Reading reference data from {}".format(ref))
    ctx.obj = {}
    ctx.obj['refdata'] = load_ref(ref)
    ctx.obj['outdir'] = outdir
    ctx.obj['libdir'] = libdir
    ctx.obj['pipeline'] = None
    ctx.obj['runner'] = get_runner(runner_name, cores)
    ctx.obj['jobdb'] = jobdb
    ctx.obj['dot_file'] = dot_file
    ctx.obj['cores'] = cores
    ctx.obj['scratch'] = scratch

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


def load_ref(ref):
    basepath = os.path.dirname(ref)
    items = ["bwaIndex", "chrsizes", "clinvar", "cosmic", "dbSNP", "exac", "genesGenePred", "genesGtf", "vep_dir",
             "genesGtfGenesOnly", "icgc", "qdnaseq_background", "reference_dict", "reference_genome", "cnvkit-ref",
             "msisites", "targets-bed-slopped20", "targets-interval_list", "targets-interval_list-slopped20"]
    with open(ref, 'r') as fh:
        refjson = json.load(fh)

        def make_paths_absolute(d):
            for k, v in d.items():
                if isinstance(v, dict):
                    make_paths_absolute(v)
                else:
                    if k in items and d[k]:
                        d[k] = os.path.join(basepath, v)
            return d

        refjson_abs = make_paths_absolute(refjson)
        return refjson_abs


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
cli.add_command(hopkins_cmd)
