import json
import logging
import os
import signal
import subprocess
import sys
import time

import click
from pypedream import runners

from autoseq.pipeline.alasccapipeline import AlasccaPipeline

__author__ = 'dankle'


@click.command()
@click.option('--ref', default='/proj/b2010040/private/nobackup/autoseq-genome/autoseq-genome.json',
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
@click.argument('sample', type=click.File('r'))
def main(sample, ref, outdir, runner_name, loglevel, jobdb, dot_file, cores, tmpdir, debug):
    setup_logging(loglevel)
    logging.debug("Reading sample config from {}".format(sample))
    sampledata = json.load(sample)
    logging.debug("Reading reference data from {}".format(ref))
    refdata = json.load(ref)

    # if tmpdir is set, transfer all raw data and all previously processed data for this sample there
    # and run the pipeline there. Transfer back when completed.
    sampledata_to_use = sampledata
    outdir_to_use = outdir
    if tmpdir:
        # sampledata_to_use = transfer_data_to_local(sampledata, tmpdir, outdir)
        new_outdir = "{}/{}".format(tmpdir, sampledata['REPORTID'])

        rsync_dirs(new_outdir, outdir)

    # autoseq_analysis.sampledata = new_sampledata
    # autoseq_analysis.setup_pipeline(outdir=new_outdir)

    if jobdb:
        mkdir(os.path.dirname(jobdb))

    # sampledata, refdata, outdir, maxcores=1, scriptdir="/tmp/.pypedream", debug=False
    runner = get_runner(runner_name, cores)
    p = AlasccaPipeline(sampledata=sampledata_to_use, refdata=refdata, outdir=outdir_to_use,
                        maxcores=cores, debug=debug, runner=runner, jobdb=jobdb, dot_file=dot_file)

    def capture_sigint(sig, frame):
        """
        Capture ctrl-c (or SIGINT sent in other ways).
        1. update remote log
        :param sig:
        :param frame:
        :return:
        """
        logging.info("Stopping jobs...")
        p.stop()

    signal.signal(signal.SIGINT, capture_sigint)
    signal.signal(signal.SIGTERM, capture_sigint)

    logging.debug("outdir = {}".format(p.outdir))

    logging.debug("Writing dot file: {}".format(dot_file))

    # start main analysis
    logging.info("Starting AlasccaPipeline.")
    p.start()

    logging.info("Waiting for AlasccaPipeline to finish.")
    while p.is_alive():
        logging.debug("Waiting for AutoseqPipeline")
        time.sleep(5)

    # return_code from run_pipeline() will be != 0 if the pipeline fails
    sys.exit(p.exitcode)


def transfer_data_to_local(sampledata, tmpdir, final_outdir):
    rawdata_dir = "{}/raw/".format(tmpdir)
    new_outdir = "{}/{}".format(tmpdir, sampledata['REPORTID'])
    mkdir(rawdata_dir)
    mkdir(new_outdir)

    logging.debug("Syncing final output dir")
    rsync_dirs(final_outdir, new_outdir)

    logging.debug("Fetching raw data")
    new_sampledata = fetch_raw_data(sampledata, rawdata_dir)
    return new_sampledata


def mkdir(dir):
    """ Create a directory if it doesn't exist
    :param dir: dir to create
    """
    if not os.path.isdir(dir):
        try:
            os.makedirs(dir)
        except OSError:
            logging.error("Couldn't create directory {}".format(dir))
    else:  # if dir already exists, do nothing
        pass


def rsync_dirs(src, target):
    """
    Sync two a source dir to target
    :param src:
    :param target:
    :return:
    """
    src = os.path.expandvars(os.path.expanduser(src))
    target = os.path.expandvars(os.path.expanduser(target))
    rsync_results_cmd = ['rsync', '--stats', '--size-only', '--delete', '-avP', src + "/", target + "/"]
    logging.debug("Running rsync with command")
    logging.debug("{}".format(" ".join(rsync_results_cmd)))
    subprocess.check_call(rsync_results_cmd, stderr=open('/dev/null', 'w'), stdout=open('/dev/null', 'w'))


def fetch_raw_data(sampledata, rawdata_dir):
    """
    Copy data for a single report to a directory
    """
    items_to_copy = ["PANEL_TUMOR_FQ1", "PANEL_TUMOR_FQ2", "PANEL_NORMAL_FQ1", "PANEL_NORMAL_FQ2",
                     "WGS_TUMOR_FQ1", "WGS_TUMOR_FQ2", "WGS_NORMAL_FQ1", "WGS_NORMAL_FQ2",
                     "RNASEQ_FQ1", "RNASEQ_FQ2", "RNASEQCAP_FQ1", "RNASEQCAP_FQ2"]

    for item in items_to_copy:
        if sampledata[item] is not None and sampledata[item] is not "NA":
            new_fqs = []
            for f in sampledata[item]:
                local_file = os.path.abspath("{}/{}".format(rawdata_dir, os.path.expanduser(f)))
                new_fqs.append(local_file)
                rsync_file(f, local_file)
            sampledata[item] = new_fqs
    return sampledata


def rsync_file(src, target):
    """
    Sync src file to target
    :param src:
    :param target:
    """
    src = os.path.expandvars(os.path.expanduser(src))
    target = os.path.expandvars(os.path.expanduser(target))
    logging.info("Copying {} to local work dir".format(src))
    rsync_command = ["rsync", "--stats", "--size-only", src, target]
    mkdir(os.path.dirname(target))
    logging.debug("Executing {}".format(rsync_command))
    subprocess.check_call(rsync_command, stderr=open('/dev/null', 'w'), stdout=open('/dev/null', 'w'))


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


if __name__ == "__main__":
    sys.exit(main())
