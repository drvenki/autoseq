#!/usr/bin/env python
import argparse
import json
import logging
import sys

from autoseq.util.path import normpath
from autoseq.util.readpair import Readpair
from autoseq.util.report import Report


def main():
    opts = parse_cli()
    setup_logging(opts.loglevel)
    logging.info("Using template file {}".format(opts.template))
    logging.info("Using prefix {}".format(opts.prefix))

    readpairs = Readpair.fromDir(normpath(opts.fastq_dir))

    reports = Report.fromFile(opts.reports)
    for report in reports:
        jsonfn = normpath(opts.prefix) + report.REPORTID + ".json"
        logging.info("Writing json to {}".format(jsonfn))
        jsonfh = open(jsonfn, 'w')
        json.dump(report.to_dict(readpairs), jsonfh, indent=4, sort_keys=True)


def parse_cli():
    """
    Parse command line argument
    :return: a Namespace object from argparse.parse_args()
    """
    ap = argparse.ArgumentParser()
    ap.add_argument('--reports', help="report to convert", action="store", required=True)
    ap.add_argument('--prefix', help="output json file", action="store", default="~/tmp/report_")
    ap.add_argument('--template', help="template json file", action="store", default="template.json")
    ap.add_argument('--fastq_dir', help="dir containing raw data in lib/file_1.fastq.gz structure",
                    default="/proj/b2010040/INBOX/exomes")
    ap.add_argument("--loglevel", help="level of logging", default='INFO', type=str,
                    choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'])
    return ap.parse_args()


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
