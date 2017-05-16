import json
import logging
import os
import sys
import time

import click

from autoseq.pipeline.liqbio import LiqBioPipeline
from autoseq.util.clinseq_barcode import extract_clinseq_barcodes
from autoseq.util.path import mkdir


@click.command()
@click.argument('sample', type=click.File('r'))
@click.pass_context
def liqbio(ctx, sample):
    logging.info("Running Liquid Biopsy pipeline")
    logging.info("Sample is {}".format(sample))

    logging.debug("Reading sample config from {}".format(sample))
    sampledata = json.load(sample)

    if ctx.obj['jobdb']:
        mkdir(os.path.dirname(ctx.obj['jobdb']))

    ctx.obj['pipeline'] = LiqBioPipeline(sampledata=sampledata,
                                         refdata=ctx.obj['refdata'],
                                         outdir=ctx.obj['outdir'],
                                         libdir=ctx.obj['libdir'],
                                         maxcores=ctx.obj['cores'],
                                         runner=ctx.obj['runner'],
                                         jobdb=ctx.obj['jobdb'],
                                         dot_file=ctx.obj['dot_file'],
                                         scratch=ctx.obj['scratch'])

    # start main analysis
    ctx.obj['pipeline'].start()
    #
    logging.info("Waiting for pipeline to finish.")
    while ctx.obj['pipeline'].is_alive():
        logging.debug("Waiting for LiqBioPipeline")
        time.sleep(5)

    # # return_code from run_pipeline() will be != 0 if the pipeline fails
    sys.exit(ctx.obj['pipeline'].exitcode)


@click.command()
@click.option('--outdir', required=True, help="directory to write config files")
@click.argument('barcodes-filename', type=str)
@click.pass_context
def liqbio_prepare(ctx, outdir, barcodes_filename):
    logging.info("Extracting clinseq barcodes from input file: " + barcodes_filename)
    clinseq_barcodes = extract_clinseq_barcodes(barcodes_filename)

    logging.info("Generating sample dictionaries from the input clinseq barcodes.")

    sample_dicts = make_sample_dicts(libraries=clinseq_barcodes)
    for sdid in sample_dicts:
        fn = "{}/{}.json".format(outdir, sdid)
        with open(fn, 'w') as f:
            json.dump(sample_dicts[sdid], f, sort_keys=True, indent=4)


# FIXME: Change this function to implement the new sample dictionary structure. Also, break into modular pieces.
def make_sample_dicts(libraries):
    sdids = set([lib['sdid'] for lib in libraries])
    dicts = {}
    for sdid in sdids:
        logging.debug("Creating config file for SDID {}".format(sdid))
        panel_t_lib = [lib['library_id'] for lib in libraries if
                       lib['sdid'] == sdid and lib['type'] == "T" and lib['capture_id'] != "WGS"]
        panel_n_lib = [lib['library_id'] for lib in libraries if
                       lib['sdid'] == sdid and lib['type'] == "N" and lib['capture_id'] != "WGS"]
        panel_p_libs = [lib['library_id'] for lib in libraries if
                        lib['sdid'] == sdid and lib['type'] == "CFDNA" and lib['capture_id'] != "WGS"]

        wgs_t_lib = [lib['library_id'] for lib in libraries if
                     lib['sdid'] == sdid and lib['type'] == "T" and lib['capture_id'] == "WGS"]
        wgs_n_lib = [lib['library_id'] for lib in libraries if
                     lib['sdid'] == sdid and lib['type'] == "N" and lib['capture_id'] == "WGS"]
        wgs_p_libs = [lib['library_id'] for lib in libraries if
                      lib['sdid'] == sdid and lib['type'] == "CFDNA" and lib['capture_id'] == "WGS"]

        def fix_lib(lib):
            """lib is a vector of libraries.
            If it has no contents, return None
            If it has a single element, return it
            If it has more than one element, raise error
            """
            if lib == []:
                return None
            if len(lib) == 1:
                return lib[0]
            if len(lib) > 1:
                raise ValueError("Too many libs for SDID {}. Expected 1, got {}".format(sdid, lib))

        panel_t_lib = fix_lib(panel_t_lib)
        panel_n_lib = fix_lib(panel_n_lib)
        wgs_t_lib = fix_lib(wgs_t_lib)
        wgs_n_lib = fix_lib(wgs_n_lib)

        d = {'sdid': sdid,
             'panel': {
                 'T': panel_t_lib,
                 'N': panel_n_lib,
                 'CFDNA': panel_p_libs
             },
             'wgs': {
                 'T': wgs_t_lib,
                 'N': wgs_n_lib,
                 'CFDNA': wgs_p_libs
             }
        }
        dicts[sdid] = d

    return dicts
