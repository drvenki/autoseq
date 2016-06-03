import logging
import os
import re

from autoseq.util.path import normpath

logger = logging.getLogger(__name__)
def get_prep_kit_name_from_id(prep_id):
    prep_kit_lookup = {"BN": "BIOO_NEXTFLEX",
                       "KH": "KAPA_HYPERPREP",
                       "TD": "THRUPLEX_DNASEQ",
                       "TP": "THRUPLEX_PLASMASEQ",
                       "TF": "THRUPLEX_FD",
                       "TS": "TRUSEQ_RNA",
                       "NN": "NEBNEXT_RNA",
                       "VI": "VILO_RNA"}

    shortname = prep_id[0:2]
    return prep_kit_lookup[shortname]


def get_capture_kit_name_from_id(capture_id):
    capture_kit_loopkup = {"CS": "clinseq_v3_targets",
                           "CZ": "clinseq_v4",
                           "EX": "EXOMEV3",
                           "EO": "EXOMEV1",
                           "RF": "fusion_v1",
                           "CC": "core_design",
                           "CD": "discovery_coho",
                           "CB": "big_design",
                           "AL": "alascca_targets",
                           "TT": "test-regions"
                           }
    if capture_id == 'WGS':
        return 'lowpass_wgs'
    shortname = capture_id[0:2]
    return capture_kit_loopkup[shortname]


def get_libdict(library_id):
    """Create a dictionary from a library_id with the individual fields.

    :rtype: dict[str,str]
    """
    elems = library_id.split("-")
    elem_nms = ['sdid', 'type', 'sample_id', 'prep_id', 'capture_id']
    if len(elems) == 6 and elems[0] == 'P':
        elems = elems[1:6]
        elems[0] = "P-{}".format(elems[0])
    d = dict(zip(elem_nms, elems))
    d['prep_kit_name'] = get_prep_kit_name_from_id(d['prep_id'])
    d['capture_kit_name'] = get_capture_kit_name_from_id(d['capture_id'])
    d['library_id'] = library_id
    return d


def find_fastqs(library, libdir):
    """Find fastq files for a given library id in a given direcory.

        Returns a tuple with two lists:
    (['foo_1.fastq.gz', 'bar_1.fastq.gz'], # read 1
     ['foo_2.fastq.gz', 'bar_2.fastq.gz'])

    Supports the following file naming convenstions:
    *_1.fastq.gz / *_2.fastq.gz
    *_1.fq.gz / *_2.fq.gz
    *R1_nnn.fastq.gz / *R2_nnn.fastq.gz

    :rtype: tuple[str,str]
    """

    regex_fq1 = '(.+)(_1\.fastq.gz|_1\.fq.gz|R1_\d{3}.fastq.gz)'
    regex_fq2 = '(.+)(_2\.fastq.gz|_2\.fq.gz|R2_\d{3}.fastq.gz)'

    d = normpath(os.path.join(libdir, library))
    logger.debug("Looking for fastq files for library {library} in {libdir}".format(library=library, libdir=libdir))

    fq1s = []
    fq2s = []

    for f in os.listdir(d):
        match1 = re.search(regex_fq1, f)
        if match1:
            fq1s.append("".join(match1.groups()))
        match2 = re.search(regex_fq2, f)
        if match2:
            fq2s.append("".join(match2.groups()))

    fq1s.sort()
    fq2s.sort()

    logging.debug("Found {}".format((fq1s, fq2s)))
    return fq1s, fq2s
