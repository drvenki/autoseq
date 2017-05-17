import collections
from autoseq.util.orderform import parse_orderform

# Fields defining a unique library capture item:
UniqueCapture = collections.namedtuple(
    'UniqueCapture',

    ['project',
    'sdid',
    'sample_type',
    'sample_id',
    'library_kit_id',
    'capture_kit_id']
)


# A named tuple representation of a clinseq barcode:
ClinseqBarcode = collections.namedtuple(
    'ClinseqBarcode',

    ['project',
     'sdid',
     'sample_type',
     'sample_id',
     'prep_id',
     'capture_id']
)


def parse_clinseq_barcode(clinseq_barcode):
    """
    Parses and validates the specified clinseq barcode, and produces a corresponding
    ClinseqBarcode namedtuple.

    :param clinseq_barcode: A clinseq barcode string 
    :return: ClinseqBarcode named tuple
    """

    if not clinseq_barcode_is_valid(clinseq_barcode):
        raise ValueError("Invalid clinseq barcode: " + clinseq_barcode)

    return ClinseqBarcode(parse_project(clinseq_barcode),
                          parse_sdid(clinseq_barcode),
                          parse_sample_type(clinseq_barcode),
                          parse_sample_id(clinseq_barcode),
                          parse_prep_id(clinseq_barcode),
                          parse_capture_id(clinseq_barcode))


def extract_unique_capture(clinseq_barcode):
    """
    Parses the specified clinseq barcode and produces a corresponding
    UniqueCapture representing the unique library capture corresponding
    to this clinseq barcode.

    :param clinseq_barcode_tuple: A clinseq barcode string
    :return: UniqueCapture named tuple
    """

    clinseq_barcode_tuple = parse_clinseq_barcode(clinseq_barcode)

    return UniqueCapture(clinseq_barcode_tuple.project,
                         clinseq_barcode_tuple.sdid,
                         clinseq_barcode_tuple.sample_type,
                         clinseq_barcode_tuple.sample_id,
                         extract_kit_id(clinseq_barcode_tuple.prep_id),
                         extract_kit_id(clinseq_barcode_tuple.capture_id))


def parse_project(clinseq_barcode):
    """
    Extract the project string from the specified clinseq barcode.

    :param clinseq_barcode: Dash-delimited clinseq barcode string.
    :return: The project field from the input string.
    """

    return clinseq_barcode.split("-")[0]


def compose_sample_str(capture):
    """
    Produce a string for a unique library capture item.

    :param capture: A named tuple identifying a unique sample library capture.
    :return: A dash-delimted string of the fields uniquely identifying the capture.
    """
    return "{}-{}-{}-{}".format(capture.sample_type,
                                capture.sample_id,
                                capture.library_kit_id,
                                capture.capture_kit_id)


def parse_sample_type(clinseq_barcode):
    """
    Extract the sample type from the clinseq barcode.

    :param clinseq_barcode: Dash-delimited clinseq barcode string.
    :return: The sample type field from the input string.
    """

    return clinseq_barcode.split("-")[3]


def parse_sample_id(clinseq_barcode):
    """
    Extract the sample ID from the clinseq barcode.

    :param clinseq_barcode: Dash-delimited clinseq barcode string.
    :return: The sample ID field from the input string.
    """

    return clinseq_barcode.split("-")[4]


def parse_sdid(clinseq_barcode):
    """
    Extract the SDID from the clinseq barcode, including the "P-" prefix.

    :param clinseq_barcode: Dash-delimited clinseq barcode string. 
    :return: The SDID field from the input string.
    """

    return "-".join(clinseq_barcode.split("-")[1:3])


def parse_prep_id(clinseq_barcode):
    """
    Extract the library prep ID (the entire string - not just the kit ID) from
    the specified clinseq barcode.

    :param clinseq_barcode: Dash-delimited clinseq barcode string. 
    :return: Library prep ID string.
    """

    return clinseq_barcode.split("-")[5]


def parse_capture_id(clinseq_barcode):
    """
    Extract the library capture ID (the entire string - not just the capture kit ID)
    from the specified clinseq barcode.

    :param clinseq_barcode: Dash-delimited clinseq barcode string. 
    :return: Library capture ID string.
    """

    return clinseq_barcode.split("-")[6]


def extract_kit_id(kit_string):
    """
    Extract the kit type from a specified kit string (either library or capture kit).

    :param kit_string: A string indicating a kit type, where the first two letters indicate the kit ID. 
    :return: The kit ID, comprising the first two letters.
    """

    if len(kit_string) < 3:
        raise ValueError("Invalid kit string: " + kit_string)

    return kit_string[:2]


def clinseq_barcode_is_valid(clinseq_barcode):
    """
    Test the structure of the specified clinseq barcode for validity.

    :param clinseq_barcode: The input clinseq barcode.
    :return: True if the barcode has valid structure, False otherwise.
    """

    # FIXME: Need to implement more stringent checking here.
    fields = clinseq_barcode.split("-")
    if len(fields) != 7:
        return False

    return True


def extract_clinseq_barcodes(input_filename):
    """
    Extrat clinseq barcodes from the specified input file:

    :param input_filename: Either a .txt listing clinseq barcodes one per line,
    or a .xlsx order form file containing the barcodes.

    :return: A list of (not-yet validated) dash-delimited clinseq barcodes.
    """

    toks = input_filename.split(".")
    if len(toks) < 1:
        raise ValueError("Invalid clinseq barcodes input filename: " + input_filename)

    if toks[-1] == "txt":
        [line.strip() for line in open(input_filename)]
    elif toks[-1] == "xlsx":
        return parse_orderform(input_filename)
    else:
        raise ValueError("Invalid clinseq barcodes file type: " + input_filename)


def extract_sample_dict_from_clinseq_barcodes(clinseq_barcodes):
    """
    Extracts a sample dictionary barcode to a nested sample dictionary structure.

    :param clinseq_barcode: A valid clinseq barcode string.
    :return: A dictionary containing the 
    """

    # Extract ClinseqBarcode named tuples from the input strings, also checking
    # them for validity:
    clinseq_barcode_tups = [parse_clinseq_barcode(clinseq_barcode)
                            for clinseq_barcode in clinseq_barcodes]

    # Munge the list of named tuples into a dictionary of the desired structure:
    return convert_barcodes_to_sampledict(clinseq_barcode_tups)


def create_scaffold_sampledict(sdids):
    """
    Generate a scaffold sample dictionary, with SDIDs as keys and clinseq barcode information
    dictionaries as values.

    :param sdids: List of SDID strings
    :return: Dictionary with SDID keys and empty clinseq barcode information dictionaries as values.
    """

    scaffold_dict = {}
    for sdid in sdids:
        curr_clinseq_barcode_info = {"sdid":sdid, "N": [], "T": [], "CFDNA": []}
        scaffold_dict[sdid] = curr_clinseq_barcode_info

    return scaffold_dict


def populate_clinseq_barcode_info(clinseq_barcode_info, clinseq_barcode_tup):
    """
    Populates the specified clinseq barcode information item with the information
    in the specified clinseq barcode.

    :param clinseq_barcode_info: A dictionary containing clinseq barcode information
    for a single SDID.
    :param clinseq_barcode_tup: A ClinseqBarcode named tuple.
    """

    # Append this clinseq barcode XXX CONTINUE HERE: FIGURE OUT HOW TO DO THIS; SHOULD I
    # HAVE BEEN WORKING WITH A STRING ALL ALONG? WHY DID I CONVERT TO A TUPLE IN THE FIRST PLACE?
    # IF NOT, THEN I NEED ANOTHER UTILITY FUNCTION TO CONVERT BACK FROM TUPLE TO STRING.


def convert_barcodes_to_sampledict(clinseq_barcode_tups):
    """
    Munges the list of named tuples into a dictionary linking from SDID to clinseq barcode
    information for that individual.

    :param clinseq_barcode_tups: A list of ClinseqBarcode named tuples.
    :return: A dictionary with the required structure.
    """

    # Extract set of unique SDIDs from the specified clinseq barcodes:
    sdids = set([clinseq_barcode_tup.sdid for clinseq_barcode_tup in clinseq_barcode_tups])

    # Create a scaffold dictionary from those SDIDs:
    sdid_to_clinseq_barcode_info = create_scaffold_sampledict(sdids)

    for clinseq_barcode_tup in clinseq_barcode_tups:
        populate_clinseq_barcode_info(sdid_to_clinseq_barcode_info[clinseq_barcode_tup.sdid],
                                      clinseq_barcode_tup)

    return sdid_to_clinseq_barcode_info

# XXX PROBABLY JUST DELETE ALL THIS:
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
