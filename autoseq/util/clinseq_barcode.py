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


def extract_unique_capture(clinseq_barcode):
    """
    Parses the specified clinseq barcode and produces a corresponding
    UniqueCapture representing the unique library capture corresponding
    to this clinseq barcode.

    :param clinseq_barcode_tuple: A clinseq barcode string
    :return: UniqueCapture named tuple
    """

    if not clinseq_barcode_is_valid(clinseq_barcode):
        raise ValueError("Invalid clinseq barcode: " + clinseq_barcode)

    project = parse_project(clinseq_barcode)
    sdid = parse_sdid(clinseq_barcode)
    sample_type = parse_sample_type(clinseq_barcode)
    sample_id = parse_sample_id(clinseq_barcode)
    prep_id = parse_prep_id(clinseq_barcode)
    capture_id = parse_capture_id(clinseq_barcode)

    return UniqueCapture(project,
                         sdid,
                         sample_type,
                         sample_id,
                         extract_kit_id(prep_id),
                         extract_kit_id(capture_id))


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


def validate_clinseq_barcodes(clinseq_barcodes):
    """Checks all the specified clinseq barcodes for validity and raises an
    Exception if any are not.
    
    :param clinseq_barcodes: List of clinseq barcode strings.
    """

    # Check the input clinseq barcodes for validity:
    for clinseq_barcode in clinseq_barcodes:
        if not clinseq_barcode_is_valid(clinseq_barcode):
            raise ValueError("Invalid clinseq barcode: " + clinseq_barcode)


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


def populate_clinseq_barcode_info(clinseq_barcode_info, clinseq_barcode):
    """
    Populates the specified clinseq barcode information item with the information
    in the specified clinseq barcode.

    :param clinseq_barcode_info: A dictionary containing clinseq barcode information
    for a single SDID.
    :param clinseq_barcode: A validated clinseq barcode string.
    """

    # Append this clinseq barcode to the relevant field in the specified clinseq barcode
    # info dictionary:
    sample_type = parse_sample_type(clinseq_barcode)
    clinseq_barcode_info[sample_type].append(clinseq_barcode)


def convert_barcodes_to_sampledict(clinseq_barcodes):
    """
    Coverts the specified clinseq barcode strings into a dictionary linking
    from SDID to clinseq barcode information for that individual.

    :param clinseq_barcode_tups: A list of ClinseqBarcode named tuples.
    :return: A dictionary with the required structure.
    """

    # Extract set of unique SDIDs from the specified clinseq barcodes:
    sdids = set([parse_sdid(clinseq_barcode) for clinseq_barcode in clinseq_barcodes])

    # Create a scaffold dictionary from those SDIDs:
    sdid_to_clinseq_barcode_info = create_scaffold_sampledict(sdids)

    for clinseq_barcode in clinseq_barcodes:
        populate_clinseq_barcode_info(sdid_to_clinseq_barcode_info[parse_sdid(clinseq_barcode)],
                                      clinseq_barcode)

    return sdid_to_clinseq_barcode_info
