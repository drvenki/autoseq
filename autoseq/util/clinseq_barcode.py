import collections
from autoseq.util.orderform import parse_orderform

# Fields defining a unique library capture item:
UniqueCapture = collections.namedtuple(
    'sample_type',
    'sample_id',
    'library_kit_id',
    'capture_kit_id'
)


def parse_capture_tuple(clinseq_barcode):
    """
    Convenience function for use in the context of joint panel analysis.

    Extracts the sample type, sample ID, library prep ID, and capture kit ID,
    from the specified clinseq barcode.

    :param clinseq_barcode: List of one or more clinseq barcodes 
    :return: (sample type, sample ID, capture kit ID) named tuple
    """
    return UniqueCapture(parse_sample_type(clinseq_barcode),
                         parse_sample_id(clinseq_barcode),
                         parse_prep_kit_id(clinseq_barcode),
                         parse_capture_kit_id(clinseq_barcode))


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

    return clinseq_barcode.split("-")[1:3]


def parse_prep_kit_id(clinseq_barcode):
    """
    Extract the library prep kit code from the clinseq barcode.

    :param clinseq_barcode: Dash-delimited clinseq barcode string.
    :return: The library prep. kit code extracted from the library prep
    field of the input string.
    """
    return clinseq_barcode.split("-")[5][:2]


def parse_capture_kit_id(clinseq_barcode):
    """
    Extract the capture kit code from the clinseq barcode.

    :param clinseq_barcode: Dash-delimited clinseq barcode string.
    :return: The capture kit code extracted from the panel capture
    field of the input string.
    """
    return clinseq_barcode.split("-")[6][:2]


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

    :return: A list of validated dash-delimite clinseq barcodes.
    """

    toks = input_filename.split(".")
    if len(toks) < 1:
        raise ValueError("Invalid clinseq barcodes input filename: " + input_filename)

    if toks[-1] == "txt":
        return [line.strip() for line in open(input_filename)
                if clinseq_barcode_is_valid(line.strip())]
    elif toks[-1] == "xlsx":
        return parse_orderform(input_filename)
    else:
        raise ValueError("Invalid clinseq barcodes file type: " + input_filename)
