import collections

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
