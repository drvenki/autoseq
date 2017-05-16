import logging

from openpyxl import load_workbook
from clinseq_barcode import *


def parse_orderform(order_form_filename):
    """
    Parse clinseq barcodes from the specified order form.

    :param order_form_filename: An excel spreadsheet filename.
    :return: List of clinseq barcodes extracted from the order form.
    """

    workbook = load_workbook(order_form_filename)
    worksheet = workbook.worksheets[0]
    first_idx = None

    libraries = []

    for idx in range(1, 10000):
        cell_value = worksheet.cell(column=1, row=idx).value
        if cell_value == "<SAMPLE ENTRIES>":
            first_idx = idx

        # no need to parse after end tag
        if cell_value == "</SAMPLE ENTRIES>":
            break

        if not cell_value:
            continue

        if first_idx is not None and idx > first_idx:
            clinseq_barcode = str(worksheet.cell(column=1, row=idx).value)

            logging.debug("Parsed library {}".format(clinseq_barcode))
            sample_type = parse_sample_type(clinseq_barcode)
            if sample_type not in ['T', 'N', 'CFDNA']:
                logging.warning(
                    "Unexpected library type detected: {} for library {}".format(sample_type, clinseq_barcode))

            if not clinseq_barcode_is_valid(clinseq_barcode):
                raise ValueError("Invalid clinseq barcode: " + clinseq_barcode)

            libraries.append(clinseq_barcode)

    return libraries
