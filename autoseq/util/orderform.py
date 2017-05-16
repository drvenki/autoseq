import logging

from openpyxl import load_workbook
from clinseq_barcode import *

def parse_orderform(xlsx):
    # orderform_name = stripsuffix(os.path.basename(xlsx), ".xlsx")
    workbook = load_workbook(xlsx)
    worksheet = workbook.worksheets[0]
    first_idx = None

    libraries = []

    for idx in range(1, 1000):
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

            # FIXME: Validate the clinseq barcodes here.

            logging.debug("Parsed library {}".format(clinseq_barcode))
            sample_type = parse_sample_type(clinseq_barcode)
            if sample_type not in ['T', 'N', 'CFDNA']:
                logging.warning(
                    "Unexpected library type detected: {} for library {}".format(sample_type, clinseq_barcode))
            libraries.append(clinseq_barcode)

    return libraries
