import logging
import os

from openpyxl import load_workbook

from autoseq.util.library import get_libdict
from autoseq.util.path import stripsuffix


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
            library_name = str(worksheet.cell(column=1, row=idx).value)

            lib = get_libdict(library_name)
            logging.debug("Parsed library {}".format(lib))
            if lib['type'] not in ['T', 'N', 'P']:
                logging.warning(
                    "Unexpected library type detected: {} for library {}".format(lib['type'], lib['library_id']))
            libraries.append(lib)

    return libraries
