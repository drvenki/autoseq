
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
                           "AL": "alascca_targets"
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
