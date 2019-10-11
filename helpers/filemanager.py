#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
#  File-management system

import pandas as pd

from constants.other import *
from constants.system import FILE_MANAGER, DATA_FOLDER


def _manager():
    with open(FILE_MANAGER)as f:
        return pd.read_csv(f, delimiter=",")


def _get_deliminator(symbol):
    if symbol == "t":
        return "\t"
    elif symbol == "c":
        return ","
    elif symbol == "s":
        return None
    else:
        raise Exception("No such deliminator symbol found : {}".format(symbol))


def _parse(filename, deliminator):
    with open(filename) as f:

        while True:
            loc = f.tell()
            k = f.readline()
            dm = _get_deliminator(deliminator)

            skip = len(k.split(dm)) > 1
            if dm is None:
                skip = True
            if not k.startswith("#") and skip:
                break
        f.seek(loc)
        if dm is None:
            return pd.read_csv(f, delim_whitespace=True)

        return pd.read_csv(f, delimiter=_get_deliminator(
            deliminator))


def get(term: str, database: str, status: str = "original"):
    if status != "original":
        print("You are working with file status :{}".format(status))
    m = _manager()
    m = m[(m[FILE_TERM] == term)
          & (m[FILE_DATASET] == database)
          & (m[FILE_STATUS] == status)]

    if len(m) != 1:
        raise Exception(
            "No such dataset or multiple results for {} in {}".format(term,
                                                                      database))

    path = "{}{}{}".format(DATA_FOLDER,
                           m[FILE_FOLDER].values[0],
                           m[FILE_NAME].values[0])

    return _parse(path, m[FILE_DELIMITER].values[0])


def zfin_expression():
    return get("expression", "zfin")


def zfin_ontology():
    return get("ontology", "zfin")


def zfin_anatomy():
    return get("anatomy", "zfin")


def zfin_pattern():
    return get("pattern", "zfin")


def zfin_assay():
    return get("assay", "zfin")


def zfin_wt_lines():
    return get("fish", "zfin")


def expression_atlas():
    return get("atlas", "atlas")


def biomart_go():
    return get("go", "biomart")


def biomart_genes():
    return get("genes", "biomart")


def string_links(status: str = "original"):
    return get("links", "string", status)


def string_actions(status: str = "original"):
    return get("actions", "string", status)


def string_info():
    return get("info", "string")


def run():
    d = string_info()
    print(d.columns)
