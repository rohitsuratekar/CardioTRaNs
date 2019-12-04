#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
#  File-management system

import pandas as pd
import dask.dataframe as dask
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


def _parse(filename, deliminator, use_dask):
    dm = _get_deliminator(deliminator)
    if use_dask:
        if dm is None:
            return dask.read_csv(filename, delim_whitespace=True)
        return dask.read_csv(filename, deliminator=dm)
    with open(filename) as f:
        while True:
            loc = f.tell()
            k = f.readline()
            skip = len(k.split(dm)) > 1
            if dm is None:
                skip = True
            if not k.startswith("#") and skip:
                break
        f.seek(loc)
        if dm is None:
            return pd.read_csv(f, delim_whitespace=True)
        return pd.read_csv(f, delimiter=dm)


def get(term: str, database: str, status: str = "original", use_dask=False):
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

    return _parse(path, m[FILE_DELIMITER].values[0], use_dask)


def get_string(term: str, status, organism, use_dask=False):
    if organism == ORG_ZEBRAFISH:
        return get(term, "string", status, use_dask)
    elif organism == ORG_HUMAN:
        return get(f"h{term}", "string", status, use_dask)
    elif organism == ORG_MOUSE:
        return get(f"m{term}", "string", status, use_dask)
    elif organism == ORG_RAT:
        return get(f"r{term}", "string", status, use_dask)
    else:
        raise Exception(f"Unknown organism {organism}")


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


def zfin_ortho(organism):
    if organism == ORG_HUMAN:
        return get("hortho", "zfin")
    elif organism == ORG_MOUSE:
        return get("mortho", "zfin")
    else:
        raise Exception(f"Unknown organism {organism}")


def expression_atlas():
    return get("atlas", "atlas")


def biomart_go():
    return get("go", "biomart")


def biomart_genes():
    return get("genes", "biomart")


def string_links(status: str = "original", organism: str = ORG_ZEBRAFISH,
                 use_dask=False):
    return get_string("links", status, organism, use_dask)


def string_actions(status: str = "original", organism: str = ORG_ZEBRAFISH,
                   use_dask=False):
    return get_string("actions", status, organism, use_dask)


def string_info(organism: str = ORG_ZEBRAFISH):
    return get_string("info", "original", organism)


def hpa_pathology():
    return get("pathology", "hpa")


def run():
    # d = string_links("temp")
    with open("data/string/links.temp") as f:
        for line in f:
            print(repr(line))
