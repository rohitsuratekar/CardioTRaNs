#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
#  ZFIN file parsers

import pandas as pd

from constants.system import FILE_ZFIN_DATA


def _get_file_list(term):
    with open(FILE_ZFIN_DATA) as f:
        m = pd.read_csv(f, delimiter="\t")
    m = m[m["term"] == term]
    if len(m) != 1:
        raise Exception(
            "No file or multiple files found with term {}".format(term))
    else:
        return "data/zfin/{}".format(m["file_name"].values[0])


def _parse_csv(filename: str):
    with open(filename) as f:
        while True:
            loc = f.tell()
            k = f.readline()
            if not k.startswith("#") and len(k.split("\t")) > 1:
                break
        f.seek(loc)
        return pd.read_csv(f, delimiter="\t")


def get_anatomy():
    return _parse_csv(_get_file_list("anatomy"))


def get_expression():
    return _parse_csv(_get_file_list("expression"))


def get_ontology():
    return _parse_csv(_get_file_list("ontology"))


def run():
    f = get_ontology()
    print(f.head())
