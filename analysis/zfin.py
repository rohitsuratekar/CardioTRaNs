#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
#   Analysis related to ZFIN datasets

import yaml

from helper.wrappers import ZData
from constants.zfin import *


def _get_files():
    with open("config.yml") as f:
        data = yaml.load(f, Loader=yaml.SafeLoader)
    return data["zfin"]


def find_positive_control_genes():
    # Positive control genes for 24-72 hpf time patch
    files = _get_files()
    exp = files["expression"]
    xpat = files["xpat"]
    ont = files["ontology"]
    z = ZData(exp, xpat, ont)
    z.add_start_stage(24)
    z.add_end_stage(72)
    z.add_structures(["heart", "cardi"])
    genes = set(z.dataframe[ZFIN_EXP_GENE_SYMBOL].values)
    print(f"Total of {len(genes)} : {genes}")


def find_negative_control_genes():
    files = _get_files()
    exp = files["expression"]
    xpat = files["xpat"]
    ont = files["ontology"]
    lines = files["fish_lines"]
    f_xpat = files["fish_xpat"]
    z = ZData(exp, xpat, ont)
    z.add_start_stage(20)
    z.add_end_stage(30)
    z.add_structures(["heart", "cardi"])
    z.search_non_expressed(f_xpat, lines)
    genes = set(z.dataframe[ZFIN_EXP_GENE_SYMBOL].values)
    print(f"Total of {len(genes)} : {genes}")


def run():
    find_negative_control_genes()
