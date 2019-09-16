#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
#  Other file parsers of the project

from constants.system import FILE_EXPRESSION_ATLAS
import pandas as pd
from helpers.parsers.zfin import get_ontology
from constants.zfin import *


# Expression Atlas
def get_expression_atlas():
    with open(FILE_EXPRESSION_ATLAS) as f:
        while True:
            loc = f.tell()
            k = f.readline()
            if not k.startswith("#"):
                break
        f.seek(loc)
        data = pd.read_csv(f, delimiter="\t")
    return data


def run():
    import matplotlib.pyplot as plt
    g = get_expression_atlas()
    print(g)
