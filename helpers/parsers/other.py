#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 19/08/19, 4:01 PM
#
#  Copyright (c) 2019.
#
# All other parsers will go here

import pandas as pd

from constants.system import *


def get_expression_atlas():
    """
    Parse the Expression atlas data (TPM) data
    :return: Pandas DataFrame
    """
    with open(DATA_FOLDER + FILE_EXPRESSION_ATLAS) as f:
        skip = f.tell()
        while f.readline().startswith("#"):
            skip = f.tell()
        f.seek(skip)
        data = pd.read_csv(f, delimiter="\t", engine='python')
        return data


def run():
    d = get_expression_atlas()
    print(d.max())
