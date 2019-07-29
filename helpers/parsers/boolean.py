#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 29/07/19, 10:00 AM
#
#  Copyright (c) 2019.
#
# All Boolean data access functions will go here

from pathlib import Path

import pandas as pd

from constants.boolean import *
from constants.ngs import *
from constants.system import DATA_FOLDER, FILE_BOOL_RNA_SEQ
from helpers.parsers.ngs import get_rna_seq_data


def _make_boolean_file(hour: int):
    """
    Generates smaller version of .tsv file so that it can be loaded quickly.
    This will include only genes present in INTERESTED_GENES constant
    :param hour: Hours post fertilization
    """
    d = get_rna_seq_data(hour)
    d = d[d[COL_STRING_TIE_2_GENE_NAME].isin(INTERESTED_GENES)]
    with open("{}{}{}.tsv".format(DATA_FOLDER, FILE_BOOL_RNA_SEQ, hour
                                  ), "w") as f:
        d.to_csv(f, sep="\t", index=False)


def get_boolean_expression(hour: int, override: bool = False) -> pd.DataFrame:
    """
    Generates files used in Boolean Modelling if needed and returns smaller
    pandas DataFrame which can be easily loaded in memory for further analysis.

    :param hour: Hour Post Fertilization
    :param override: If True, then file will be regenerated
    :return: Pandas DataFrame
    """
    path = "{}{}{}.tsv".format(DATA_FOLDER, FILE_BOOL_RNA_SEQ, hour)
    if not Path(path).is_file() or override:
        _make_boolean_file(hour)
    with open(path) as f:
        return pd.read_csv(f, delimiter="\t")


def run():
    get_boolean_expression(24)
