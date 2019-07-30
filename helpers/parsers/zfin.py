#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 29/07/19, 4:52 PM
#
#  Copyright (c) 2019.
#
# ZFIN file parsing

import pandas as pd

from constants.system import *


def _get_file_data(filename: str) -> pd.DataFrame:
    """
    Gets filename and prepossess it to remove unnecessary headers
    :param filename: Name of the file without data folder
    :return: Pandas DataFrame of the same data
    """
    with open(DATA_FOLDER + filename) as f:
        k = f.readline()
        if len(k.split("\t")) != 1:
            f.seek(0)

        df = pd.read_csv(f, delimiter="\t")
        df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
        return df


def get_wt_expression() -> pd.DataFrame:
    """
    Get wild type expression data
    :return: Pandas DataFrame
    """
    return _get_file_data(FILE_ZFIN_EXPRESSION)


def get_anatomy_items() -> pd.DataFrame:
    """
    Gets anatomy items data
    :return: Pandas DataFrame
    """
    return _get_file_data(FILE_ZFIN_ANATOMY_ITEMS)


def get_anatomy_relations() -> pd.DataFrame:
    """
    Gets anatomy relationships between the anatomy items
    :return: Pandas DataFrame
    """
    return _get_file_data(FILE_ZFIN_ANATOMY_RELATIONSHIP)


def get_stage_ontology() -> pd.DataFrame:
    """
    Gets gene ontology file
    :return: Pandas DataFrame
    """
    return _get_file_data(FILE_ZFIN_STAGE_ONTOLOGY)


def get_xpat_anatomy() -> pd.DataFrame:
    """
    Gets XPAT anatomy file
    :return: Pandas DataFrame
    """
    return _get_file_data(FILE_ZFIN_XPAT_STAGE)


def run():
    k = get_xpat_anatomy()
    print(k.columns)
