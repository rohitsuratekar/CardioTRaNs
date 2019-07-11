"""
Project: CardioTrans
Author: Rohit Suratekar
Year: 2019

This file should contain all functions related to Boolean analysis
"""

import pandas as pd

from constants.boolean_constants import *
from constants.ngs_constants import COL_STRING_TIE_2_GENE_NAME
from helpers.parsers import get_boolean_genes


def get_state(gene: str, df: pd.DataFrame):
    """
    Get the expression state of given gene in the dataframe
    :param gene: name of the gene
    :param df: Pandas DataFrame with BOOL_OBJ_COL column assigned
    :return: True if gene is expressed. [note: change the condition
    of expression from the model properties]
    """
    ob = df[df[COL_STRING_TIE_2_GENE_NAME] == gene][BOOL_OBJ_COL].values
    if len(ob) != 1:
        raise Exception(
            "Either zero or multiple results found for {}".format(gene))
    return ob[0].is_expressed


def get_initial_condition():
    data = get_boolean_genes()
    for g in INTERESTED_GENES:
        print(get_state(g, data))


def run():
    get_initial_condition()
