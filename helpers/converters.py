"""
Project: CardioTrans
Author: Rohit Suratekar
Year: 2019

This converter function.
"""

import pandas as pd

from helpers.parsers import get_zfin_stage_ontology
from models.ngs_models import STGeneExpression
from models.zfin_models import ZFINExpression


def convert_to_zfin_object(data: pd.DataFrame) -> list:
    """
    Converts Pandas dataframe into ZFINExpression objects
    :param data: Pandas.DataFrame
    :return: List of ZFINExpression objects
    """

    exp_list = []

    if data.empty:
        return exp_list

    stages = {s.name: s.id for s in get_zfin_stage_ontology()}

    for d in data.iterrows():
        row = []
        for i in d[1]:
            row.append(str(i))
        exp_list.append(ZFINExpression(row, stages))

    return exp_list


def convert_to_ngs_object(data: pd.DataFrame) -> list:
    """
    Converts the pandas DataFrames into the STGeneExpression objects
    :param data: pandas.DataFrame
    :return: list of STGeneExpression objects
    """
    exp_list = []

    if data.empty:
        return exp_list

    for ind, value in data.iterrows():
        row = []
        for v in value:
            row.append(str(v).strip())
        exp_list.append(STGeneExpression(row))

    return exp_list
