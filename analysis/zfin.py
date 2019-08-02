#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 30/07/19, 3:22 PM
#
#  Copyright (c) 2019.
#
# All functions related to the zfin data analysis will go here

from constants.zfin import *
from helpers.parsers.zfin import *


def filter_by_genes(data: pd.DataFrame, genes: list,
                    allow_substring: bool = False):
    """
    Filters pandas DataFrames with gene symbols
    :param data: pd.DataFrame
    :param genes: List of Gene Symbol
    :param allow_substring: If True, substring matches will also be returned
    :return: pd.DataFrame after filtering
    """

    if data.empty:
        return data

    if allow_substring:
        return data[
            data[COL_EXP_1_GENE_SYMBOL].str.contains(
                "|".join(genes).lower(), )].reset_index(drop=True)
    return data[data[COL_EXP_1_GENE_SYMBOL].isin(genes)].reset_index(
        drop=True)


def filter_by_gene(data: pd.DataFrame, gene: str,
                   allow_substring: bool = False):
    """
    Filters pandas DataFrames with gene symbols
    :param data: pd.DataFrame
    :param gene: Gene Symbol
    :param allow_substring: If True, substring matches will also be returned
    :return: pd.DataFrame after filtering
    """

    if data.empty:
        return data

    if allow_substring:
        return data[
            data[COL_EXP_1_GENE_SYMBOL].str.contains(
                gene.strip().lower())].reset_index(drop=True)
    return data[data[COL_EXP_1_GENE_SYMBOL] == gene.strip()].reset_index(
        drop=True)


def filter_by_start_time(data: pd.DataFrame, hour: float,
                         tolerance: float = None,
                         tolerance_left: float = None,
                         tolerance_right: float = None) -> pd.DataFrame:
    """
    Filters data based on the start time.

    Note: If tolerance value is provided, tolerance_left and tolerance_right
    will not be used

    :param data: Pandas DataFrame of expression data
    :param hour: Hour value to filter
    :param tolerance: If provided, values around hour with given tolerance
    will also be included
    :param tolerance_left: If provided, values BEFORE hour with given
    tolerance_left will be included
    :param tolerance_right: If provided, values AFTER hour with given
    tolerance_right will be included
    :return: Pandas DataFrame after filtering according to above rules
    """
    if data.empty:
        return data

    # Setup tolerance values
    if tolerance_right is None:
        tolerance_right = 0
    if tolerance_left is None:
        tolerance_left = 0
    if tolerance is not None:
        tolerance_left = tolerance
        tolerance_right = tolerance

    # Get stage values
    ontology = get_stage_ontology()
    ontology = ontology[ontology[COL_ONT_4_BEGIN_HOUR] >= hour -
                        tolerance_left]
    ontology = ontology[ontology[COL_ONT_4_BEGIN_HOUR] <= hour
                        + tolerance_right]
    ontology = ontology[COL_ONT_3_STAGE_NAME].values
    d = data[data[COL_EXP_7_START_STAGE].isin(ontology)].reset_index(drop=True)
    return d


def filter_by_end_time(data: pd.DataFrame, hour: float,
                       tolerance: float = None,
                       tolerance_left: float = None,
                       tolerance_right: float = None) -> pd.DataFrame:
    """
    Filters data based on the end time.

    Note: If tolerance value is provided, tolerance_left and tolerance_right
    will not be used

    :param data: Pandas DataFrame of expression data
    :param hour: Hour value to filter
    :param tolerance: If provided, values around hour with given tolerance
    will also be included
    :param tolerance_left: If provided, values BEFORE hour with given
    tolerance_left will be included
    :param tolerance_right: If provided, values AFTER hour with given
    tolerance_right will be included
    :return: Pandas DataFrame after filtering according to above rules
    """
    if data.empty:
        return data

    # Setup tolerance values
    if tolerance_right is None:
        tolerance_right = 0
    if tolerance_left is None:
        tolerance_left = 0
    if tolerance is not None:
        tolerance_left = tolerance
        tolerance_right = tolerance

    # Get stage values
    ontology = get_stage_ontology()
    ontology = ontology[ontology[COL_ONT_5_END_HOUR] >= hour -
                        tolerance_left]
    ontology = ontology[ontology[COL_ONT_5_END_HOUR] <= hour
                        + tolerance_right]
    ontology = ontology[COL_ONT_3_STAGE_NAME].values
    d = data[data[COL_EXP_8_END_STAGE].isin(ontology)].reset_index(drop=True)
    return d


def filter_by_structure(data: pd.DataFrame, structures: list,
                        exact: bool = True):
    """
    Filters data set with 'super structure' values in the column

    :param data: Pandas DataFrame with expression data
    :param structures: List of structures
    :param exact: If True, exact name match will be done, else substring
    will be searched against structure names
    :return: Pandas DataFrame
    """
    if data.empty:
        return data

    if not exact:
        s = [str(x).lower() for x in structures]
        s = "|".join(s)

        d = data[data[COL_EXP_4_SUPER_STR_NAME].str.contains(s)].reset_index(
            drop=True)

        return d

    return data[data[COL_EXP_4_SUPER_STR_NAME].isin(structures)].reset_index(
        drop=True)


def run():
    d = get_wt_expression()
    d = filter_by_structure(d, ["heart", "cardi"], exact=False)
    d = filter_by_start_time(d, 24, tolerance=0)
    d = filter_by_gene(d, "nkx2.5")
    print(d)
