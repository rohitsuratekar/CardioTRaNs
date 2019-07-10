"""
Project: CardioTrans
Author: Rohit Suratekar
Year: 2019

All NGS analysis functions should go here.
"""

import pandas as pd

from constants.ngs_constants import *
from constants.system_constants import *
from constants.boolean_constants import *
from helpers.parsers import get_rna_seq_data


def filter_with_gene(data: pd.DataFrame, gene: str,
                     allow_substring: bool = False) -> pd.DataFrame:
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
            data[COL_STRING_TIE_2_GENE_NAME].str.contains(
                gene.strip().lower())].reset_index(drop=True)
    return data[data[COL_STRING_TIE_2_GENE_NAME] == gene.strip()].reset_index(
        drop=True)


def filter_with_tpm(data: pd.DataFrame, tpm_cutoff: float) -> pd.DataFrame:
    """
    Filters pandas DataFrames with TPM as a cutoff
    :param data: pd.DataFrame
    :param tpm_cutoff: cutoff value
    :return: pd.DataFrame
    """
    if data.empty:
        return data

    return data[data[COL_STRING_TIE_9_TPM] > tpm_cutoff].reset_index(drop=True)


def filter_with_coverage(data: pd.DataFrame,
                         coverage_cutoff: float) -> pd.DataFrame:
    """
    Filters pandas DataFrames with coverage as a cutoff
    :param data: pd.DataFrame
    :param coverage_cutoff: cutoff value
    :return: pd.DataFrame
    """
    if data.empty:
        return data

    return data[data[COL_STRING_TIE_7_COVERAGE] > coverage_cutoff].reset_index(
        drop=True)


def filter_with_chromosome(data: pd.DataFrame,
                           chromosomes: object) -> pd.DataFrame:
    """
    Filter pandas DataFrame with chromosome number
    :param data: pd.DaraFrame
    :param chromosomes: list or integer
    :return: pd.DataFrame
    """
    if data.empty:
        return data

    if type(chromosomes) != list:
        chromosomes = [chromosomes]

    temp_df = None
    for c in chromosomes:
        if temp_df is None:
            temp_df = data[data[COL_STRING_TIE_3_REFERENCE] == str(c)]
        else:
            temp_df = temp_df.append(
                data[data[COL_STRING_TIE_3_REFERENCE] == str(c)],
                ignore_index=True)

    return temp_df.reset_index(drop=True)


def get_chromosome_with_gene(data: pd.DataFrame, gene: str) -> pd.DataFrame:
    """
    Returns all genes and their expression on the same chromosome as the query gene
    :param data: pd.DataFrame
    :param gene: Query gene
    :return: pd.DataFrame
    """
    data_temp = filter_with_gene(data, gene.strip())
    if len(data_temp) != 1:
        raise Exception(
            "Something went wrong. Total of '{}' data points found for "
            "current gene.".format(len(data_temp)))
    data = filter_with_chromosome(data,
                                  [data_temp[
                                       COL_STRING_TIE_3_REFERENCE].values[0]])
    return data


def save_average_data():
    """
    This saves average RNA-seq data into separate file so that can be easily
    accessible while performing further analysis.
    """
    data = get_rna_seq_data()
    with open(DATA_FOLDER + FILE_RNA_SEQ, 'w') as f:
        header = []
        for a in data:
            header.append(a)

        header = "\t".join(header)
        print(header, file=f)

        for gene in INTERESTED_GENES:
            for i, d in filter_with_gene(data, gene).iterrows():
                print("\t".join([str(x) for x in d.values]), file=f)


def run():
    save_average_data()
