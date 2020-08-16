#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#

import pandas as pd

from helper.wrappers import OutRNA
from helper.additional import get_runs, get_file_names


def average_dataframes(dfs,
                       index_col: str = OutRNA.GENE_NAME,
                       joining_method: str = "outer",
                       tpm_col: str = OutRNA.TPM) -> pd.DataFrame:
    """
    Takes list of Pandas.DataFrames as an input and averages its TPM values

    If you have created DataFrames using 'OutRNA' wrapper, your DataFrame
    will already have one default column as 'gene_name'. If you are using
    third party or manually created DataFrames then you should have some
    common column name among themselves and provide it with 'index_col'
    parameter.

    Additional Note: If you are using your own DataFrame, please remove the
    duplicates before sending them to this function. Otherwise it will give
    error while concatenating the DataFrames.

    :param joining_method: Method to join the DataFrames (default: outer)
    :param dfs: Pandas.DataFrame
    :param index_col: Column which will be used for indexing.
    :param tpm_col: Name of column which will be used for averaging
    :return: First DataFrame with tpm_col replaced with average TPM
    """
    tpms = [x.set_index(index_col)[[tpm_col]] for x in dfs]
    tpms = pd.concat(tpms, join=joining_method, axis=1)
    tpms[OutRNA.AVG_TPM] = tpms.mean(axis=1)
    new_df = pd.concat([dfs[0].set_index(index_col), tpms[[
        OutRNA.AVG_TPM]]], axis=1, join=joining_method)
    new_df = new_df.reset_index()
    del new_df[tpm_col]
    new_df[tpm_col] = new_df[OutRNA.AVG_TPM]
    del new_df[OutRNA.AVG_TPM]
    return new_df


def run():
    runs = get_runs("PRJNA492280")
    method = "salmon"
    files = get_file_names(runs, method, 72, "wt")
    dfs = []
    for f in files:
        dfs.append(OutRNA(f, method).dataframe)

    df = average_dataframes(dfs)
    print(df)
