#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
#   Analysis done on the DESeq2 output files

import pandas as pd

from constants.outputs import *
from helper.wrappers import OutRNA


def extract_constant_genes(df: pd.DataFrame,
                           name_col: str = OutRNA.GENE_NAME) -> list:
    df = df[df[DESEQ2_PADJ] <= 0.05]
    df = df[df[DESEQ2_LOG2_CHANGE] > -1]
    df = df[df[DESEQ2_LOG2_CHANGE] < 1]
    return df[name_col].values


def check_expression(df: pd.DataFrame, genes: list,
                     name_col: str = OutRNA.GENE_NAME):
    df = df[df[name_col].isin(genes)]
    df = df[[name_col, DESEQ2_LOG2_CHANGE, DESEQ2_PADJ]]
    print(df)


def run():
    filename = "/mnt/windows/Enigma/Zebrafish/data/deseq2/mutant/salmon/"
    filename += "salmon_tbx5a_vs_wt.csv"
    df = OutRNA(filename, "deseq2").dataframe
    w = extract_constant_genes(df)
    print(w)
