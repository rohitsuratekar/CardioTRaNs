#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
#   All analysis related to finding the genes

import pandas as pd

from analysis.deseq2 import extract_constant_genes
from analysis.rnaseq import average_dataframes
from helper.wrappers import OutRNA
from helper.additional import get_file_names, get_runs
from scipy.stats.mstats import gmean


def find_geometric_mean(deseq2: pd.DataFrame,
                        rnaseq: list,
                        tpm_col: str = OutRNA.TPM,
                        name_col: str = OutRNA.GENE_NAME):
    w = extract_constant_genes(deseq2, name_col=name_col)
    for i, df in enumerate(rnaseq):
        df = df[df[name_col].isin(w)]
        df = df[tpm_col].values
        gm = gmean(df)
        print(f"Geometric Mean Sample {i + 1}: {gm} ")

    print(f"Total Genes: {len(w)}")


def _example_frames(genotype):
    runs = get_runs("PRJNA492280")
    method = "salmon"
    files = get_file_names(runs, method, 72, genotype)
    dfs = []
    for f in files:
        dfs.append(OutRNA(f, method).dataframe)
    return dfs


def run():
    genotype = "gata5"
    filename = "/mnt/windows/Enigma/Zebrafish/data/deseq2/mutant/salmon/"
    filename += f"salmon_{genotype}_vs_wt.csv"
    de_df = OutRNA(filename, "deseq2").dataframe
    rna_df = _example_frames(genotype)
    find_geometric_mean(de_df, rna_df)
