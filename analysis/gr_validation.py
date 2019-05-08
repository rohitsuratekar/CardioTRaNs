"""
CardioTRaNs 2019
Author: Rohit Suratekar

Functions related to the validation of ZFIN data with DESeq data of
Pawlak et. al. 2019.
"""

from pandas import DataFrame

from analysis.basic import get_genes_for_hours, get_genes_from_structure
from helpers.parsers import *


def compare(zfin_df: DataFrame, gr_df: DataFrame):
    gr_filtered = gr_df[COL_GR_S1_GENE_SYMBOL].dropna()
    gr_filtered = list(gr_filtered.reset_index(drop=True))
    zfin_df = zfin_df[COL_EXP_1_GENE_SYMBOL].values
    extra_genes = []
    for gene in list(zfin_df):
        if gene not in gr_filtered:
            extra_genes.append(gene)

    for x in set(extra_genes):
        print(x)

    print(len(set(extra_genes)))


def run():
    zfin_data = get_expression_data_frame()
    zfin_data = get_genes_from_structure(zfin_data, ["heart", "cardi"])
    # zfin_data = get_genes_for_hours(zfin_data, 0, 30000)
    gr_data = get_all_gr_expressed_data()
    compare(zfin_data, gr_data)
