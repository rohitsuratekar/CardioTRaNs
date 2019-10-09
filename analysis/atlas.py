#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Analysis related to the Expression Atlas

import functools
import operator

from constants.biomart import *
from constants.other import *
from helpers.atlas_parser import atlas_data
from helpers.filemanager import biomart_go


def isolate_by_time(data, start_time: float, end_time: float,
                    tpm_cutoff: float,
                    strict: bool = True,
                    not_expressed: bool = False,
                    ignore: list = None):
    """
    Isolates the Atlas dataframe based on various criteria

    :param data: Pandas DataFrame for Expression Atlas
    :param start_time: Starting time point (hpf)
    :param end_time: Ending time point (hpf)
    :param tpm_cutoff: TPM Value above which gene will be considered as
    expressed
    :param strict: If True, all other time points should have opposite state
     of gene expression
    :param not_expressed: If True, it will look for gene which are NOT
    expressed in given time range
    :param ignore: List of time points which will be ignore if analysis is
    running in 'strict' mode
    :return: Pandas DataFrame after filtering
    """

    if ignore is None:
        ignore = []

    d = data.set_index([EXP_ATLAS_GENE_NAME, EXP_ATLAS_GENE_ID]).fillna(0)
    for time in d.columns:
        if start_time <= time <= end_time:
            if not_expressed:
                d = d[d[time] < tpm_cutoff]
            else:
                d = d[d[time] >= tpm_cutoff]
        else:
            if strict and time not in ignore:
                if not_expressed:
                    d = d[d[time] >= tpm_cutoff]
                else:
                    d = d[d[time] < tpm_cutoff]

    d["temp"] = d.sum(axis=1)
    d = d.sort_values(by="temp", ascending=False)
    del d["temp"]
    d = d.reset_index()
    return d


def general_statistics():
    """
    Generates basic statistics for expression atlas
    :return:
    """
    d = atlas_data()
    print("Number of genes : {}".format(len(d)))
    print("Number of time points : {}".format(len(d.columns) - 2))
    del d[EXP_ATLAS_GENE_ID]
    del d[EXP_ATLAS_GENE_NAME]


def go_distribution(by, start_time, end_time, tpm_cutoff, strict=False):
    go = biomart_go()
    atlas = atlas_data()
    atlas = (isolate_by_time(atlas,
                             start_time,
                             end_time,
                             tpm_cutoff,
                             strict=strict)[EXP_ATLAS_GENE_ID].values)

    go = (go[go[BIOMART_GENE_ID]
          .isin(atlas)][[by, BIOMART_GENE_NAME]]
          .dropna()
          .groupby(by)
          .agg([(BIOMART_GENE_NAME, lambda x: [y for y in x]),
                ("count", "count")])
          .sort_values((BIOMART_GENE_NAME, "count"), ascending=False))

    # Flatten the multi index columns
    go.columns = go.columns.map(lambda x: x[1] if len(x[1]) != 0 else x[0])
    go = go.reset_index()
    return go


def get_annotated_genes(start, end, tpm):
    d = go_distribution(BIOMART_GO_NAME, start, end, tpm, strict=True)
    d = d[BIOMART_GENE_NAME].values
    d = list(set(functools.reduce(operator.iconcat, d, [])))
    print(d)


def run():
    get_annotated_genes(5, 12, 1)
