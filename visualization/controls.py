#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# All control visualizations

import matplotlib.pyplot as plt
import numpy as np
from SecretColors import Palette

from constants.other import *
from constants.ngs import *
from helpers.atlas_parser import atlas_data
from helpers.ngs_parser import compare_stringtie_salmon

p = Palette()


def plot_tpm_disparity(genotype, time_series, bioproject):
    string = "string"
    colors = p.get_color_list

    for i, time in enumerate(time_series):
        data = compare_stringtie_salmon(genotype, time, bioproject)
        data[string] = data[OUTPUT_STRING_TIE] / (
                data[OUTPUT_SALMON] + data[OUTPUT_STRING_TIE])
        data = data.dropna()

        v = plt.violinplot(data[string].values, positions=[i], showmeans=True)

        for body in v["bodies"]:
            body.set_facecolor(colors[i])
            body.set_alpha(0.9)

        v["cmaxes"].set_color(p.black())
        v["cbars"].set_color(p.black())
        v["cmins"].set_color(p.black())
        v["cmeans"].set_color(p.black())

    plt.xticks(range(len(time_series)), time_series)
    plt.xlim([-1, len(time_series)])
    plt.axhline(0.5, color=p.gray(), ls="--", zorder=0)
    plt.ylabel(" TPM $( \\frac{StringTie}{StringTie + Salmon} )$")
    plt.xlabel("hours post fertilization")
    plt.title("'{}' in Bioproject {}".format(genotype, bioproject))
    plt.tight_layout()
    plt.show()


def _get_gene_expression(data, gene):
    data = data[data[STRING_GENE_NAME] == gene]
    if len(data) != 1:
        raise Exception(
            "Either gene {} has multiple values or data does not "
            "exists".format(gene))
    return data


def check_tpm_values_per_gene(genotype, timeseries, bioproject, genes):
    values = []
    for time in timeseries:
        data = compare_stringtie_salmon(genotype, time, bioproject)
        temp = []
        for g in genes:
            tpms = _get_gene_expression(data, g)
            temp.append([tpms[OUTPUT_STRING_TIE].values[0],
                         tpms[OUTPUT_SALMON].values[0]])

        values.append(temp)

    values = np.asarray(values)
    for v in values:
        plt.imshow(v)

    plt.colorbar()
    plt.show()


def check_atlas_expression_data():
    d = atlas_data()
    del d[EXP_ATLAS_GENE_ID]
    d = d.set_index(EXP_ATLAS_GENE_NAME).fillna(0)
    print(d)


def run():
    check_atlas_expression_data()
