#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# All control visualizations

import matplotlib.pyplot as plt
from SecretColors import Palette

from constants.ngs import *
from constants.system import GENOTYPE_WT
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


def run():
    plot_tpm_disparity(GENOTYPE_WT, [24, 48, 72], BIO_PROJECT_WINATA_LAB)
