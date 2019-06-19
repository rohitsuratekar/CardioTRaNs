"""
CardioTRaNs 2019
Author: Rohit Suratekar

All functions related to analysis of RNA-seq data
"""

import matplotlib
import matplotlib.pylab as plt
import numpy as np
from SecretColors.palette import Palette, ColorMap
from matplotlib.cm import ScalarMappable

from constants import *
from helpers.parsers import get_rna_seq_gene_expression_data
from models.ngs import *


def check_coverage():
    p = Palette()
    data_tpm = []
    data_coverage = []
    for d in get_rna_seq_gene_expression_data():
        if 20000 > d.tpm:
            data_tpm.append(d.tpm)
            data_coverage.append(d.coverage)

    plt.hist(data_tpm, 100, color=p.blue(shade=50), label="TPM")
    plt.hist(data_coverage, 100, color=p.red(shade=50), label="Coverage")
    plt.yscale("log")
    plt.axhline(1, linestyle="--", color=p.black())
    plt.xlabel("Transcripts Per Million (Blue)/ Coverage (Red) \nClipped at "
               "20000")
    plt.legend(loc=0)
    plt.ylabel("Frequency (Log Scale)")
    plt.show()


def convert_format(st: str) -> str:
    return st.replace("e", " $x10^{") + "}$"


def check_expressing_genes():
    p = Palette()
    genes = {}
    names = []
    for d in get_rna_seq_gene_expression_data():
        if d.tpm > 0:
            genes[d.gene_name] = d.tpm
            names.append(d.gene_name)
    names = sorted(names)
    values = np.asanyarray([genes[x] for x in names])
    min_v = min(values)
    values = np.log10(values / min_v)
    data = []
    for i, v in enumerate(values):
        data.append([v, v])

    c_colors = p.blue(no_of_colors=5)
    c_colors = [x for x in reversed(c_colors)]
    c_colors.extend(p.red(no_of_colors=5))
    c = ColorMap(matplotlib, p).from_list(c_colors)
    dt = plt.pcolor(data, cmap=c)
    plt.gca().set_xticks([])
    t = np.linspace(min(values), max(values), 10)
    color_bar = plt.colorbar(dt, ticks=t)
    tick_labels = [convert_format("{0:.2e}".format(pow(10, x) * min_v)) for x
                   in t]
    color_bar.ax.set_yticklabels(tick_labels)
    plt.show()


def map_to_chromosome(chr_no: int, max_size: int = 100):
    p = Palette()
    data = []
    for d2 in get_rna_seq_gene_expression_data():
        d = d2  # type:STGeneExpression
        if str(d.reference) == str(chr_no):
            data.append(d)

    # Get farthest position
    max_pos = max([x.end for x in data])
    chromosome = [0] * max_size
    for d in data:
        d.start = d.start * max_size / max_pos
        d.end = d.end * max_size / max_pos
        chromosome[int(d.start)] += d.tpm

    pd = []
    for p1 in chromosome:
        pd.append([p1, p1])

    c_colors = p.blue(no_of_colors=5)
    c_colors = [x for x in reversed(c_colors)]
    c = ColorMap(matplotlib, p).from_list(c_colors)
    dt = plt.pcolor(pd, cmap=c)
    color_bar = plt.colorbar(dt)
    y_ticks = np.linspace(0, max_size, 10)
    plt.xticks([])
    plt.yticks(y_ticks,
               ["{:.2f} MB".format(x * max_pos / (max_size * 1000000)) for x in
                y_ticks])
    plt.show()


def map_all_chromosomes_bar(min_exp_percentage=0.01):
    p = Palette()
    data = get_rna_seq_gene_expression_data(True)
    max_pos = data[COL_STRING_TIE_6_END].max()
    y_ticks = data[COL_STRING_TIE_3_REFERENCE].unique()

    # Sort All Chromosome names according to their number
    y_ticks_temp1 = []  # Chromosome numbers
    y_ticks_temp2 = []  # Additional Chromosomes like mitochondrial
    for x in y_ticks:
        try:
            y_ticks_temp1.append(int(x))
        except ValueError:
            y_ticks_temp2.append(x)

    y_ticks_temp1 = sorted(y_ticks_temp1)
    y_ticks = [str(x) for x in y_ticks_temp1]
    y_ticks.extend(y_ticks_temp2)
    data[COL_STRING_TIE_5_START] = data[COL_STRING_TIE_5_START] * 100 / max_pos
    data[COL_STRING_TIE_6_END] = data[COL_STRING_TIE_6_END] * 100 / max_pos

    c_colors = [p.cyan(shade=50), p.white()]
    c = ColorMap(matplotlib, p).from_list(c_colors)

    for ind, cr in enumerate(y_ticks):
        temp_d = data[data[COL_STRING_TIE_3_REFERENCE] == cr]
        temp_end = temp_d[COL_STRING_TIE_6_END].max()
        temp_d[COL_STRING_TIE_9_TPM] /= temp_d[COL_STRING_TIE_9_TPM].max()
        y_loc = (ind, 0.8)
        plt.broken_barh([(0, temp_end)], y_loc, color=p.cyan(shade=80))
        for i, row in temp_d.iterrows():
            start = row[COL_STRING_TIE_5_START]
            end = row[COL_STRING_TIE_6_END]
            if row[COL_STRING_TIE_9_TPM] > min_exp_percentage:
                plt.broken_barh([(start, end - start)], y_loc,
                                color=c(row[COL_STRING_TIE_9_TPM]))

    plt.yticks([x + 0.4 for x in range(len(y_ticks))],
               ["chr {}".format(x) for x in y_ticks])
    sm = ScalarMappable(cmap=c, norm=plt.Normalize(0, 100))
    sm.set_array([])
    cbar = plt.colorbar(sm)
    cbar.set_label("TPM % per chromosome")
    x_ticks = np.linspace(0, 100, 10)
    plt.xticks(x_ticks,
               ["{:.2f} MB".format(x * max_pos / (100 * 1000000)) for x in
                x_ticks])
    plt.xlabel("Position")
    plt.title("TPM greater than {} % of largest per chromosome".format(
        min_exp_percentage))
    plt.show()


def check_gene(name:str):
    for d in get_rna_seq_gene_expression_data():
        if str(d.gene_name).lower() == name.strip().lower():
            print(d.data)
            print(d.tpm)


def run():
    check_gene("tbx5b")
