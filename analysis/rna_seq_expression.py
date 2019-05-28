"""
CardioTRaNs 2019
Author: Rohit Suratekar

All functions related to analysis of RNA-seq data
"""
import matplotlib
import matplotlib.pylab as plt
import numpy as np
from SecretColors.palette import Palette, ColorMap

from helpers.parsers import get_rna_seq_gene_expression_data


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

    c = ColorMap(matplotlib, p).calm()
    dt = plt.pcolor(data, cmap=c)
    plt.gca().set_xticks([])
    t = np.linspace(min(values), max(values), 10)
    color_bar = plt.colorbar(dt, ticks=t)
    tick_labels = [convert_format("{0:.2e}".format(pow(10, x) * min_v)) for x
                   in t]
    color_bar.ax.set_yticklabels(tick_labels)
    plt.show()


def run():
    check_expressing_genes()
