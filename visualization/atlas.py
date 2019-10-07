#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Atlas visualizations

import copy

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from SecretColors import Palette, ColorMap

from constants.other import *
from helpers.atlas_parser import atlas_data


def _expressed_genes(genes, tpm_cutoff):
    g = copy.copy(genes)
    g[g < tpm_cutoff] = 0
    g[g >= tpm_cutoff] = 1
    return sum(g)


def tpm_threshold_analysis():
    p = Palette(show_warning=False)
    colors = p.red(no_of_colors=10, starting_shade=20)
    cmap = ColorMap(matplotlib, p).from_list(colors)
    d = atlas_data()
    del d[EXP_ATLAS_GENE_ID]
    del d[EXP_ATLAS_GENE_NAME]
    d = d.fillna(0)
    tpm_range = np.linspace(0, 10, 100)

    norm = matplotlib.colors.Normalize(vmin=min(d.columns), vmax=max(
        d.columns))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    for t in d.columns:
        temp = []
        for tpm in tpm_range:
            g = _expressed_genes(d[t].values, tpm)
            temp.append(g * 100 / len(d))
        plt.plot(tpm_range, temp, color=cmap(t / max(d.columns)))

    plt.xlim(0, 10)
    plt.ylim(0, 100)
    plt.colorbar(sm, label="Time (hpf)")
    plt.ylabel("Percentage of genes expressed")
    plt.xlabel("TPM Threshold")
    plt.grid(axis="both")
    plt.tight_layout()
    plt.show()


def plot_tpm_values():
    p = Palette()
    d = atlas_data()
    del d[EXP_ATLAS_GENE_NAME]
    del d[EXP_ATLAS_GENE_ID]
    d = d.fillna(0).values

    plt.hist(d.flatten(), bins=100, color=p.teal())
    plt.yscale("log")
    plt.ylabel("Frequency (on log scale)")
    plt.xlabel("TPM (from all time points)")
    plt.tight_layout()
    plt.show()


def run():
    plot_tpm_values()
