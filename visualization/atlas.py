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
from SecretPlots import *

from analysis.atlas import isolate_by_time, get_go_terms
from constants.biomart import *
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


def plot_time_range(start_time, end_time, tpm_cutoff):
    p = Palette(show_warning=False)
    data = atlas_data()
    data = isolate_by_time(data, start_time=start_time, end_time=end_time,
                           tpm_cutoff=tpm_cutoff).head(40)
    print("Total gene found : {}".format(len(data)))
    del data[EXP_ATLAS_GENE_ID]
    data = data.set_index(EXP_ATLAS_GENE_NAME)
    cp = (ColorPlot(data.values)
          .change_orientation("y")
          .add_x_ticklabels(data.columns.values, rotation=45, ha="left")
          .add_y_ticklabels(data.index.values)
          .add_x_label("Hours post fertilization")
          .add_x_padding(0, 0)
          .add_y_padding(0, 0)
          .add_on_color(p.amber())
          .add_off_color(p.blue(shade=60))
          .invert_y()
          )
    cp.draw()
    cp.ax.xaxis.set_ticks_position("top")
    cp.ax.xaxis.set_label_position("top")
    cp.show(tight=True)


def plot_targeted_go(start, end, tpm):
    by = BIOMART_GO_NAME
    p = Palette()
    no_of_genes = 30
    go = get_go_terms(GO_PROCESS, start, end, tpm_cutoff=tpm).head(no_of_genes)

    bar = (BarPlot(go["count"].values)
           .change_orientation("y")
           .invert_y()
           .add_y_ticklabels(go[by].values)
           .add_x_label("Number of expressed genes (> {} TPM)".format(tpm))
           .add_colors([p.red()])
           )
    bar.draw()
    title = ""
    if len(go) > no_of_genes:
        title = "Top {} ".format(no_of_genes)
    title += "'{}' for genes expressed\nonly between {} -{} hpf".format(by,
                                                                        start,
                                                                        end)
    bar.ax.set_title(title)
    bar.show()


def run():
    plot_targeted_go(5, 10, 1)
