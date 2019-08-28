#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 02/08/19, 11:14 AM
#
#  Copyright (c) 2019.
#
# All visualization related to ZFIN database

import matplotlib
import matplotlib.pylab as plt
import numpy as np
from SecretColors import Palette, ColorMap
from matplotlib.patches import Patch

from analysis.zfin import *
from constants.zfin import *
from helpers.parsers.zfin import *
from constants.boolean import *


def visualize_genes_expression(data: pd.DataFrame,
                               genes: list,
                               at_time: list = None):
    """
    Visualize the gene expression present in the ZFIN database.

    :param data: Data to visualize
    :param genes: List of genes which will be filtered and plotted
    :param at_time: You can provide list of time points at which final plot
    will be shown
    """

    p = Palette()
    on_color = p.peach(shade=40)
    off_color = p.violet(shade=70)
    cmap = ColorMap(matplotlib, p).from_list([off_color, on_color],
                                             is_qualitative=True)
    d = filter_by_genes(data, genes)

    stages = get_stage_ontology()
    stages = stages[[COL_ONT_3_STAGE_NAME, COL_ONT_4_BEGIN_HOUR]]
    stages = stages[stages[COL_ONT_3_STAGE_NAME] != "Unknown"]

    d = d.groupby([COL_EXP_1_GENE_SYMBOL, COL_EXP_7_START_STAGE]).agg(
        {COL_EXP_0_GENE_ID: "count"})
    d = d.reset_index()
    exp = []
    for g in genes:
        temp = d[d[COL_EXP_1_GENE_SYMBOL] == g]
        k = stages.merge(temp, left_on=COL_ONT_3_STAGE_NAME,
                         right_on=COL_EXP_7_START_STAGE,
                         how="outer")

        k[COL_EXP_0_GENE_ID] = k[COL_EXP_0_GENE_ID] / k[COL_EXP_0_GENE_ID]
        k = k.fillna(0)
        exp.append(k[COL_EXP_0_GENE_ID].values)

    exp = np.asarray(exp)
    # _, ax = plt.subplots(figsize=(15, 10))
    _, ax = plt.subplots()

    time_labels = stages[COL_ONT_4_BEGIN_HOUR].values
    if at_time is not None:
        time_index = []
        time_labels = []
        for t in at_time:
            try:
                time_index.append(
                    list(stages[COL_ONT_4_BEGIN_HOUR].values).index(t))
                time_labels.append(t)
            except ValueError:
                # If given time is not available in the database, get the
                # nearest time point

                temp_time = 0
                for k in stages[COL_ONT_4_BEGIN_HOUR].values:
                    if t > k:
                        temp_time = k
                print("Given value {} not found in the database, "
                      "hence nearest {} value is used".format(t, temp_time))
                time_index.append(
                    list(stages[COL_ONT_4_BEGIN_HOUR].values).index(temp_time))
                time_labels.append(temp_time)

        exp = exp[:, time_index]

    plt.pcolormesh(exp, cmap=cmap, edgecolors=p.gray(shade=70),
                   linewidth=0.05, vmax=1, vmin=0)

    ax.set_xticks(np.arange(0, len(time_labels)) + 0.5)
    ax.set_yticks(np.arange(0, len(genes)) + 0.5)
    ax.set_xticklabels(time_labels, rotation=45, ha="left")
    ax.set_yticklabels([str(x).capitalize() for x in genes])
    plt.xlabel("Hours Post Fertilization")
    legend_elements = [Patch(facecolor=on_color, label='ON'),
                       Patch(facecolor=off_color, label='OFF')]
    plt.gca().legend(handles=legend_elements, loc='upper center',
                     bbox_to_anchor=(0.5, -0.05), fancybox=True, ncol=2)
    plt.gca().set_aspect(1.2)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    plt.tight_layout()
    plt.savefig(PLOT_FOLDER + "zfin_gene_expression.png", dpi=300, type="png")
    plt.show()


def plot_structure_wise_genes(data: pd.DataFrame, top: int = 10):
    p = Palette()
    d = data.groupby(COL_EXP_4_SUPER_STR_NAME).count()[[COL_EXP_0_GENE_ID]]
    d = d.sort_values(by=COL_EXP_0_GENE_ID,
                      ascending=False).reset_index().head(top)

    names = d[COL_EXP_4_SUPER_STR_NAME].values
    values = d[COL_EXP_0_GENE_ID].values
    ind = range(len(names))
    plt.barh(ind, values, color=p.blue())
    plt.yticks(ind, names)
    plt.ylabel("ZFIN Super Structures")
    plt.title("Top {} structures".format(top))
    plt.xlabel("Number of genes expressed")
    plt.tight_layout()
    plt.savefig(PLOT_FOLDER + "zfin_structure.png", dpi=300, type="png")
    plt.show()


def zfin(data, genes: list, time: list, structure: list):
    # TODO : Check properly
    p = Palette()
    on_color = p.peach(shade=40)
    off_color = p.violet(shade=70)
    cmap = ColorMap(matplotlib, p).from_list([off_color, on_color],
                                             is_qualitative=True)

    data = filter_by_genes(data, genes)
    data = filter_by_structure(data, structure, exact=False)
    exp = []
    for t in time:
        temp = filter_by_start_time(data, t)
        all_genes = list(set(temp[COL_EXP_1_GENE_SYMBOL].values))
        k = []
        for g in genes:
            k.append(g in all_genes)

        exp.append([int(x) for x in k])

    exp = np.asarray(exp)

    fig, ax = plt.subplots()
    plt.pcolormesh(exp.T, edgecolor=p.gray(shade=70),
                   linewidth=0.05, cmap=cmap, vmin=0, vmax=1)
    ax.set_aspect(1.2)
    ax.set_xticks(np.arange(0, len(time)) + 0.5)
    ax.set_yticks(np.arange(0, len(genes)) + 0.5)
    ax.set_xticklabels(["{}".format(x) for x in time], rotation=45,
                       ha="left")
    ax.set_yticklabels([str(x).capitalize() for x in genes])
    plt.xlabel("Hours Post Fertilization")
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    # ax.set_title("RNA-seq Boolean Conditions", y=1.08)
    legend_elements = [Patch(facecolor=on_color, label='ON'),
                       Patch(facecolor=off_color, label='OFF or N/A')]
    ax.legend(handles=legend_elements, loc='upper center',
              bbox_to_anchor=(0.5, -0.05), fancybox=True, ncol=2,
              title="with structure {}".format(structure))

    plt.tight_layout()

    plt.savefig(PLOT_FOLDER + "zfin_gene_expression.png", dpi=300,
                type="png")
    plt.show()


def run():
    d = get_wt_expression()
    zfin(d, INTERESTED_GENES, [24, 48, 72], ["heart", "cardi"])
