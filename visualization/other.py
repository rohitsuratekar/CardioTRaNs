#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 27/08/19, 9:09 PM
#
#  Copyright (c) 2019.

from constants.boolean import *
from constants.other import *
from helpers.parsers.other import get_expression_manual
import numpy as np
import matplotlib.pylab as plt
import matplotlib
from SecretColors import Palette, ColorMap
from constants.system import PLOT_FOLDER
from matplotlib.patches import Patch


def plot_manual_genes(genes: list, time: list):
    p = Palette()
    on_color = p.peach(shade=40)
    off_color = p.violet(shade=70)
    cmap = ColorMap(matplotlib, p).from_list([off_color, on_color],
                                             is_qualitative=True)

    data = get_expression_manual()
    data = data[data[COL_MANUAL_2_GENE_NAME].isin(genes)]
    data = data[data[COL_MANUAL_6_TIME].isin(time)]

    exp = []

    for t in time:
        k = []
        for g in genes:
            temp = data[data[COL_MANUAL_2_GENE_NAME] == g]
            temp = temp[temp[COL_MANUAL_6_TIME] == t]
            if len(temp) > 0:
                k.append(1)
            else:
                k.append(0)
        exp.append(k)

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
              title="Manual")

    plt.tight_layout()

    plt.savefig(PLOT_FOLDER + "manual_expression.png", dpi=300,
                type="png")
    plt.show()


def run():
    plot_manual_genes(INTERESTED_GENES, [24, 48, 72])
