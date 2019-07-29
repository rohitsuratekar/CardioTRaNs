#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 29/07/19, 2:40 PM
#
#  Copyright (c) 2019.
#
# Plots and visualization related to the Boolean Model

from collections import defaultdict

import matplotlib
import matplotlib.pylab as plt
import numpy as np
from SecretColors import Palette, ColorMap
from matplotlib.patches import Patch

from analysis.boolean import is_expressed
from constants.boolean import INTERESTED_GENES


def boolean_expression_pattern(genes: list, exp_hours: list):
    p = Palette()
    on_color = p.blue()
    off_color = p.red()
    cmap = ColorMap(matplotlib, p).from_list([off_color, on_color],
                                             is_qualitative=True)
    exp = defaultdict(list)
    for g in genes:
        for h in exp_hours:
            exp[h].append(is_expressed(g, h))

    full = []
    for h in exp_hours:
        full.append(exp[h])
    fig, ax = plt.subplots()
    full = np.asarray(full).astype(int)
    plt.pcolormesh(full.T, edgecolors=p.blue(shade=10), linewidth=0.1,
                   cmap=cmap)
    ax.set_aspect('0.4')
    ax.set_xticks(np.arange(0, len(exp_hours)) + 0.5)
    ax.set_yticks(np.arange(0, len(genes)) + 0.5)
    ax.set_xticklabels(["{} hpf".format(x) for x in exp_hours])
    ax.set_yticklabels([str(x).capitalize() for x in genes])
    ax.xaxis.set_ticks_position('top')
    # ax.set_title("RNA-seq Boolean Conditions", y=1.08)
    legend_elements = [Patch(facecolor=on_color, label='ON'),
                       Patch(facecolor=off_color, label='OFF')]
    ax.legend(handles=legend_elements, loc='upper center',
              bbox_to_anchor=(0.5, -0.05), fancybox=True, ncol=2)
    plt.tight_layout()
    plt.show()


def run():
    boolean_expression_pattern(INTERESTED_GENES, [20, 24, 48, 72])
