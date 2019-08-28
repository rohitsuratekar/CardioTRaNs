#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 19/08/19, 4:17 PM
#
#  Copyright (c) 2019.
#
# Data obtained from the EBI Expression atlas
# https://www.ebi.ac.uk/gxa/experiments/E-ERAD-475/Results

import matplotlib
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from SecretColors import Palette, ColorMap
from matplotlib.patches import Patch

from constants.other import *
from constants.system import PLOT_FOLDER
from constants.zfin import *
from helpers.parsers.other import get_expression_atlas
from helpers.parsers.zfin import get_stage_ontology


def _normalized_expression_atlas():
    stages = get_stage_ontology()
    data = get_expression_atlas()
    stages[COL_ONT_3_STAGE_NAME] = stages[COL_ONT_3_STAGE_NAME].str.replace(
        ":", " ")
    stages[COL_ONT_3_STAGE_NAME] = stages[COL_ONT_3_STAGE_NAME].str.lower()

    new_names = {}
    for d in data.columns:
        t = stages[stages[COL_ONT_3_STAGE_NAME] == d]
        if len(t) > 0:
            new_names[d] = t[COL_ONT_4_BEGIN_HOUR].values[0]
        else:
            if d == "zygote":
                new_names[d] = 0.0
            elif d == "larval protruding mouth":
                new_names[d] = 72.0
            else:
                new_names[d] = d

    data = data.rename(columns=new_names)
    return data


def _plot_dataframe(data: pd.DataFrame, genes: list, tpm_cutoff: float):
    p = Palette()
    on_color = p.peach(shade=40)
    off_color = p.violet(shade=70)
    cmap = ColorMap(matplotlib, p).from_list([off_color, on_color],
                                             is_qualitative=True)

    data = data.fillna(0)
    data[data < tpm_cutoff] = 0
    data[data >= tpm_cutoff] = 1

    fig, ax = plt.subplots()
    plt.pcolormesh(data.values, edgecolor=p.gray(shade=70),
                   linewidth=0.05, cmap=cmap, vmin=0, vmax=1)
    ax.set_aspect(1.2)
    ax.set_xticks(np.arange(0, len(data.columns)) + 0.5)
    ax.set_yticks(np.arange(0, len(genes)) + 0.5)
    ax.set_xticklabels(["{}".format(x) for x in data.columns], rotation=45,
                       ha="left")
    ax.set_yticklabels([str(x).capitalize() for x in genes])
    plt.xlabel("Hours Post Fertilization (whole animal)")
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    legend_elements = [Patch(facecolor=on_color, label='ON'),
                       Patch(facecolor=off_color, label='OFF')]
    ax.legend(handles=legend_elements, loc='upper center',
              bbox_to_anchor=(0.5, -0.05), fancybox=True, ncol=2,
              title="with TPM > {} ".format(tpm_cutoff))

    plt.tight_layout()

    plt.savefig(PLOT_FOLDER + "expression_atlas.png", dpi=300, type="png")
    plt.show()


def check_expression_pattern(genes: list, tpm_cutoff: float = 1):
    """
    Checks expression pattern in the expression atlas data
    :param genes: List of genes to check the expression
    :param tpm_cutoff: TPM Cutoff used to keep genes off or on
    """

    data = _normalized_expression_atlas()
    data = data[data[COL_EXP_ATLAS_2_GENE_NAME].isin(genes)]
    del data[COL_EXP_ATLAS_1_GENE_ID]
    data = data.set_index(COL_EXP_ATLAS_2_GENE_NAME)
    _plot_dataframe(data, genes, tpm_cutoff)


def include_expression(data,
                       include_from: float,
                       include_till: float,
                       ignore_others: bool = False,
                       exclude_mode: bool = False,
                       tpm_cutoff: float = 1):
    try:
        del data[COL_EXP_ATLAS_1_GENE_ID]
    except KeyError:
        # Ignore the post processed data
        pass
    data = data.set_index(COL_EXP_ATLAS_2_GENE_NAME)
    # Check if all time points are available
    times = []
    for t in data.columns.values:
        if include_from <= t <= include_till:
            times.append(t)

    data = data.fillna(0)

    for t in data.columns.values:
        if t in times:
            if exclude_mode:
                data = data[data[t] < tpm_cutoff]
            else:
                data = data[data[t] >= tpm_cutoff]
        else:
            if not ignore_others:
                if exclude_mode:
                    data = data[data[t] >= tpm_cutoff]
                else:
                    data = data[data[t] < tpm_cutoff]

    temp = "temp"
    data[temp] = data.sum(axis=1)

    data = data.sort_values(by=temp, ascending=False)
    del data[temp]
    return data.reset_index(drop=False)


def exclude_expression(data, exclude_from: float, exclude_till: float,
                       ignore_others: bool = True,
                       tpm_cutoff: float = 1):
    return include_expression(data, include_from=exclude_from,
                              include_till=exclude_till,
                              exclude_mode=True,
                              ignore_others=ignore_others,
                              tpm_cutoff=tpm_cutoff)


def run():
    data = _normalized_expression_atlas()
    # data = exclude_expression(data, exclude_from=20 , exclude_till=80)
    data = include_expression(data, include_from=80, include_till=130)
    data = data.set_index(COL_EXP_ATLAS_2_GENE_NAME)
    _plot_dataframe(data.head(15), list(data.head(15).index), 1)
