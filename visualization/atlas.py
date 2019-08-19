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
from SecretColors import Palette, ColorMap
from matplotlib.patches import Patch

from constants.boolean import *
from constants.other import *
from constants.system import PLOT_FOLDER
from constants.zfin import *
from helpers.parsers.other import get_expression_atlas
from helpers.parsers.zfin import get_stage_ontology


def check_expression_pattern(genes: list, tpm_cutoff: float = 1):
    """
    Checks expression pattern in the expression atlas data
    :param genes: List of genes to check the expression
    :param tpm_cutoff: TPM Cutoff used to keep genes off or on
    """
    p = Palette()
    on_color = p.peach(shade=40)
    off_color = p.violet(shade=70)
    cmap = ColorMap(matplotlib, p).from_list([off_color, on_color],
                                             is_qualitative=True)

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
    data = data[data[COL_EXP_ATLAS_2_GENE_NAME].isin(genes)]
    del data[COL_EXP_ATLAS_1_GENE_ID]
    data = data.set_index(COL_EXP_ATLAS_2_GENE_NAME)
    data = data.fillna(0)
    data[data < tpm_cutoff] = 0
    data[data >= tpm_cutoff] = 1
    fig, ax = plt.subplots()
    plt.pcolormesh(data.values, edgecolor=p.gray(shade=70),
                   linewidth=0.05, cmap=cmap)
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


def run():
    check_expression_pattern(CONTROL_GENES, 2)
