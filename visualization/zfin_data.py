"""
CardioTRaNs 2019
Author: Rohit Suratekar

All functions related to visualization of ZFIN data-sets
"""

from collections import Counter

import matplotlib.pylab as plt
from SecretColors.palette import Palette
from matplotlib_venn import venn2_unweighted

from analysis.basic import *
from helpers.parsers import *


def plot_expression_factors_wise(no_of_factors: int = 20):
    data = get_expression_data_frame()
    fc = Counter()
    p = Palette()
    for d in data[COL_EXP_1_GENE_SYMBOL]:
        fc.update({d})

    names = []
    values = []
    for t in fc.most_common(no_of_factors):
        names.append(t[0])
        values.append(t[1])

    ind = range(len(names))

    plt.barh(ind, values, color=p.violet())
    plt.yticks(ind, names)
    plt.xlabel("Number of expression data available in ZFIN")
    plt.ylabel("Gene Symbol")
    plt.title("Top {} factors with most expression data in ZFIN".format(
        no_of_factors))
    plt.show()


def plot_expression_structure_wise(no_of_structures):
    data = get_expression_data_frame()
    fc = Counter()
    p = Palette()
    for d in data[COL_EXP_4_SUPER_STR_NAME]:
        fc.update({d})

    names = []
    values = []
    for t in fc.most_common(no_of_structures):
        names.append(t[0])
        values.append(t[1])

    ind = range(len(names))

    plt.barh(ind, values, color=p.pink())
    plt.yticks(ind, names)
    plt.xlabel("Number of expression data available in ZFIN")
    plt.ylabel("Super Structure Names")
    plt.title("Top {} super structures with most expression data in "
              "ZFIN".format(no_of_structures))
    plt.show()


def plot_expression_stage_wise():
    data = get_expression_data_frame()
    stage = parse_stage_ontology()
    stage = {x.name: x.start for x in stage}
    fc = Counter()
    p = Palette()
    for d in data[COL_EXP_7_START_STAGE]:
        fc.update({d})

    names = []
    values = []
    time = []

    for t in fc.most_common():
        names.append(t[0])
        values.append(t[1])
        time.append(stage[t[0]])

    ind = range(len(names))

    all_zipped = zip(time, names, values)
    all_zipped = sorted(all_zipped)
    all_zipped = zip(*all_zipped)
    time, names, values = [list(x) for x in all_zipped]

    colors = [p.blue(shade=50)] * len(names)
    high_light_color = p.red(shade=40)
    high_light_items = [16.0, 24.0, 48.0]
    for a in high_light_items:
        colors[time.index(a)] = high_light_color

    for i in range(len(names)):
        names[i] = "{} ({})".format(names[i], time[i])
    plt.barh(ind, values, color=colors)
    plt.yticks(ind, names)

    for a in high_light_items:
        plt.gca().get_yticklabels()[time.index(a)].set_color(p.red(shade=50))

    plt.xlabel("Number of expression data available in ZFIN")
    plt.ylabel("Developmental Start Stage (hpf)")
    plt.title("Stage wise expression data in ZFIN")
    plt.show()


def plot_with_filters(structures: list, star_time: float, end_time: float,
                      no_of_genes: int = 20):
    # Get all data
    data = get_expression_data_frame()
    # Filter according to hours
    filtered_data = get_genes_for_hours(data, start_hour=star_time,
                                        end_hour=end_time)
    filtered_data = get_genes_from_structure(filtered_data, structures)

    c = Counter(list(filtered_data[COL_EXP_1_GENE_SYMBOL]))

    names = []
    values = []

    for k in c.most_common(no_of_genes):
        names.append(k[0])
        values.append(k[1])

    ind = range(len(names))
    plt.barh(ind, values, color=Palette().red(shade=40))
    plt.yticks(ind, names)
    plt.xlabel("Number of expression data available in ZFIN")
    plt.ylabel("Gene Symbol")
    # plt.title("Top {} genes between {} - {} hpf\n(structures "
    #           "with filtering terms {})".format(no_of_genes,
    #                                             star_time, end_time,
    #                                             structures))
    plt.show()


def plot_zfin_gr_venn():
    p = Palette()
    zfin_data = get_expression_data_frame()
    zfin_data = get_genes_from_structure(zfin_data, ["heart", "cardi"])
    zfin_data = get_genes_for_hours(zfin_data, 48, 48)
    gr_data = get_all_gr_expressed_data()
    set1 = set(gr_data[COL_GR_S1_GENE_SYMBOL].values)
    set2 = set(zfin_data[COL_EXP_1_GENE_SYMBOL].values)

    c = venn2_unweighted([set1, set2], ('DESeq2', 'ZFIN'))
    c.get_patch_by_id('10').set_color(p.cerulean())
    c.get_patch_by_id('11').set_color(p.red())
    c.get_patch_by_id('01').set_color(p.blue())
    plt.show()


def run():
    plot_with_filters(structures=["heart", "cardi"], star_time=24, end_time=24)
