"""
CardioTRaNs 2019
Author: Rohit Suratekar

All functions related to visualization of NGS data-sets
"""
import matplotlib
import matplotlib.pylab as plt
import numpy as np
from SecretColors.palette import Palette, ColorMap
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

from analysis.ngs_analysis import *
from constants.boolean_constants import INTERESTED_GENES


def visualize_single_chromosome(chromosome_no, genes: list,
                                max_size: int = 100):
    chrome_y_loc = 0
    chrome_height = 0.3
    gene_onset = 0.3
    p = Palette()

    # Prepare data with filtering
    data = get_rna_seq_gene_data()
    data = filter_with_chromosome(data, chromosome_no)
    # Store maximum location
    max_pos = data[COL_STRING_TIE_6_END].max()
    # Normalize data
    data[COL_STRING_TIE_6_END] = data[
                                     COL_STRING_TIE_6_END] * max_size / max_pos
    data[COL_STRING_TIE_5_START] = data[
                                       COL_STRING_TIE_5_START] * max_size / max_pos

    # Prepare gene locations (if any)
    gene_list = []
    if len(genes) == 0:
        genes = [x for x in INTERESTED_GENES]

    for g in genes:
        d = filter_with_gene(data, g)
        if len(d) == 1:
            gene_list.append(
                [d[COL_STRING_TIE_2_GENE_NAME].values[0],
                 d[COL_STRING_TIE_5_START].values[0]])

    # Prepare colormap
    c_colors = [p.black(), p.lime(), p.lime(shade=30)]
    c = ColorMap(matplotlib, p).from_list(c_colors)

    y_loc = (chrome_y_loc, chrome_height)

    # Plot
    plt.figure()
    ax1 = plt.subplot(211, frameon=False)

    plt.yticks([])
    plt.xticks([])

    # plt.broken_barh([(0, data[COL_STRING_TIE_6_END].max())], y_loc,
    #                 color=p.white())
    for i, row in data.iterrows():
        start = row[COL_STRING_TIE_5_START]
        end = row[COL_STRING_TIE_6_END]
        plt.broken_barh([(start, end - start)], y_loc,
                        color=c(row[COL_STRING_TIE_9_TPM]))

    x_ticks = np.linspace(0, 100, 10)
    plt.ylim(0, chrome_y_loc + chrome_height + gene_onset + 0.2)

    ax_divider = make_axes_locatable(ax1)
    # define size and padding of axes for colorbar
    cax = ax_divider.append_axes('top', size='8%')

    norm = matplotlib.colors.Normalize(vmin=0, vmax=max_pos)
    sm = plt.cm.ScalarMappable(cmap=c, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ticks=np.linspace(0, max_pos, 10),
                 orientation="horizontal", cax=cax)
    plt.title("TPM normalized to highest on given chromosome")

    ax2 = plt.subplot(212, frameon=False)
    for g in gene_list:
        ax1.text(g[1], chrome_y_loc + chrome_height + gene_onset + 0.01, g[0],
                 ha='center')
        ax1.plot([g[1], g[1]], [chrome_y_loc + chrome_height,
                                chrome_y_loc + chrome_height + gene_onset],
                 linestyle="--")
        ax2.axvline(g[1], linestyle="--")

    plt.hist(data[COL_STRING_TIE_5_START], 200, color=p.blue(shade=30))
    plt.ylabel("No of genes")
    plt.xticks(x_ticks,
               ["{:.2f} MB".format(x * max_pos / (100 * 1000000)) for x in
                x_ticks])
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel("Chromosome No: {}".format(chromosome_no))
    plt.show()


def map_all_chromosomes_bar(min_exp_percentage=0.01):
    p = Palette()
    data = get_rna_seq_gene_data()
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


def map_given_genes(genes: list):
    p = Palette()
    data = get_rna_seq_gene_data()
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
            if row[COL_STRING_TIE_2_GENE_NAME] in genes:
                plt.broken_barh([(start, end - start)], y_loc,
                                color=c(row[COL_STRING_TIE_9_TPM]))

    plt.yticks([x + 0.4 for x in range(len(y_ticks))],
               ["chr {}".format(x) for x in y_ticks])
    x_ticks = np.linspace(0, 100, 10)
    plt.xticks(x_ticks,
               ["{:.2f} MB".format(x * max_pos / (100 * 1000000)) for x in
                x_ticks])
    plt.xlabel("Position")
    plt.title("Gene Mapping")
    plt.show()


def run():
    map_given_genes(INTERESTED_GENES)
