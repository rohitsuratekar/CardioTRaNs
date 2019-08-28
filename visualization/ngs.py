#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 29/07/19, 10:00 AM
#
#  Copyright (c) 2019.
#
# All functions related to NGS analysis will go here

import matplotlib
import matplotlib.pylab as plt
import numpy as np
from SecretColors import Palette, ColorMap
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from scipy.stats import gaussian_kde

from analysis.ngs import *
from constants.boolean import *
from constants.system import GENOTYPE_WT
from helpers.parsers.ngs import get_rna_seq_data, get_mapping_property


def visualize_single_chromosome(hour: int, chromosome_no, genes: list,
                                tpm_cutoff=0, max_size: int = 100):
    """
    Visualize single chromosome by showing gene locations on that chromosome
    :param hour: Hour post fertilization
    :param chromosome_no: Chromosome number
    :param genes: Genes to visualize (if any)
    :param tpm_cutoff: This cut off will be used to plot number of genes
    histogram at the bottom
    :param max_size: Maximum width in arbitrary units
    """

    chrome_y_loc = 0
    chrome_height = 0.3
    gene_onset = 0.3
    p = Palette()

    # Prepare data with filtering
    data = get_rna_seq_data(hour)
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
                 linestyle="--", color=p.gray())
        ax2.axvline(g[1], linestyle="--", color=p.gray())

    data = data[data[COL_STRING_TIE_9_TPM] > tpm_cutoff]
    plt.hist(data[COL_STRING_TIE_5_START], 200, color=p.blue(shade=40))
    plt.ylabel("No of genes (with TPM > {})".format(tpm_cutoff))
    plt.xticks(x_ticks,
               ["{:.2f} MB".format(x * max_pos / (100 * 1000000)) for x in
                x_ticks])
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel("Chromosome No: {} at {} hpf".format(chromosome_no, hour))
    plt.show()


def map_all_chromosomes_bar(hour: int, tpm_percentage_cutoff=0.01):
    """
    Plots all 25 chromosomes plus mitochondrial chromosomes

    :param hour: Hours post fertilization
    :param tpm_percentage_cutoff: Only genes greater than this cutoff will be
    considered for plotting. Each cutoff is calculated individually for
    every chromosome.
    """
    p = Palette()
    data = get_rna_seq_data(hour)
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
            if row[COL_STRING_TIE_9_TPM] > tpm_percentage_cutoff:
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
    plt.title("TPM greater than {} % of largest per chromosome at {} "
              "hpf".format(tpm_percentage_cutoff, hour))
    plt.show()


def map_given_genes(hour: int, genes: list, genotype: str):
    """
    Shows location of given gene on the chromosome plot
    :param hour: Hours post fertilization
    :param genes: List of genes to map
    :param genotype: Genotype
    """
    p = Palette()
    data = get_rna_seq_data(hour, genotype)
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


def plot_fpkm_tpm_density(hour: int, genotype: str):
    """
    Plots the FPKM and TPM density distributions
    """
    p = Palette()
    data = get_rna_seq_data(hour, genotype)
    d = data[data[COL_STRING_TIE_8_FPKM] > 0][COL_STRING_TIE_8_FPKM].apply(
        lambda x: np.log2(x))

    d2 = data[data[COL_STRING_TIE_9_TPM] > 0][COL_STRING_TIE_9_TPM].apply(
        lambda x: np.log2(x))

    density1 = gaussian_kde(d)
    density2 = gaussian_kde(d2)

    _, r1, _ = plt.hist(d.values, 50, color=p.blue(shade=50), density=True)
    _, r2, _ = plt.hist(d2.values, 50, color=p.red(shade=50), alpha=0.8,
                        density=True)
    plt.axvline(0, linestyle="--", color=p.black())
    plt.plot(r1, density1(r1), linewidth=3, color=p.blue(shade=60),
             label="FPKM")
    plt.plot(r2, density2(r2), linewidth=3, color=p.red(shade=60),
             label="TPM")
    plt.legend(loc=0)
    plt.xlabel("$Log_2$(Entity)")
    plt.ylabel("Normalized Density")
    plt.title("Average for {} hpf".format(hour))
    plt.show()


def plot_single_gene_expression(gene: str):
    """
    Shows expression parameters for given gene
    :param gene: Name of genes
    """
    p = Palette()
    hours = [24, 48, 72]
    data = [get_rna_seq_data(x) for x in hours]
    tpm = [filter_with_gene(x, gene)[COL_STRING_TIE_9_TPM].values[0] for x in
           data]
    fpkm = [filter_with_gene(x, gene)[COL_STRING_TIE_8_FPKM].values[0] for x in
            data]
    coverage = [filter_with_gene(x, gene)[COL_STRING_TIE_7_COVERAGE].values[0]
                for x in data]

    ind = np.asarray(range(len(hours)))
    w = 0.25
    plt.bar(ind, tpm, width=w, color=p.blue(shade=40), label="TPM")
    plt.bar(ind + w, fpkm, width=w, color=p.blue(shade=60), label="FPKM")
    plt.bar(ind + 2 * w, coverage, width=w, color=p.blue(shade=80),
            label="Coverage")
    plt.xticks(ind + w, ["{} hpf".format(x) for x in hours])
    plt.ylabel("Average Quantity")
    plt.legend(loc=0)
    plt.title("{} expression".format(gene))
    plt.show()


def plot_expression_series(genes: list, genotype: str):
    """
    Show how given gene changes their TPM values over the developmental time
    point
    :param genes: List of genes
    :param genotype: Genotype of the sample
    """
    p = Palette()
    hours = [24, 48, 72]
    h_data = [get_rna_seq_data(x, genotype) for x in hours]
    data = {}

    for g in genes:
        k = []
        for h in h_data:
            k.append(filter_with_gene(h, g)[COL_STRING_TIE_9_TPM].values[0])
        data[g] = k

    ind = range(len(hours))
    # colors = p.random(no_of_colors=len(genes), shade=50)
    colors = [p.blue(shade=60), p.red(shade=50), p.yellow(shade=30),
              p.orange(shade=40), p.aqua()]
    for i, g in enumerate(genes):
        plt.plot(ind, data[g], label=str(g).capitalize(), color=colors[i],
                 marker='o',
                 linestyle="--")
    plt.xticks(ind, ["{} hpf".format(x) for x in hours])
    plt.grid(True, axis="y", color=p.gray(shade=25), linestyle='-')
    plt.gca().patch.set_facecolor(p.gray(shade=10))
    plt.ylabel("Average TPM")
    plt.legend(loc=0)
    plt.show()


def chromosome_activity(chromosome_number: int, tpm_cutoff=0, divisions=100):
    """
    Plots gene expression activity of given chromosome in the form of TPM in
    given bins

    :param chromosome_number: Chromosome number
    :param tpm_cutoff: All genes below this TPM cut off will be discarded
    :param divisions: Number of divisions in which chromosome will be
    divided or binned to check total TPM values
    """
    p = Palette()
    hours = [24, 48, 72]
    colors = [p.cyan(), p.yellow(shade=40), p.magenta()]

    dv = np.linspace(0, 100, divisions)
    temp_name = "temp"
    max_pos = 0

    def _format_column(n):
        return min(dv, key=lambda x: abs(x - n))

    def _format_data(h):
        genes1 = {}
        d = get_rna_seq_data(h)
        d = filter_with_chromosome(d, chromosome_number)

        for g in INTERESTED_GENES:
            f = filter_with_gene(d, g)
            if len(f) == 1:
                genes1[g] = f[COL_STRING_TIE_5_START].values[0]

        d = d[d[COL_STRING_TIE_9_TPM] > tpm_cutoff]
        max_loc = d[COL_STRING_TIE_6_END].max()
        d[COL_STRING_TIE_5_START] = d[COL_STRING_TIE_5_START] * 100 / max_loc
        d[temp_name] = d[COL_STRING_TIE_5_START].apply(_format_column)
        d = d.groupby(by=temp_name).sum()
        d = d.reset_index()
        return d, max_loc, genes1

    genes = None
    for i, hr in enumerate(hours):
        data, max_pos, genes = _format_data(hr)
        plt.plot(data[temp_name], data[COL_STRING_TIE_9_TPM], marker="o",
                 label="{} hpf".format(hr), color=colors[i])

    if genes is not None:
        for k in genes:
            loc = genes[k] * 100 / max_pos
            plt.axvline(loc, linestyle="--", color=p.gray())

            plt.text((genes[k] / max_pos), 0.9, k,
                     horizontalalignment='center',
                     verticalalignment='center',
                     bbox=dict(facecolor=p.gray(shade=10)),
                     transform=plt.gca().transAxes)

    x_ticks = np.linspace(0, 100, 10)
    plt.xticks(x_ticks,
               ["{:.2f} MB".format(x * max_pos / (100 * 1000000)) for x in
                x_ticks])
    plt.legend(loc=0)
    plt.ylabel("TPM")
    plt.title("Chromosome {} (Divisions = {}, TPM > {})".format(
        chromosome_number, divisions, tpm_cutoff))
    plt.show()


def mapping_stats(name: str):
    """
    Creates box plot for given mapping property
    :param name: EXACT name of the property
    """
    p = Palette()
    data = get_mapping_property(name)
    data = [float(x.replace("%", "")) for x in data]
    plt.boxplot(data, patch_artist=True,
                medianprops=dict(color=p.black()),
                boxprops=dict(facecolor=p.green(), color=p.green()))
    plt.xlabel(name)
    plt.xticks([])
    plt.show()


def run():
    # mapping_stats("Uniquely mapped reads % ")
    # visualize_single_chromosome(24, 14, INTERESTED_GENES, tpm_cutoff=1)
    # plot_expression_series(BASE_GENES, GENOTYPE_WT)
    plot_fpkm_tpm_density(72, GENOTYPE_WT)
    # map_given_genes(24, BASE_GENES, GENOTYPE_WT)
    # chromosome_activity(1)
    # test()
