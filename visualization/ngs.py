import matplotlib
import matplotlib.pylab as plt
import numpy as np
from SecretColors.palette import Palette, ColorMap
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from scipy.stats import gaussian_kde

from analysis.ngs import *
from constants.boolean import *
from helpers.parsers.ngs import get_rna_seq_data


def visualize_single_chromosome(hour: int, chromosome_no, genes: list,
                                max_size: int = 100):
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

    plt.hist(data[COL_STRING_TIE_5_START], 200, color=p.blue(shade=30))
    plt.ylabel("No of genes")
    plt.xticks(x_ticks,
               ["{:.2f} MB".format(x * max_pos / (100 * 1000000)) for x in
                x_ticks])
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel("Chromosome No: {} at {} hpf".format(chromosome_no, hour))
    plt.show()


def map_all_chromosomes_bar(hour: int, min_exp_percentage=0.01):
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
    plt.title("TPM greater than {} % of largest per chromosome at {} "
              "hpf".format(min_exp_percentage, hour))
    plt.show()


def map_given_genes(hour: int, genes: list):
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


def plot_fpkm_tpm_density(hour: int):
    """
    Plots the FPKM and TPM density distributions
    """
    p = Palette()
    data = get_rna_seq_data(hour)
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


def plot_expression_series(genes: list):
    p = Palette(allow_gray_shades=True)
    hours = [24, 48, 72]
    h_data = [get_rna_seq_data(x) for x in hours]
    data = {}

    for g in genes:
        k = []
        for h in h_data:
            k.append(filter_with_gene(h, g)[COL_STRING_TIE_9_TPM].values[0])
        data[g] = k

    ind = range(len(hours))
    colors = p.random(no_of_colors=len(genes), shade=50)
    for i, g in enumerate(genes):
        plt.plot(ind, data[g], label=g, color=colors[i], marker='o',
                 linestyle="--")
    plt.xticks(ind, ["{} hpf".format(x) for x in hours])
    plt.ylabel("Average TPM")
    plt.legend(loc=0)
    plt.show()


def run():
    # visualize_single_chromosome(24, 14, INTERESTED_GENES)
    plot_expression_series(BASE_GENES)
    # plot_single_gene_expression("nkx2.5")
