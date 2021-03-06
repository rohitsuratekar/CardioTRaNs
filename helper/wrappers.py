#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#

import itertools

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SecretColors import Palette
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA

from analysis.deseq2 import ma_plot, volcano_plot
from constants.outputs import *
from helper.resolver import NameResolver

palette = Palette(show_warning=False)


def single_ma(nr: NameResolver):
    method = "star"
    lab = "hills"
    condition = [30, 48]

    _, ax = plt.subplots()
    filename = nr.deseq2_results(method,
                                 lab=lab,
                                 condition=condition,
                                 lfc=True)

    ma_plot(filename,
            color=palette.cerulean(shade=60),
            accent_color=palette.green_light(),
            ax=ax,
            marker=".")

    ax.set_xlim(0, 1000)
    ax.set_ylim(-15, 15)
    ax.set_title(f"{condition[0]} v {condition[1]} ({method}, {lab} lab)")
    plt.tight_layout()
    plt.savefig("ma.png", dpi=300)
    plt.show()


def count_pca(nr: NameResolver):
    all_conditions = sorted([30, 48, 72])
    lab = "hills"
    methods = ["stringtie"]
    combinations = list(itertools.combinations(all_conditions, 2))
    dfs = []
    r = {}
    for pair in combinations:
        for method in methods:
            f = nr.deseq2_vst(method, lab=lab, condition=list(pair))
            df = pd.read_csv(f).set_index("gene_id")

            col_change = {}
            for c in df.columns:
                col_change[c] = f"{c}_{method}"
            df = df.rename(columns=col_change)
            cols = df.columns
            r[f"{pair[0]}_{method}"] = list(cols[: int(len(cols) / 2)])
            r[f"{pair[1]}_{method}"] = list(cols[int(len(cols) / 2):])
            dfs.append(df)

    colors = []
    columns = []
    labels = {}
    time_colors = {}
    for i, t in enumerate(all_conditions):
        time_colors[str(t)] = palette.get_color_list[i]

    for color, key in zip(palette.get_color_list[:len(r)], r):
        labels[key] = color
        columns.extend(r[key])
        # colors.extend([color] * len(r[key]))
        colors.extend([time_colors[key.split("_")[0]]] * len(r[key]))

    df = pd.concat(dfs, axis=1, join="outer", sort=False)
    df = df.loc[:, ~df.columns.duplicated()]
    df = df[columns].fillna(0)

    pca = PCA()
    a1 = pca.fit(df.values)
    var = a1.components_
    percentages = [round(x, 2) for x in a1.explained_variance_ratio_ * 100]
    plt.scatter(var[0, :], var[1, :], color=colors)
    plt.xlabel(f"PC1: {percentages[0]} %")
    plt.ylabel(f"PC2: {percentages[1]} %")

    hdl = []
    for key in time_colors:
        hdl.append(Patch(color=time_colors[key],
                         label=f"{key} hpf"))

    plt.legend(handles=hdl, loc=0)
    plt.title(f"Counts with {methods} ({lab} lab)")
    plt.tight_layout()
    plt.savefig("pca.png", dpi=300)
    plt.show()


def get_fold_change_distribution(nr: NameResolver):
    lab = "winata"
    methods = ["star", "salmon", "kallisto"]
    condition = [24, 72]
    colors = palette.get_color_list
    p_value = 0.05
    for m, c in zip(methods, colors):
        fn = nr.deseq2_results(method=m, lab=lab, condition=condition,
                               lfc=True)
        df = pd.read_csv(fn)
        data = df
        data = data[data[DESEQ2_PADJ] < p_value]
        data = data[DESEQ2_LOG2_CHANGE].values
        density = gaussian_kde(data)
        xs = np.linspace(-20, 20, 1000)
        plt.plot(xs, density(xs), color=c, label=m, lw=2)
        plt.fill_between(xs, density(xs), color=c, alpha=0.1)
    plt.axvline(0, ls="--", color=palette.black())
    plt.ylabel(f"Probability Density (p-value<{p_value})")
    plt.xlabel("Log$_2$ Fold Change")
    plt.legend(loc=0)
    plt.title(
        f"{condition[0]} v {condition[1]} ({lab} lab)")
    plt.xlim(-20, 20)
    plt.grid(axis="both", zorder=0, ls=":")
    plt.ylim(0)
    plt.tight_layout()
    plt.savefig("density.png", dpi=300)
    plt.show()


def mean_sd_plot(nr: NameResolver, *, vst: bool):
    lab = "winata"
    condition = [24, 48]
    method = "kallisto"

    fn = nr.deseq2_counts(method=method, lab=lab, condition=condition)
    if vst:
        fn = nr.deseq2_vst(method=method, lab=lab, condition=condition)

    data = pd.read_csv(fn).set_index("gene_id")
    data["temp"] = data.sum(axis=1)
    plt.scatter(data.mean(axis=1), data.std(axis=1), marker=".",
                c=data["temp"])
    plt.xlabel("Mean")
    plt.ylabel("Standard Deviation")
    clb = plt.colorbar()
    clb.ax.set_title("Total Counts")
    plt.title(f"{condition[0]} vs {condition[1]} ({method}, {lab} lab)")
    plt.tight_layout()
    plt.savefig("mean_sd.png", dpi=300)
    plt.show()


def plot_volcano(nr: NameResolver):
    lab = "winata"
    method = "salmon"
    condition = [24, 48]
    fn = nr.deseq2_results(method=method, lab=lab,
                           condition=condition,
                           lfc=True)
    volcano_plot(fn)
    plt.ylim(0, 100)
    plt.xlim(-15, 15)
    plt.axvline(0, ls="--", color=palette.black())
    plt.title(f"{condition[0]} vs {condition[1]} ({method}, {lab} lab)")
    plt.tight_layout()
    plt.savefig("volcano.png", dpi=300)
    plt.show()


def plot_all_volcano(nr: NameResolver):
    methods = ["salmon", "kallisto", "stringtie"]
    conditions = [[30, 48], [48, 72], [30, 72]]
    max_cols = 3
    lab = "hills"

    txt_cond = []
    data = []
    for m in methods:
        txt_cond.extend([f"{x[0]}v{x[1]} {m}" for x in conditions])
        for con in conditions:
            fn = nr.deseq2_results(m, lab=lab, condition=con, lfc=True)
            data.append(fn)

    condition_array = txt_cond
    dfs = data

    rows = int(np.ceil(len(condition_array) / max_cols))
    offset = max_cols * rows - len(condition_array)
    temp = np.zeros(len(condition_array) + offset)
    if max_cols > len(condition_array):
        max_cols = len(condition_array)
        temp = np.zeros(len(condition_array))

    temp = temp.reshape(-1, max_cols)
    gs = gridspec.GridSpec(temp.shape[0], temp.shape[1])
    fig = plt.figure()
    for ind, condition in enumerate(condition_array):
        ax = fig.add_subplot(gs[ind])  # type:plt.Axes
        volcano_plot(dfs[ind], ax=ax,
                     add_labels=False)
        ax.set_title(condition)
        ax.set_xlim(-10, 10)
        ax.set_ylim(0, 50)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none',
                    top=False,
                    bottom=False,
                    left=False,
                    right=False)
    plt.ylabel("-Log$_{10}$ Adj P value")
    plt.xlabel("Log$_2$ Fold change")
    plt.tight_layout()
    plt.savefig("all_volcano.png", dpi=300)
    plt.show()


def run():
    nr = NameResolver("config.json")
    plot_volcano(nr)
