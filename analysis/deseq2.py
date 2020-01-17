#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
# All analysis related to DESeq2

import matplotlib.gridspec as gridspec
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from SecretColors import Palette
from matplotlib.lines import Line2D

from constants.methods import *
from helper import ConfigParser, FileManager


def get_data(cm: ConfigParser, method, conditions, *, lab) -> list:
    fm = FileManager(cm)
    dfs = []
    for ar in conditions:
        res = fm.deseq2_results(method,
                                lab=lab,
                                conditions=ar,
                                lfc=True)
        dfs.append(res)
    return dfs


def ma_plot(df: pd.DataFrame):
    p = Palette()
    data = df.fillna(0)
    data["color"] = p.blue()
    data.loc[data[DESEQ2_PADJ] < 0.05, "color"] = p.red()
    plt.scatter(data[DESEQ2_BASE_MEAN],
                data[DESEQ2_LOG2_CHANGE],
                color=data["color"])
    plt.axhline(0, ls="--", color="k")
    plt.xlim(0, 5000)
    plt.ylabel("Log$_2$ Fold change")
    plt.xlabel("Average normalize counts")
    plt.show()


def volcano_plot(df: pd.DataFrame,
                 p_value: float = 0.05,
                 color=Palette().blue(),
                 significant=Palette().red(),
                 ax: plt.Axes = None,
                 add_legend: bool = True,
                 add_labels: bool = True,
                 show_plot: bool = True):
    p = Palette()
    if ax is None:
        _, ax = plt.subplots()
    data = df.fillna(0)
    data["color"] = color
    data.loc[data[DESEQ2_PADJ] < p_value, "color"] = significant

    data["temp"] = data[DESEQ2_PADJ].map(lambda x: -np.log10(x))
    ax.scatter(data[DESEQ2_LOG2_CHANGE],
               data["temp"], color=data["color"], marker=".")

    ax.axvline(-1, ls="--", color=p.gray(), label="2 Fold change")
    ax.axvline(1, ls="--", color=p.gray())
    if add_legend:
        le = [
            Line2D([], [], marker=".",
                   color=significant,
                   label=f"Adj. P value < {p_value}",
                   linestyle="None"),
            Line2D([], [], linestyle="--", color=p.gray(),
                   label="2 Fold change")
        ]
        ax.legend(handles=le, loc=0)

    if add_labels:
        ax.set_ylabel("-Log$_{10}$ Adj P value")
        ax.set_xlabel("Log$_2$ Fold change")

    ax.set_xlim(-10, 10)

    if show_plot:
        plt.show()


def plot_all_volcano(condition_array: list, dfs: list, max_cols=3):
    rows = int(np.ceil(len(condition_array) / max_cols))
    offset = max_cols * rows - len(condition_array)
    temp = np.zeros(len(condition_array) + offset)
    if max_cols > len(condition_array):
        max_cols = len(condition_array)
    temp = temp.reshape(-1, max_cols)
    gs = gridspec.GridSpec(temp.shape[0], temp.shape[1])
    fig = plt.figure()
    for ind, condition in enumerate(condition_array):
        ax = fig.add_subplot(gs[ind])  # type:plt.Axes
        volcano_plot(dfs[ind], ax=ax,
                     show_plot=False,
                     add_legend=False,
                     add_labels=False)
        ax.set_title(condition)

    plt.tight_layout()
    plt.show()


def run(cm: ConfigParser):
    methods = ["star"]
    conditions = [[24, 48], [48, 72], [24, 72]]

    txt_cond = []
    data = []
    for m in methods:
        txt_cond.extend([f"{x[0]}v{x[1]} {m}" for x in conditions])
        data.extend(get_data(cm, m, conditions, lab="winata"))

    plot_all_volcano(txt_cond, data, max_cols=2)
