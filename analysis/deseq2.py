#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
# All analysis related to DESeq2

from typing import Union

import matplotlib.gridspec as gridspec
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from SecretColors import Palette
from matplotlib.patches import Patch

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


def volcano_plot(data: Union[str, pd.DataFrame],
                 ax: plt.Axes = None,
                 lfc_col: str = DESEQ2_LOG2_CHANGE,
                 pad_col: str = DESEQ2_PADJ,
                 threshold: float = 0.05,
                 use_threshold: bool = True,
                 color_col: str = None,
                 normal_color: Union[str, tuple] = Palette().blue(),
                 sign_color: Union[str, tuple] = Palette().red(),
                 scatter_options: dict = None,
                 add_labels: bool = True
                 ):
    if isinstance(data, str):
        d = pd.read_csv(data)
    elif isinstance(data, pd.DataFrame):
        d = data.copy(deep=True)
    else:
        raise TypeError("Unknown format. Either provide path to the data "
                        "file or Pandas Dataframe of the same.")

    # Sanity check
    try:
        t = d[lfc_col], d[pad_col]  # Temporary variable to check
    except KeyError as e:
        raise KeyError(f"Please check your input data. Your input data "
                       f"should have columns for log fold change "
                       f"'{DESEQ2_LOG2_CHANGE}' "
                       f"and adjusted p value '{DESEQ2_PADJ}' column. If "
                       f"your column names are different, please provide "
                       f"them with argument 'lfc_col' and 'pad_col'"
                       f"") from e

    d = d.fillna(0)  # type: pd.DataFrame
    # Generate temporary column name
    temp = f"temp_col{np.random.randint(1000, 9000)}"
    colors = f"temp_color{np.random.randint(1000, 8000)}"

    d[temp] = d[pad_col].map(lambda x: -np.log10(x))
    d[colors] = normal_color

    # Use different color for values below threshold
    if use_threshold:
        d.loc[d[pad_col] < threshold, colors] = sign_color

    # If explicit color column is given, use those colors
    if color_col is not None:
        colors = color_col

    # If no axes is given, generate default one
    if ax is None:
        _, ax = plt.subplots()

    opts = {
        "marker": "."
    }
    if scatter_options is not None:
        opts = {**opts, **scatter_options}

    ax.scatter(d[lfc_col],
               d[temp], color=d[colors],
               **opts)

    if add_labels:
        ax.set_ylabel("-Log$_{10}$ Adj P value")
        ax.set_xlabel("Log$_2$ Fold change")


def plot_all_volcano(condition_array: list, dfs: list, max_cols=3):
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
    plt.show()


def multi_plots(cm: ConfigParser):
    methods = ["star"]
    conditions = [[24, 48], [48, 72]]

    txt_cond = []
    data = []
    for m in methods:
        txt_cond.extend([f"{x[0]}v{x[1]} {m}" for x in conditions])
        data.extend(get_data(cm, m, conditions, lab="winata"))

    plot_all_volcano(txt_cond, data)


def pca_plot(cm):
    p = Palette()
    from sklearn.decomposition import PCA
    conditions = [24, 48]

    pca = PCA()
    fm = FileManager(cm)
    df = fm.deseq2_counts("star", lab="winata", conditions=conditions)
    a1 = pca.fit(df.values[:, 1:])
    var = a1.components_
    percentages = [round(x, 2) for x in a1.explained_variance_ratio_ * 100]

    colors = [p.red()] * int((len(df.columns) - 1) / 2)
    colors.extend([p.blue()] * int((len(df.columns) - 1) / 2))
    lc = [Patch(label=f"{conditions[0]} hpf", color=p.red()),
          Patch(label=f"{conditions[1]} hpf", color=p.blue())]
    plt.legend(handles=lc, loc=0)
    plt.scatter(var[0, :], var[1, :], color=colors)
    plt.xlabel(f"PC1: {percentages[0]} %")
    plt.ylabel(f"PC2: {percentages[1]} %")
    plt.show()


def convert_to_symbol(cm: ConfigParser, df: pd.DataFrame) -> pd.DataFrame:
    print(df)


def run(cm: ConfigParser):
    dfs = get_data(cm, "star", conditions=[[24, 48]], lab="winata")
    convert_to_symbol(cm, dfs[0])
