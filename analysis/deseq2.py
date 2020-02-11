#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
# All DESeq2 related analysis

from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SecretColors import Palette

from constants.outputs import *
from helper import NameResolver

palette = Palette(show_warning=False)


def _validate(df: pd.DataFrame, *args):
    try:
        for a in args:
            m = df[a]
    except KeyError as e:
        raise KeyError(f"Your file should contain following column {args}. If "
                       f"you "
                       f"are using modified file, please specify these "
                       f"column names as an argument in the function. "
                       f"Please check documentation for more details") from e


def volcano_plot(filename: str, *,
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
    d = pd.read_csv(filename)
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


def ma_plot(filename: str, *,
            count_col: str = DESEQ2_BASE_MEAN,
            padj_col: str = DESEQ2_PADJ,
            lfc_col: str = DESEQ2_LOG2_CHANGE,
            color_col: str = None,
            color: str = palette.gray(shade=30),
            accent_color: str = palette.red(),
            threshold: float = 0.05,
            use_threshold: bool = True,
            ax: plt.Axes = None,
            decorations: bool = True,
            **kwargs):
    data = pd.read_csv(filename)
    _validate(data, count_col, padj_col, lfc_col)
    if color_col is None:
        color_col = f"clr{np.random.randint(10000, 90000)}"
        data[color_col] = color
        if use_threshold:
            data.loc[data[padj_col] < threshold, color_col] = accent_color

    if ax is None:
        _, ax = plt.subplots()

    ax.scatter(data[count_col], data[lfc_col], color=data[color_col], **kwargs)

    if decorations:
        ax.axhline(0, ls="--", color=palette.black())
        ax.set_ylabel("Log$_2$ Fold Change")
        ax.set_xlabel("Average Normalized Counts")


def ma_diff_plot(filename1: str, filename2: str, *,
                 count_col: str = DESEQ2_BASE_MEAN,
                 padj_col: str = DESEQ2_PADJ,
                 lfc_col: str = DESEQ2_LOG2_CHANGE,
                 gene_col: str = "gene_id"):
    data1 = pd.read_csv(filename1)
    data2 = pd.read_csv(filename2)
    _validate(data1, count_col, padj_col, lfc_col, gene_col)
    _validate(data2, count_col, padj_col, lfc_col, gene_col)

    dfs = [data1, data2]
    dfs = [x.set_index(gene_col) for x in dfs]
    print(dfs[0].loc["ENSDARG00000102141"][count_col])
    print(dfs[1].loc["ENSDARG00000102141"][count_col])


def run():
    nr = NameResolver("config.json")
    ma_plot(nr.deseq2_results("salmon",
                              lab="winata",
                              condition=[24, 48],
                              lfc=True),
            marker=".")

    plt.xlim(0, 1000)
    plt.show()
