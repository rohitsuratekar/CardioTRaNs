#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# These functions will analyse file output from R


import matplotlib.pyplot as plt
import pandas as pd
from SecretColors import Palette

from constants.biomart import *
from constants.other import *
from helpers.filemanager import r_out, biomart_genes


def get_data() -> pd.DataFrame:
    data = r_out()
    genes = biomart_genes()
    data = data.set_index(R_GENES)
    genes = genes.set_index(BIOMART_GENE_ID)[[BIOMART_GENE_NAME]]
    genes = genes.drop_duplicates()
    data = (
        data.join(genes, how="outer")
            .rename_axis(BIOMART_GENE_ID)
            .reset_index()
    )
    return data


def ma_plot():
    data = get_data()
    p = Palette()
    plt.scatter(data[R_BASE_MEAN].values,
                data[R_LOG2_FOLD_CHANGE].values,
                color=p.red(),
                marker=".")
    plt.axhline(0, ls="--", color=p.gray(shade=70))
    plt.show()


def run():
    ma_plot()
