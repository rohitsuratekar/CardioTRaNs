#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
#  All expression pattern visualization will go here

import pandas as pd
from SecretPlots import BooleanPlot


def _plot_expression(data: pd.DataFrame, threshold):
    b = BooleanPlot(data.values.T, threshold=threshold)
    b.x_labels = data.columns
    b.y_labels = data.index
    b.add_legend()
    b.show()


def run():
    d = {
        "a": [0, 1, 2, 8],
        "b": [1, 20, 0, 0]
    }
    d = pd.DataFrame(d)
    _plot_expression(d, 2)
