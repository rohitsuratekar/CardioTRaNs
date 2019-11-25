#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
#  Boolean Modelling analysis


from SecretPlots import BooleanPlot

from models.biology import *
from models.network import *


def base_network():
    nkx2_5 = Gene("Nkx2.5", True)
    gata4 = Gene("Gata4", True)
    tbx5a = Gene("Tbx5a", True)
    mef2c = Gene("Mef2ca", True)
    hand2 = Gene("Hand2", True)

    nkx2_5.input(OR(gata4, tbx5a))
    tbx5a.input(COPY(nkx2_5))
    gata4.input(COPY(tbx5a))
    mef2c.input(OR(nkx2_5, gata4, tbx5a))
    hand2.input(OR(gata4, mef2c))

    t = TRN()
    t.add(nkx2_5, gata4, tbx5a, mef2c, hand2)
    return t


def test_network():
    a = Gene("a", False)
    b = Gene("b", True)
    c = Gene("c", True)

    c.input(OR(a, OR(b, a)))

    t = TRN()
    t.add(c, b, a)
    return t


def visualize(base: TRN, iterations: int):
    base.update(iterations)
    base.print_network()
    pattern, names = base.expression_pattern()
    b = BooleanPlot(pattern, 1)
    (b
     .change_orientation("y")
     .invert_y()
     .add_y_midlines(color="w", alpha=0.5)
     .add_x_midlines(color="w", alpha=0.5)
     .add_y_ticklabels(names)
     .add_x_label("Iterations")
     .add_x_padding(0, 0)
     .add_y_padding(0, 0)
     .change_aspect_ratio(1)
     .add_legends(loc='center left', bbox_to_anchor=(1, 0.5))
     .show(tight=True))


def run():
    tn = base_network()
    visualize(tn, 2)
