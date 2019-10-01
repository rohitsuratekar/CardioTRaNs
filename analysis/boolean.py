#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
#  Boolean Modelling analysis


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

    c.input(NOT(c))

    t = TRN()
    t.add(c)
    return t


def run():
    base = test_network()
    base.update(3)
    base.print_expression_pattern()
