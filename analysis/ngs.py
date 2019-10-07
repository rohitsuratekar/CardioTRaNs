#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# NGS Analysis

from constants.ngs import *
from constants.system import *
from helpers.ngs_parser import get_string_tie_data, average_data


def test():
    d = get_string_tie_data(GENOTYPE_WT, 24, BIO_PROJECT_WINATA_LAB)
    d = average_data(d, OUTPUT_STRING_TIE)
    print(d)


def run():
    test()
