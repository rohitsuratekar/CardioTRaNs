#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 29/07/19, 2:10 PM
#
#  Copyright (c) 2019.
#
# All Boolean related analysis functions will go here

from analysis.ngs import filter_with_gene
from constants.boolean import *
from constants.ngs import *
from helpers.parsers.boolean import get_boolean_expression


def is_expressed(gene: str, hour: int):
    """
    Checks if given gene is expressed at given time point
    :param gene: Name of the gene
    :param hour: Hours Post Fertilization
    :return: True or False depending on TPM_CUT_OFF and FPKM_CUT_OFF
    """
    d = get_boolean_expression(hour)
    d = filter_with_gene(d, gene)
    if len(d[COL_STRING_TIE_9_TPM].values) != 1:
        raise Exception("Some problem with retrieving '{}' at {} hpf. Either "
                        "there is no such gene or multiple instances of the "
                        "same genes are found. Also check your "
                        "INTERESTED_GENES constant and re-run the "
                        "parser if needed".format(gene, hour))

    return (d[COL_STRING_TIE_9_TPM].values[0] > TPM_CUT_OFF) and \
           (d[COL_STRING_TIE_8_FPKM].values[0] > FPKM_CUT_OFF)


def run():
    for g in BASE_GENES:
        k = is_expressed(g, 72)
        print("{} \t {}".format(g, k))
