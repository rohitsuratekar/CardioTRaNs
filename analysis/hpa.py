#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Analysis related to Human Protein Atlas
# https://www.proteinatlas.org/about/download

from helpers.filemanager import hpa_pathology
from constants.hpa import *


def get_data():
    d = hpa_pathology()
    print(d[d[HPA_PATHOLOGY_GENE_NAME] == "CHL1"][HPA_PATHOLOGY_CANCER])


def run():
    get_data()
