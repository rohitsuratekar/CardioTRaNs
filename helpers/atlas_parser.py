#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Functions related to the expression atlas


from constants.zfin import *
from helpers.filemanager import expression_atlas, zfin_ontology


def atlas_data():
    """
    :return: Normalized Expression Atlas data in the form of Pandas
    DataFrames. All the column names will be changed to hpf values of the
    respective developmental stages
    """
    stages = zfin_ontology()
    data = expression_atlas()
    stages[ZFIN_ONT_STAGE_NAME] = stages[ZFIN_ONT_STAGE_NAME].str.replace(
        ":", " ")
    stages[ZFIN_ONT_STAGE_NAME] = stages[ZFIN_ONT_STAGE_NAME].str.lower()

    new_names = {}
    for d in data.columns:
        t = stages[stages[ZFIN_ONT_STAGE_NAME] == d]
        if len(t) > 0:
            new_names[d] = t[ZFIN_ONT_BEGIN_HOUR].values[0]
        else:
            if d == "zygote":
                new_names[d] = 0.0
            elif d == "larval protruding mouth":
                new_names[d] = 72.0
            else:
                new_names[d] = d

    data = data.rename(columns=new_names)

    return data


def run():
    pass
