#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
#   Analysis related to Expression atlas data
#   https://www.ebi.ac.uk/gxa/experiments/E-ERAD-475/Results

import pandas as pd
import yaml

from constants.zfin import *


def _get_file_data():
    with open("config.yml") as f:
        data = yaml.load(f, Loader=yaml.SafeLoader)
    return data


def _get_df(atlas: str, zfin: str):
    df = pd.read_csv(atlas, sep="\t")
    onto = pd.read_csv(zfin, sep="\t")
    del df["Gene ID"]
    df = df.rename(
        columns={"zygote": "zygote 1-cell",
                 "larval protruding mouth": "larval protruding-mouth"})
    df = df.melt(id_vars="Gene Name")
    df = df.rename(columns={"variable": ZFIN_ONT_BEGIN_HOUR, "value": "TPM"})
    onto_map = dict(zip(onto[ZFIN_ONT_STAGE_NAME],
                        onto[ZFIN_ONT_BEGIN_HOUR]))
    onto_map = {x.lower().replace(":", " "): onto_map[x] for x in onto_map}
    df[ZFIN_ONT_BEGIN_HOUR] = df[ZFIN_ONT_BEGIN_HOUR].map(
        lambda x: onto_map[x])
    df = df.fillna(0)
    return df


def get_expressing_genes(hour):
    data = _get_file_data()
    atlas_file = data["atlas"]
    zfin_file = data["zfin"]["ontology"]
    df = _get_df(atlas_file, zfin_file)
    df = df[df[ZFIN_ONT_BEGIN_HOUR] == hour]
    print(df)


def run():
    get_expressing_genes(72)
