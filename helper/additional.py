#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
#   These are helper functions I made for myself. User should have their own

import os

import yaml


def get_runs(project_id) -> dict:
    with open("runs.yml") as f:
        data = yaml.load(f, Loader=yaml.SafeLoader)
    if project_id not in data.keys():
        raise KeyError(f"'runs.yml' file do not have runs from {project_id}. "
                       f"Available projects: {list(data.keys())}")
    return data[project_id]


def get_file_names(runs: dict, method: str, time: int, genotype: str):
    exp = []
    for key in runs:
        if runs[key]["time"] == time and runs[key]["genotype"] == genotype:
            exp.append(key)

    if len(exp) == 0:
        raise ValueError(f"No runs found for the combination time :{time} "
                         f"and genotype : {genotype} in current runs")

    with open("config.yml") as f:
        data = yaml.load(f, Loader=yaml.SafeLoader)

    if method not in data["rnaseq"].keys():
        raise KeyError(f"RNASeq method '{method}' not found in 'config.yml'")

    base_folder = data["rnaseq"][method]
    files = []
    for key in exp:
        name = f"{base_folder}/{key}/"
        if method == "salmon":
            name += "quant.sf"
        elif method == "stringtie":
            name += f"{key}_gene_expression.tsv"
        elif method == "kallisto":
            name += "abundance.tsv"
        elif method == "star":
            name += f"{key}_Aligned.sortedByCoord.out.bam"
        else:
            raise KeyError(f"Method '{method}' is not supported yet.")

        if os.path.exists(name):
            files.append(name)
        else:
            raise FileNotFoundError(f"Output file {name} not found")

    return files
