#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
#   All functions related to the STRING database analysis will go here

import subprocess
import warnings
from typing import List
from io import StringIO
import yaml
import pandas as pd

from constants.ourputs import STRING_COLUMNS


def _get_files() -> dict:
    with open("config.yml") as f:
        data = yaml.load(f, Loader=yaml.SafeLoader)
        return data["string"]


def _scan_content(filename: str,
                  content: str,
                  exact: bool = False) -> List[str]:
    opts = ["grep"]
    if exact:
        opts.append("-w")
    opts.append(content)
    opts.append(filename)

    results = subprocess.run(opts, capture_output=True,
                             universal_newlines=True)
    if results.returncode != 0:
        raise KeyError(f"Unable to find '{content}' in {filename}")
    return results.stdout.strip().split("\n")


def extract_name(info_file: str, gene: str, is_id: bool = False):
    """
    Check the gene name and returns its STRING ID or Vice Versa (by
    providing is_id = True)
    To avoid returning multiple genes which share common part of gene_id.
    Whenever is_id is true, exact match will be searched

    >>> extract_name("filename", "nkx2.5") # Returns Gene ID of nkx2.5
    Following retuns Gene name with ID "7955.ENSARP00044"
    >>> extract_name("filename", "7955.ENSARP00044", is_id=True)

    :param is_id: if True, searches Gene ID and returns Gene name
    :param info_file: Full path of STRING input file (usually such file will
    have name as 'XXXX.protein.info.v11.0.txt')
    :param gene: Name of the gene (e.g. nkx2.5) or ID (e.g.:
    7955.ENSDARP00000118756)
    :return: gene ID or gene Name used in other STRING files
    """
    if not isinstance(gene, str):
        raise ValueError("You should pass exactly one gene name as a 'String'")

    lines = _scan_content(info_file, gene, is_id)
    loc = int(is_id)
    name = lines[0].split("\t")[loc]
    lines = [x.split("\t")[loc] for x in lines]
    if len(lines) > 1:
        warnings.warn(f"Multiple gene names found for given query of gene"
                      f" '{gene}'. Retuning first result from {lines}")
    return name


def get_interactions(link_file: str, gene_id: str):
    links = _scan_content(link_file, gene_id, True)
    links = StringIO("\n".join(links))
    df = pd.read_csv(links, delim_whitespace=True, header=None,
                     names=STRING_COLUMNS)
    print(df.columns)


def run():
    files = _get_files()
    info = files["info"]
    actions = files["actions"]
    links = files["links"]
    nkx = extract_name(info, "slc35a5")
    get_interactions(links, nkx)
