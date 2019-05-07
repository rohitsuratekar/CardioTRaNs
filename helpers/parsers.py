"""
CardioTRaNs 2019
Author: Rohit Suratekar

All file parsers and converters.
IMPORTANT: Check the column order in ZFIN files whenever you are using
following parser functions. Following code is written based on column order
present in April 2019.

ZFA:0100000 entry was not available in FILE_ANATOMY_ITEMS, hence, it was
added to the data-set while parsing the file
"""

import os

import pandas as pd

from constants import *
from local.operations import LocalDatabase
from models.zfin import *


def parse_stage_ontology() -> list:
    """
    Parse ZFIN stage ontology file (FILE_STAGE_ONTOLOGY) and converts it into
    list of ZFINStages objects
    :return: List of ZFINStages objects
    """
    models = []
    with open(DATA_FOLDER + FILE_STAGE_ONTOLOGY) as f:
        for line in f:
            items = line.strip().split("\t")
            if len(items) == 5:
                if items[0].startswith("ZDB"):  # Skips Header
                    models.append(ZFINStages(items))
    return models


def parse_anatomy_relations() -> dict:
    """
    Parse anatomy relation file (FILE_ANATOMY_RELATIONSHIP) and returns the
    list of items
    :return: List of parent-child-relation information with child as a key
    """
    relations = {}
    with open(DATA_FOLDER + FILE_ANATOMY_RELATIONSHIP) as f:
        for line in f:
            items = line.strip().split("\t")
            if len(items) == 3:
                if items[0].startswith("ZFA"):  # Skips Header
                    if items[2] == "part of":
                        relations[items[1]] = items
    return relations


def parse_zfin_anatomy_items() -> list:
    """
    Parse ZFIN anatomy file (FILE_ANATOMY_ITEMS) and convert into list of
    ZFINAnatomy objects
    :return: List of ZFINAnatomy object
    """
    models = []
    relations = parse_anatomy_relations()
    with open(DATA_FOLDER + FILE_ANATOMY_ITEMS) as f:
        for line in f:
            items = line.strip().split("\t")
            if len(items) == 4:
                if items[0].startswith("ZFA"):  # Skips Header
                    models.append(ZFINAnatomy(items, relations))

    # Additional item for missing entries
    models.append(ZFINAnatomy(["ZFA:0100000", "zebrafish anatomical entity",
                               "None", "None"], relations))
    return models


def parse_wt_expression():
    """
    Parse wild type expression data and converts it into list of
    ZFINExpression objects
    :return: List of ZFINExpression objects
    """
    stages = {}
    for s in parse_stage_ontology():
        stages[s.name] = s.id

    exp = []
    with open(DATA_FOLDER + FILE_WILD_TYPE_EXPRESSION) as f:
        for line in f:
            items = line.strip().split("\t")
            if len(items) == 15:
                if items[0].startswith("ZDB"):  # Skips Header
                    exp.append(ZFINExpression(items, stages))

    return exp


def get_expression_data_frame() -> pd.DataFrame:
    """
    :return: Pandas DataFrame of expression data

    IMPORTANT: Use proper number for "skiprows" variable. This will skip top
    rows from the CSV file. As of April 2019, ZFIN have first line as a date
    on which data file is downloaded and second line as a headers. Since we
    are going to use our own column names, we will skip the column row.
    However ALWAYS check if column order is correct before proceeding to use
    this function for analysis.
    """
    with open(DATA_FOLDER + FILE_WILD_TYPE_EXPRESSION) as f:
        return pd.read_csv(f, sep="\t", skiprows=2, names=COL_EXP_ALL)


def make_database():
    """
    WARNING : DELETES the existing database and generates new.
    """
    try:
        os.remove(DATABASE)
        print("Existing database removed.")
    except FileNotFoundError:
        print("Database not found...")

    print("Started making new database...")
    db = LocalDatabase()
    db.add_anatomy_items(parse_zfin_anatomy_items())
    print("Task 1 of 2 completed...")
    db.add_stages(parse_stage_ontology())
    print("Task 2 of 2 completed...")
    print("Construction of new database successful")


def get_all_gr_expressed_data() -> pd.DataFrame:
    """
    Expression data from the supplement data of article Pawlak et. al. 2019
    (Pubmed ID: 30760547)
    :return: Pandas DataFrame
    """
    return pd.read_excel(DATA_FOLDER + FILE_DE_SEQ_GFP, skiprows=1)


def run():
    get_all_gr_expressed_data()
