"""
CardioTRaNs 2019
Author: Rohit Suratekar

All file parsers and converters
"""

import os

from constants import *
from local.operations import LocalDatabase
from models.zfin import *


def parse_zfin_anatomy_items() -> list:
    """
    Parse ZFIN anatomy file (FILE_ANATOMY_ITEMS) and convert into list of
    ZFINAnatomy objects
    :return: List of ZFINAnatomy object
    """
    models = []
    with open(DATA_FOLDER + FILE_ANATOMY_ITEMS) as f:
        for line in f:
            items = line.strip().split("\t")
            if len(items) == 4:
                if items[0].startswith("ZFA"):  # Skips Header
                    models.append(ZFINAnatomy(items))

    return models


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


def make_database():
    """
    WARNING : DELETES the existing database and generates new.
    """
    os.remove(DATABASE)
    print("Existing database removed.")
    print("Started Making Database...")
    db = LocalDatabase()
    db.add_anatomy_items(parse_zfin_anatomy_items())
    print("Task 1 of 4 completed...")
    db.add_stages(parse_stage_ontology())
    print("Construction of new Database successful")


def run():
    db = LocalDatabase()
    print(db.get_all_stages())
