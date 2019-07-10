"""
Project: CardioTrans
Author: Rohit Suratekar
Year: 2019

All file parsers and related functions should go here.
"""

import os.path

import pandas as pd

from constants.boolean_constants import *
from constants.ngs_constants import *
from constants.system_constants import *
from constants.zfin_constants import COL_EXP_ALL
from models.ngs_models import STGeneExpression
from models.zfin_models import *


def get_zfin_stage_ontology() -> list:
    """
    Parse ZFIN stage ontology file (FILE_ZFIN_STAGE_ONTOLOGY) and converts
    it into list of ZFINStages objects

    :return: List of ZFINStages objects
    """
    models = []
    with open(DATA_FOLDER + FILE_ZFIN_STAGE_ONTOLOGY) as f:
        for line in f:
            items = line.strip().split("\t")
            if len(items) == 5:
                if items[0].startswith("ZDB"):  # Skips Header
                    models.append(ZFINStages(items))
    return models


def get_zfin_anatomy_relations() -> dict:
    """
    Parse anatomy relation file (FILE_ZFIN_ANATOMY_RELATIONSHIP) and returns
    the list of items

    :return: List of parent-child-relation information with child as a key
    """
    relations = {}
    with open(DATA_FOLDER + FILE_ZFIN_ANATOMY_RELATIONSHIP) as f:
        for line in f:
            items = line.strip().split("\t")
            if len(items) == 3:
                if items[0].startswith("ZFA"):  # Skips Header
                    if items[2] == "part of":
                        relations[items[1]] = items
    return relations


def get_zfin_anatomy_items() -> list:
    """
    Parse ZFIN anatomy file (FILE_ANATOMY_ITEMS) and convert into list of
    ZFINAnatomy objects

    :return: List of ZFINAnatomy object
    """
    models = []
    relations = get_zfin_anatomy_relations()
    with open(DATA_FOLDER + FILE_ZFIN_ANATOMY_ITEMS) as f:
        for line in f:
            items = line.strip().split("\t")
            if len(items) == 4:
                if items[0].startswith("ZFA"):  # Skips Header
                    models.append(ZFINAnatomy(items, relations))

    # Additional item for missing entries
    models.append(ZFINAnatomy(["ZFA:0100000", "zebrafish anatomical entity",
                               "None", "None"], relations))
    return models


def get_zfin_expression_objects():
    """
    Parse wild type expression data and converts it into list of
    ZFINExpression objects

    :return: List of ZFINExpression objects
    """
    stages = {}
    for s in get_zfin_stage_ontology():
        stages[s.name] = s.id

    exp = []
    with open(DATA_FOLDER + FILE_ZFIN_EXPRESSION) as f:
        for line in f:
            items = line.strip().split("\t")
            if len(items) == 15:
                if items[0].startswith("ZDB"):  # Skips Header
                    exp.append(ZFINExpression(items, stages))

    return exp


def get_zfin_expression_dataframe() -> pd.DataFrame:
    """
    :return: Pandas DataFrame of expression data

    IMPORTANT: Use proper number for "skiprows" variable. This will skip top
    rows from the CSV file. Check your csv file and use it accordingly.
    ALWAYS check if column order is correct before proceeding to use this
    function for analysis.
    """
    with open(DATA_FOLDER + FILE_ZFIN_EXPRESSION) as f:
        return pd.read_csv(f, delimiter="\t", skiprows=0, names=COL_EXP_ALL)


def get_specific_rna_seq(file_name: str):
    """
    Gene Expression data from specific file
    :param file_name: Name of the file
    :return: pd.Dataframe
    """
    with open(DATA_FOLDER + file_name) as f:
        return pd.read_csv(f, delimiter="\t")


def get_rna_seq_data() -> pd.DataFrame:
    """
    Gene Expression data from the TSV file output from StringTie program

    This will average out all the TPM values from the available files
    :return: pd.Dataframe
    """

    # Get all the data
    d1 = get_specific_rna_seq(FILE_RNA_SEQ_GENE_EXPRESSION1)
    d2 = get_specific_rna_seq(FILE_RNA_SEQ_GENE_EXPRESSION2)
    d3 = get_specific_rna_seq(FILE_RNA_SEQ_GENE_EXPRESSION3)

    # Create the dataframe to calculate name dictionary

    n = pd.concat([d1, d2, d3])
    n = n.drop_duplicates(COL_STRING_TIE_1_GENE_ID, keep="first")
    n = n.set_index(COL_STRING_TIE_1_GENE_ID).T.to_dict('list')
    data = {x: n[x] for x in n}

    # Delete other columns after getting gene names

    d1 = d1[[COL_STRING_TIE_1_GENE_ID, COL_STRING_TIE_9_TPM]]
    d2 = d2[[COL_STRING_TIE_1_GENE_ID, COL_STRING_TIE_9_TPM]]
    d3 = d3[[COL_STRING_TIE_1_GENE_ID, COL_STRING_TIE_9_TPM]]

    # Take only transcripts which are with highest TPM in the case of
    # multiple transcripts

    d1 = d1.groupby(COL_STRING_TIE_1_GENE_ID, as_index=False).max()
    d2 = d2.groupby(COL_STRING_TIE_1_GENE_ID, as_index=False).max()
    d3 = d3.groupby(COL_STRING_TIE_1_GENE_ID, as_index=False).max()

    # Merge all transcript to make average

    m = d1.merge(d2, on=COL_STRING_TIE_1_GENE_ID)
    m = m.merge(d3, on=COL_STRING_TIE_1_GENE_ID)
    m["avg"] = m.mean(axis=1)
    m = m.sort_values("avg", ascending=False)

    m = m[[COL_STRING_TIE_1_GENE_ID, "avg"]]

    # Add other columns.
    # Coverage and FPKM columns are removed because of averaging

    m[COL_STRING_TIE_2_GENE_NAME] = m[COL_STRING_TIE_1_GENE_ID].apply(
        lambda x: data[x][0])

    m[COL_STRING_TIE_3_REFERENCE] = m[COL_STRING_TIE_1_GENE_ID].apply(
        lambda x: data[x][1])

    m[COL_STRING_TIE_4_STRAND] = m[COL_STRING_TIE_1_GENE_ID].apply(
        lambda x: data[x][2])

    m[COL_STRING_TIE_5_START] = m[COL_STRING_TIE_1_GENE_ID].apply(
        lambda x: data[x][3])

    m[COL_STRING_TIE_6_END] = m[COL_STRING_TIE_1_GENE_ID].apply(
        lambda x: data[x][4])

    m[COL_STRING_TIE_7_COVERAGE] = 0
    m[COL_STRING_TIE_8_FPKM] = 0
    m[COL_STRING_TIE_9_TPM] = m["avg"]
    del m["avg"]
    return m


def __generate_interested_gene_data():
    """
    Generates new files with only genes required for the Boolean Modelling
    analysis.
    """
    data = get_rna_seq_data()
    with open(DATA_FOLDER + FILE_RNA_SEQ, 'w') as f:
        header = []
        for a in data:
            header.append(a)

        header = "\t".join(header)
        print(header, file=f)

        for d in INTERESTED_GENES:
            k = data[data[COL_STRING_TIE_2_GENE_NAME] == d]
            if len(k.values) == 1:
                print("\t".join([str(x) for x in k.values[0]]), file=f)
            else:
                raise Exception(
                    "Zero or multiple gene found with {}".format(d))


def get_boolean_genes() -> pd.DataFrame:
    """
    Takes the averaged data and make extra column with generated
    STGeneExpression object
    :return: pd.DataFrame
    """

    def _make_pandas_objects(df: pd.DataFrame):
        if len(df.values) != 1:
            print(df)
            raise Exception("Only 1 gene is required for this conversion")
        return STGeneExpression(df.values[0])

    if not os.path.isfile(DATA_FOLDER + FILE_RNA_SEQ):
        print("File doesn't exist and generating new")
        __generate_interested_gene_data()

    with open(DATA_FOLDER + FILE_RNA_SEQ) as f:
        data = pd.read_csv(f, delimiter="\t")

    data[BOOL_OBJ_COL] = data[COL_STRING_TIE_2_GENE_NAME].apply(
        lambda c: _make_pandas_objects(
            data[data[COL_STRING_TIE_2_GENE_NAME] == c]))
    return data


def run():
    get_boolean_genes()
