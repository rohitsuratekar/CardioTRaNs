"""
Project: CardioTrans
Author: Rohit Suratekar
Year: 2019

All zfin analysis functions should go here.
"""

from constants.zfin_constants import *
from helpers.parsers import *

# Following disables Pandas chained assignment warnings
pd.options.mode.chained_assignment = None


def filter_with_gene(data: pd.DataFrame, gene: str,
                     allow_substring: bool = False) -> pd.DataFrame:
    """
    Get all expression details from ZFIN database

    :param gene: Name of the gene
    :param data: Pandas DataFrame
    :param allow_substring: If True, substring matches will be returned
    :return: Either Pandas DataFrame or List of ZFINExpression objects
    """
    if data.empty:
        return data

    if allow_substring:
        return data[
            data[COL_EXP_1_GENE_SYMBOL].str.contains(
                gene.strip().lower())].reset_index(drop=True)

    return data[
        data[COL_EXP_1_GENE_SYMBOL] == gene.strip().lower()].reset_index(
        drop=True)


def filter_with_hours(
        data: pd.DataFrame,
        start_hour: float = 0,
        end_hour: float = 30000):
    """
    Returns Pandas.DataFrame after filtering. Start and End are not defined,
    default value of 0 and 30000 will be used respectively.

    :param data: Pandas DataFrame of FILE_ZFIN_EXPRESSION file with
    COL_EXP_ALL as column names (similar to get_zfin_expression_dataframe() parser)
    :param start_hour: Starting Developmental Hour
    :param end_hour: Ending Developmental Hour
    :return: Pandas.DataFrame after filtering conditions
    """

    if data.empty:
        return data

    # Get developmental stage details
    stages = {x.name: x for x in get_zfin_stage_ontology()}

    # Convert starting hour column into ZFINStages objects
    data[COL_EXP_7_START_STAGE] = data[COL_EXP_7_START_STAGE].map(stages)

    # Sort the rows according to start and end values.
    # IMP: We used only start stage value as a cut off
    data = data[data[COL_EXP_7_START_STAGE].apply(
        lambda v: start_hour <= v.start <= end_hour)]

    # Bring back the names of the stages
    data[COL_EXP_7_START_STAGE] = data[COL_EXP_7_START_STAGE].apply(lambda v:
                                                                    v.name)
    return data.reset_index(drop=True)


def filter_with_structure(data: pd.DataFrame, structures: list or str = None,
                          check_in_parent: bool = True) -> pd.DataFrame:
    """
    Filters pandas DataFrame with structure or list of structure
    :param data: Pandas.DataFrame
    :param structures: single structure or list of structures
    :param check_in_parent: Check in the super structure (if available)
    :return: Pandas.DataFrame
    """

    if structures is None or len(structures) == 0 or data.empty:
        return data

    parent = "parent"  # New column name
    # If provided input is not a list, convert it to list
    if type(structures) != list:
        structures = [structures]

    # All the types will be converted into string
    structures = [str(x).strip().lower() for x in structures]

    # Get all anatomy items
    anatomy_items = {x.id: x for x in get_zfin_anatomy_items()}

    # Generate one more column as parent ID column
    data[parent] = data[COL_EXP_3_SUPER_STR_ID].apply(
        lambda v: anatomy_items[v].parent_id)

    # Convert Parent IDs into structure names
    data[parent] = data[parent].apply(lambda v: anatomy_items[v].name)

    # We will search given structure name in super structure, sub structure
    # or parent column

    def check_if_exists(r):
        """
        Simple function to check if all the conditions are satisfied
        :param r: Row
        :return: True or False
        """
        s1 = r[COL_EXP_4_SUPER_STR_NAME]
        s2 = r[COL_EXP_6_SUB_STR_NAME]
        s3 = r[parent]
        for s in structures:
            if type(s1) == str:
                if s in s1:
                    return True
            if type(s2) == str:
                if s in s2:
                    return True
            if type(s3) == str and check_in_parent:
                if s in s3:
                    return True
        return False

    bool_list = []
    for index, row in data.iterrows():
        bool_list.append(check_if_exists(row))

    # Filter rows and generate new index
    data = data[bool_list].reset_index(drop=True)
    # Drop parent column which was made extra
    data = data.drop([parent], axis=1)
    return data.reset_index(drop=True)


def run():
    data = get_zfin_expression_dataframe()
    data = filter_with_hours(data, 23, 25)
    data = filter_with_structure(data, "heart")
    print(data[COL_EXP_1_GENE_SYMBOL].to_list())
