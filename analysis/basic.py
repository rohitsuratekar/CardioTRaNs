"""
CardioTRaNs 2019
Author: Rohit Suratekar

Data analysis of the ZFIN data-set.

Important Note:
You will need to generate database (local.operations.LocalDatabase) before
using few of the functions listed in this files. To make database use
`helpers.parsers.make_database()` . You should use this function once after
updating database files. Working with HDF5 data-set is faster than csv files.
However we will use regular parsing for expression analysis as our
benchmarking suggest it is faster than HDF5.
"""

from pandas import DataFrame

from helpers.parsers import *

# Following disables Pandas chained assignment warnings
pd.options.mode.chained_assignment = None


def get_genes_for_hours(
        data: DataFrame,
        start_hour: float = 0,
        end_hour: float = 30000):
    """
    Returns Pandas.DataFrame after filtering. Start and End are not defined,
    default value of 0 and 30000 will be used respectively.

    :param data: Pandas DataFrame of FILE_WILD_TYPE_EXPRESSION file with
    COL_EXP_ALL as column names (similar to get_expression_data_frame() parser)
    :param start_hour: Starting Developmental Hour
    :param end_hour: Ending Developmental Hour
    :return: Pandas.DataFrame after filtering conditions
    """

    # Get developmental stage details
    stages = {x.name: x for x in parse_stage_ontology()}

    # Convert starting hour column into ZFINStages objects
    data[COL_EXP_7_START_STAGE] = data[COL_EXP_7_START_STAGE].map(stages)

    # Sort the rows according to start and end values.
    # IMP: We used only start stage value as a cut off
    data = data[data[COL_EXP_7_START_STAGE].apply(
        lambda v: start_hour <= v.start <= end_hour)]

    # Bring back the names of the stages
    data[COL_EXP_7_START_STAGE] = data[COL_EXP_7_START_STAGE].apply(lambda v:
                                                                    v.name)
    data = data.reset_index(drop=True)  # Reset the DataFrame indices
    return data


def get_genes_from_structure(data: DataFrame, structures: list or str = None,
                             check_in_parent: bool = True):
    if structures is None or len(structures) == 0:
        return data

    parent = "parent"  # New column name
    # If provided input is not a list, convert it to list
    if type(structures) != list:
        structures = [structures]

    # All the types will be converted into string
    structures = [str(x).strip().lower() for x in structures]

    # Get local data-set (see header of this file for more information)
    db = LocalDatabase()
    # Get all anatomy items
    anatomy_items = {x.id: x for x in db.get_all_anatomy_items()}

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
    return data


def get_gr_expressed_genes(hours: int) -> list:
    if hours not in [24, 48, 72]:
        raise Exception("Data is available only for 24, 48 and 72 hours")
    df = get_all_gr_expressed_data()
    # Get only required hour genes and drop items where gene symbol is not
    # available

    # NOTE: In published file "ZFIN_ID" is actually gene symbol
    df = df[['ZFIN_ID', 'RNAseq_{}_plus_minus_log2FC'.format(hours)]].dropna()
    df = df.reset_index(drop=True)  # Reset Index
    return df['ZFIN_ID'].to_list()


def show_structures_with(substring: str):
    db = LocalDatabase()
    # Get all anatomy items
    anatomy_items = db.get_all_anatomy_items()
    for s in anatomy_items:
        if substring.lower().strip() in s.name.strip().lower():
            print("{}: \t{}".format(s.id, s.name))


def test():
    gene_name = "nkx2.5"
    data = parse_wt_expression()

    for d in data:
        if gene_name == d.gene_symbol:
            print(d.sub_structure_id)


def run():
    test()
