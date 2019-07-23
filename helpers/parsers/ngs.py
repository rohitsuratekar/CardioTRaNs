#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 19/07/19, 2:18 PM
#
#  Copyright (c) 2019.
#
#  All file parsers related to NGS data

import pandas as pd

from constants.ngs import *
from constants.system import *


def get_ngs_file_names(hour: int, all_files: bool = True):
    """Returns files based on hour of development

    This is simple function which returns file(s) names based on hour provided
    Currently supports 20, 24, 48 and 72 hpf

    >>> get_ngs_file_names(24) # Returns list of files available for 24hpf

    :param hour: Hour Post Fertilization
    :param all_files: If True, all available files in that hours are
    returned as a list. If False, only single file is returned
    :return: Single file name or list of file names
    """
    if hour == 20:
        return FILE_RNA_SEQ_20H_1
    if hour == 24:
        if all_files:
            return [FILE_RNA_SEQ_24H_1, FILE_RNA_SEQ_24H_2]
        else:
            return FILE_RNA_SEQ_24H_1
    if hour == 48:
        if all_files:
            return [FILE_RNA_SEQ_48H_1, FILE_RNA_SEQ_48H_2]
        else:
            return FILE_RNA_SEQ_48H_1
    if hour == 72:
        if all_files:
            return [FILE_RNA_SEQ_72H_1, FILE_RNA_SEQ_72H_2]
        else:
            return FILE_RNA_SEQ_72H_1

    raise Exception("File is not available for this hour. Try 20, 24, 48 or 72")


def convert_to_dataframe(files):
    """Returns Pandas DataFrame of given file(s)

    If list is provided, it will return list of DataFrames

    :param files: Name of the file(s). Either single string or list
    :return: Pandas DataFrame or list of DataFrames
    """

    def _convert(name):
        with open(DATA_FOLDER + name) as f:
            return pd.read_csv(f, delimiter="\t", names=ALL_STRING_TIE_COLS,
                               skiprows=1)

    if type(files) == str:
        return _convert(files)

    if type(files) == list:
        return [_convert(x) for x in files]

    raise Exception("Unexpected file format. It should be either single "
                    "file name or list of file names")


def average_data(data: list):
    """Averages the data

    Takes different pandas DataFrames as an input and converts them into
    single DataFrame by averaging out Coverage, FPKMs and TPMs. Strand , Ref,
    Start and End will be inherited from the first DataFrame.

    All analysis will be done on the basis of Gene IDs. Multiple Gene IDs
    will be dropped and only first occurrence will be considered for further
    analysis.

    :param data: List of DataFrames to be average out
    :return: Average Pandas DataFrame
    """

    def _take_avg(d, col):
        d_new = [
            x.sort_values(by=COL_STRING_TIE_9_TPM, ascending=False).set_index(
                COL_STRING_TIE_1_GENE_ID).groupby(
                level=COL_STRING_TIE_1_GENE_ID).first()[col] for x in
            d]
        m = pd.concat(d_new, axis=1, sort=False)
        return m

    # Extracts all gene ids and make dictionary for future use
    n = pd.concat(data)
    n = n.drop_duplicates(COL_STRING_TIE_1_GENE_ID, keep="first")
    n = n.set_index(COL_STRING_TIE_1_GENE_ID).T.to_dict('list')
    details = {x: n[x] for x in n}

    k = _take_avg(data, COL_STRING_TIE_9_TPM)
    print(k)


def get_rna_seq_data() -> pd.DataFrame:
    """
    Gene Expression data from the TSV file output from StringTie program

    This will average out all the TPM values from the available files
    :return: pd.Dataframe
    """

    # Get all the data
    d1 = (FILE_RNA_SEQ_24H_1)
    d2 = (FILE_RNA_SEQ_24H_2)
    d3 = (FILE_RNA_SEQ_20H_1)

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


def test():
    def take_avg(data, col):
        d_new = [
            x.sort_values(by="b", ascending=False).set_index("a").groupby(
                level="a").first()[col] for x in data]
        m = pd.concat(d_new, axis=1, sort=False).fillna(0).mean(axis=1)
        return m

    d1 = {"a": ["x", "y", "z", "x"], "b": [1, 1, 1, 11], "c": [1, 1, 1, 1]}
    d2 = {"a": ["x", "y", "z"], "b": [2, 2, 2], "c": [1, 1, 1]}
    d3 = {"a": ["x", "y", "zi"], "b": [3, 3, 3], "c": [1, 1, 1]}

    d1 = pd.DataFrame(d1)
    d2 = pd.DataFrame(d2)
    d3 = pd.DataFrame(d3)
    d = [d1, d2, d3]

    k1 = take_avg(d, "c")
    k2 = take_avg(d, "b")

    n = pd.concat([k1, k2], axis=1)
    n.columns = ["c", "b"]
    n.reset_index(inplace=True)
    n = n.rename(columns={"index": "a"})
    n = n[["b", "c", "a"]]

    print(n)


def run():
    d = convert_to_dataframe([FILE_RNA_SEQ_24H_1, FILE_RNA_SEQ_24H_2])
    average_data(d)
    # test()
