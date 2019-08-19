#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 19/07/19, 2:18 PM
#
#  Copyright (c) 2019.
#
#  All file parsers related to NGS data
import os

import pandas as pd

from constants.ngs import *
from constants.system import *


def _get_ngs_file_names(hour: int, genotype: str, all_files: bool = True):
    """Returns files based on hour of development

    This is simple function which returns file(s) names based on hour provided
    Currently supports 20, 24, 48 and 72 hpf

    >>> _get_ngs_file_names(24, "wt")
    >>> # Returns list of files available for 24hpf for wild type

    :param hour: Hour Post Fertilization
    :param all_files: If True, all available files in that hours are
    returned as a list. If False, only single file is returned
    :return: Single file name or list of file names
    """

    g = genotype.strip().lower()

    if hour == 20:
        if g == GENOTYPE_WT:
            return FILE_RNA_SEQ_20H_1
    if hour == 24:
        if g == GENOTYPE_WT:
            if all_files:
                return [FILE_RNA_SEQ_24H_1, FILE_RNA_SEQ_24H_2]
            else:
                return FILE_RNA_SEQ_24H_1
    if hour == 48:
        if g == GENOTYPE_WT:
            if all_files:
                return [FILE_RNA_SEQ_48H_1, FILE_RNA_SEQ_48H_2]
            else:
                return FILE_RNA_SEQ_48H_1
    if hour == 72:
        if g == GENOTYPE_WT:
            if all_files:
                return [FILE_RNA_SEQ_72H_1, FILE_RNA_SEQ_72H_2]
            else:
                return FILE_RNA_SEQ_72H_1
        elif g == GENOTYPE_GATA5:
            if all_files:
                return [FILE_MT_GATA5_72H_1, FILE_MT_GATA5_72H_2]
            else:
                return FILE_MT_GATA5_72H_1

    raise Exception(
        "File is not available for this hour ({}) and genotype ({}) "
        "combination".format(hour, genotype))


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
    will be dropped and only row with highest TPM will be taken

    Note: This will reduce number of rows if there are any multiple
    transcripts present.

    :param data: List of DataFrames to be average out
    :return: Average Pandas DataFrame
    """

    # If its single DataFrame, convert it into the list
    if type(data) != list:
        data = [data]

    # Drop the duplicates (which are generally different Transcripts) and
    # take the one with highest TRP

    new_data = [
        x.sort_values(by=COL_STRING_TIE_9_TPM, ascending=False).set_index(
            COL_STRING_TIE_1_GENE_ID).groupby(
            level=COL_STRING_TIE_1_GENE_ID).first() for x in data]

    def _take_avg(d, col):
        """
        :param d: Data
        :param col: Column to consider for averaging
        :return: Pandas DataFrame with averaged column
        """
        n_d = [x[col] for x in d]
        m = pd.concat(n_d, axis=1, sort=False).fillna(0).mean(
            axis=1).to_frame(col)
        return m

    # Take average of TPM, FPKM and Coverage
    tpm = _take_avg(new_data, COL_STRING_TIE_9_TPM)
    fpkm = _take_avg(new_data, COL_STRING_TIE_8_FPKM)
    coverage = _take_avg(new_data, COL_STRING_TIE_7_COVERAGE)

    # Delete old columns to add new averaged columns
    del new_data[0][COL_STRING_TIE_9_TPM]
    del new_data[0][COL_STRING_TIE_8_FPKM]
    del new_data[0][COL_STRING_TIE_7_COVERAGE]

    # Concatenate the new averaged columns
    f = pd.concat([new_data[0], coverage, fpkm, tpm], axis=1, sort=False)

    # Just sanity check to warn if there are any duplicates which are not
    # incorporated

    if f.isna().sum().sum() != 0:
        print("Warning: Some of the values are not common in given DataFrames")

    # Get gene ID back as a column
    f.reset_index(inplace=True)

    return f


def get_rna_seq_data(hour: int, genotype: str):
    """ Returns average RNA-seq data for given hour
    :param genotype: Genotype of the data
    :param hour: Hour post fertilization
    :return: Pandas DataFrame with averaged coverage, fpkm and tpm
    """
    return average_data(
        convert_to_dataframe(_get_ngs_file_names(hour, genotype)))


def get_mapping_property(name: str):
    """
    Extracts information from the STAR Mapping log files
    :param name: EXACT name of the property (space sensitive)
    :return: List of properties
    """
    data = []

    def __extract_property(p):
        with open(MAPPING_FOLDER + p) as fl:
            for line in fl:
                if name.lower().strip() in line.lower().strip():
                    data.append(line.split("|")[1].strip())

    files = os.listdir(MAPPING_FOLDER)
    for f in files:
        if MAPPING_FILE_SUFFIX in f:
            __extract_property(f)

    return data


def run():
    get_mapping_property("Uniquely mapped reads %")
