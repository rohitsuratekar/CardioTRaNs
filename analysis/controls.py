#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# All control analysis

from analysis.atlas import *
from constants.zfin import *
from helpers.filemanager import *


def get_ontology_info(name):
    return zfin_ontology().set_index(name)[[
        ZFIN_ONT_BEGIN_HOUR, ZFIN_ONT_END_HOUR]].to_dict("index")


def full_expression_range():
    """
    This generates widest available expression data.
    For example, if there are entries about gene X expressed from 5 to 10 as
    well as 20-30, this function will assign expression range as 5 to 30
    """
    zfin = zfin_expression()

    ont = get_ontology_info(ZFIN_ONT_STAGE_NAME)

    zfin[ZFIN_EXP_START_STAGE] = (zfin[ZFIN_EXP_START_STAGE]
                                  .map(lambda x: ont[x][ZFIN_ONT_BEGIN_HOUR]))

    zfin[ZFIN_EXP_END_STAGE] = (zfin[ZFIN_EXP_END_STAGE]
                                .map(lambda x: ont[x][ZFIN_ONT_END_HOUR]))

    zfin = zfin.sort_values(ZFIN_EXP_START_STAGE)
    zfin = zfin.groupby(ZFIN_EXP_GENE_SYMBOL).agg({
        ZFIN_EXP_START_STAGE: lambda x: min([y for y in x]),
        ZFIN_EXP_END_STAGE: lambda x: max([y for y in x])})
    zfin = zfin.reset_index()
    return zfin


def zfin_genes_not_expressed(min_time, max_time):
    """
    Gets ZFIN genes which are not expressed within the range provided
    :param min_time: minimum time
    :param max_time: maximum time
    :return: pandas DataFrame
    """
    # Get data
    xpat = zfin_pattern()
    fish = zfin_wt_lines()
    assay = zfin_assay()
    ont = get_ontology_info(ZFIN_ONT_STAGE_ID)

    # Only take entries where expression is absent
    xpat = xpat[xpat[ZFIN_XPAT_EXPRESSION_FOUND] == "f"]

    # Convert developmental terms into time points
    xpat[ZFIN_XPAT_START_STAGE] = (xpat[ZFIN_XPAT_START_STAGE]
                                   .map(lambda x: ont[x][ZFIN_ONT_BEGIN_HOUR]))

    xpat[ZFIN_XPAT_END_STAGE] = (xpat[ZFIN_XPAT_END_STAGE]
                                 .map(lambda x: ont[x][ZFIN_ONT_END_HOUR]))

    # Filter by time
    xpat = xpat[xpat[ZFIN_XPAT_START_STAGE] >= min_time]
    xpat = xpat[xpat[ZFIN_XPAT_END_STAGE] <= max_time]
    xpat = xpat[ZFIN_XPAT_EXPRESSION_ID].values

    # Get only assay terms which are filtered above
    assay = assay[assay[ZFIN_ASSAY_EXPRESSION_ID].isin(xpat)].reset_index(
        drop=True)

    # Only take wild type assay terms
    assay = assay[assay[ZFIN_ASSAY_FISH_ID].isin(fish[ZFIN_FISH_ID].values)]

    # Get gene-lists in that assay
    gene_list = list(set(assay[ZFIN_ASSAY_GENE_SYMBOL].values))

    # Get widest expression pattern of ALL genes
    cross_check = full_expression_range()

    # Get gene which are not present in given range. These will be genes
    # which have expression pattern range either lower than or higher than
    # our input times
    cross_check = cross_check[(cross_check[ZFIN_EXP_START_STAGE] > max_time)
                              | (cross_check[ZFIN_EXP_END_STAGE] < min_time)]
    cross_check = cross_check[ZFIN_EXP_GENE_SYMBOL].values

    # Return genes which are present in assay list (where expression is
    # absent) as well as cross check list (where expression pattern range is
    # outside our input time)
    return list(set(gene_list).intersection(cross_check))


def atlas_not_expressed(min_time, max_time, tpm_cutoff):
    """
    Get genes from Expression Atlas whose expression is absent between
    min - max time

    :param min_time: Minimum time to start expression
    :param max_time: Maximum time till when expression should be there
    :param tpm_cutoff: Below which genes will be considered as not expressed
    :return: Pandas DataFrame
    """
    atlas = isolate_by_time(atlas_data(), min_time, max_time,
                            tpm_cutoff=tpm_cutoff,
                            strict=False,
                            not_expressed=True)[EXP_ATLAS_GENE_NAME].values
    return atlas


def not_expressed(min_time, max_time, tpm_cutoff):
    atlas = atlas_not_expressed(min_time, max_time, tpm_cutoff)
    print(
        f"Number of atlas genes : {len(atlas)} (with TPM cutoff {tpm_cutoff})")
    zfin = zfin_genes_not_expressed(min_time, max_time)
    print(f"Number of zfin genes : {len(zfin)}")
    genes = list(set(atlas).intersection(zfin))
    print(f"Number of common genes : {len(genes)}")
    print(f"Common genes : {', '.join(genes)}")


def run():
    not_expressed(23, 73, 1)
