#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# All control analysis

from SecretColors import Palette
from SecretPlots import *

from analysis.atlas import *
from constants.zfin import *
from helpers.filemanager import *

WHOLE_ORGANISM = "ZFA:0001094"


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


def zfin_genes_expressed(times: list):
    """
    Returns zfin genes expressed
    :return:
    """
    zfin = zfin_expression()
    ont = get_ontology_info(ZFIN_ONT_STAGE_NAME)

    zfin[ZFIN_EXP_START_STAGE] = (zfin[ZFIN_EXP_START_STAGE]
                                  .map(lambda x: ont[x][ZFIN_ONT_BEGIN_HOUR]))

    zfin[ZFIN_EXP_END_STAGE] = (zfin[ZFIN_EXP_END_STAGE]
                                .map(lambda x: ont[x][ZFIN_ONT_END_HOUR]))

    # Filter by time
    zfin = zfin[(zfin[ZFIN_EXP_START_STAGE].isin(times))]
    zfin = (zfin.groupby(ZFIN_EXP_GENE_SYMBOL)
            .agg({ZFIN_EXP_START_STAGE: lambda x: len(set([y for y in x]))})
            )
    zfin = zfin[zfin[ZFIN_EXP_START_STAGE] == len(times)]
    zfin = zfin.reset_index()

    return zfin[ZFIN_ASSAY_GENE_SYMBOL].values


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

    # To be more strict, take only entries where whole organism is analyzed
    xpat = xpat[xpat[ZFIN_XPAT_ANATOMY_SUPER] == WHOLE_ORGANISM]

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


def atlas_expression(min_time, max_time, tpm_cutoff, not_expressed):
    """
    Get genes from Expression Atlas whose expression is absent between
    min - max time

    :param not_expressed: If True, it will return not expressed genes
    :param min_time: Minimum time to start expression
    :param max_time: Maximum time till when expression should be there
    :param tpm_cutoff: Below which genes will be considered as not expressed
    :return: Pandas DataFrame
    """
    atlas = isolate_by_time(atlas_data(), min_time, max_time,
                            tpm_cutoff=tpm_cutoff,
                            strict=False,
                            not_expressed=not_expressed)[
        EXP_ATLAS_GENE_NAME].values
    return atlas


def get_not_expressed(min_time, max_time, tpm_cutoff,
                      skip_print: bool = False):
    atlas = atlas_expression(min_time, max_time, tpm_cutoff, True)
    zfin = zfin_genes_not_expressed(min_time, max_time)
    if not skip_print:
        print(
            f"Number of atlas genes : {len(atlas)} (with TPM cutoff {tpm_cutoff})")
        print(f"Number of zfin genes : {len(zfin)}")

    genes = list(set(atlas).intersection(zfin))

    if not skip_print:
        print(f"Number of common genes : {len(genes)}")
        print(f"Common genes : {', '.join(genes)}")
    return genes


def get_expressed(min_time, max_time, tpm_cutoff, skip_print: bool = False):
    atlas = atlas_expression(min_time, max_time, tpm_cutoff, False)
    zfin = zfin_genes_expressed([24, 48, 72])
    if not skip_print:
        print(
            f"Number of atlas genes : {len(atlas)} (with TPM cutoff {tpm_cutoff})")
        print(f"Number of zfin genes : {len(zfin)}")

    genes = list(set(atlas).intersection(zfin))
    if not skip_print:
        print(f"Number of common genes : {len(genes)}")
        print(f"Common genes : {', '.join(genes)}")
    return genes


def tpm_effect_on_negative_control(start, end):
    values = []
    locs = list(range(5))
    for x in locs:
        atlas = atlas_expression(start, end, x, True)
        zfin = zfin_genes_not_expressed(start, end)
        values.append(len(list(set(atlas).intersection(zfin))))

    bar = BarPlot(values)
    bar = (bar
           .add_x_label("TPM Values")
           .add_x_ticklabels(locs)
           .add_y_label("No of genes found")
           .add_x_gap(0)
           )
    bar.draw()
    bar.ax.set_title("Effect of TPM on number of Negative control genes")
    bar.ax.axhline(56, ls="--", color="k")
    bar.ax.set_ylim(0, 65)
    bar.show()


def tpm_effect_on_positive_control(start, end):
    p = Palette()
    values = []
    locs = list(range(5))
    for x in locs:
        atlas = atlas_expression(start, end, x, False)
        zfin = zfin_genes_expressed([24, 48, 72])
        values.append(len(list(set(atlas).intersection(zfin))))

    bar = BarPlot(values)
    bar = (bar
           .add_x_label("TPM Values")
           .add_x_ticklabels(locs)
           .add_y_label("No of genes found")
           .add_x_gap(0)
           .add_colors([p.blue()])
           )
    bar.draw()
    bar.ax.set_title("Effect of TPM on number of positive control genes")
    bar.ax.axhline(2118, ls="--", color="k")
    bar.ax.set_ylim(0, 2300)
    bar.show()


def get_gene_categories(genes: list):
    go = biomart_go()
    go = go[go[BIOMART_GO_DOMAIN] == GO_PROCESS]
    go = (go[go[BIOMART_GENE_NAME].isin(genes)]
          .groupby(BIOMART_GO_NAME)
          .agg({BIOMART_GENE_NAME: lambda x: len(set(y for y in x))})
          .rename(columns={BIOMART_GENE_NAME: "count"})
          .reset_index())
    go = go.sort_values("count", ascending=False).reset_index()
    return go


def plot_gene_categories(genes: list):
    p = Palette()
    g = get_gene_categories(genes).head(30)
    bar = BarPlot(g["count"].values)
    (bar.change_orientation("y")
     .invert_y()
     .add_y_ticklabels(g[BIOMART_GO_NAME].values)
     .add_x_padding(0, 10)
     .add_x_label("Frequency")
     .add_colors([p.teal()])
     )
    bar.show()


def get_genes_from_categories(genes: list):
    go = biomart_go()
    go = go[go[BIOMART_GO_DOMAIN] == GO_PROCESS]
    go = (go[go[BIOMART_GENE_NAME].isin(genes)]
          .groupby(BIOMART_GO_NAME)
          .agg({BIOMART_GENE_NAME: lambda x: list(set(y for y in x))})
          .reset_index())
    go["count"] = go[BIOMART_GENE_NAME].map(lambda x: len(x))
    go = go.sort_values("count", ascending=False).reset_index()
    go = go[[BIOMART_GO_NAME, "count", BIOMART_GENE_NAME]]
    pd.set_option("display.max_colwidth", -1)
    # # print(go.head(30).to_string())
    # go = (go[(go[BIOMART_GO_NAME].str.contains("heart")) |
    #          (go[BIOMART_GO_NAME].str.contains("cardi"))])
    go = go.head(30)
    atlas = atlas_data()
    all_names = []
    for genes_list in go[BIOMART_GENE_NAME].values:
        temp = atlas[atlas[EXP_ATLAS_GENE_NAME].isin(genes_list)]
        del temp[EXP_ATLAS_GENE_ID]
        temp = temp.set_index(EXP_ATLAS_GENE_NAME)[[24, 48, 72]]
        temp["avg"] = temp.mean(axis=1)
        temp = temp.sort_values("avg", ascending=False).reset_index()
        number = 0
        k = True

        while k:
            try:
                g = temp[EXP_ATLAS_GENE_NAME].values[number]
            except IndexError:
                break
            if g not in all_names:
                all_names.append(g)
                k = False
            else:
                number += 1

    for f, n in zip(go[BIOMART_GO_NAME], all_names):
        print(f"{f}\t{n}")

    print(all_names)


def run():
    g = get_not_expressed(23, 73, 3, True)
    get_genes_from_categories(g)
