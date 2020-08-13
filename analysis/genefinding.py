#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#

import matplotlib.pyplot as plt
import pandas as pd
from SecretColors import Palette
from scipy.stats import gmean
from constants.biomart import BIOMART_GENE_ID, BIOMART_GENE_NAME
from constants.boolean import *
from constants.outputs import *
from helper.resolver import NameResolver

palette = Palette(show_warning=False)


def convert_id_to_gene(nr: NameResolver, data: pd.DataFrame):
    # Get the gene id to name
    name_data = pd.read_table(nr.id_to_gene).set_index(BIOMART_GENE_ID)
    name_data = name_data[BIOMART_GENE_NAME].to_dict()
    data["gene_id"] = data["gene_id"].apply(lambda x: name_data[x])
    return data


def count_genes(nr: NameResolver):
    lab = "winata"
    method = "star"
    condition = [24, 72]

    def __assign(x):
        if x[DESEQ2_PADJ] > 0.05:
            return palette.gray()
        if -1 < x[DESEQ2_LOG2_CHANGE] < 1:
            return palette.blue()
        else:
            return palette.red()

    # genes = list(set(INTERESTED_GENES).difference(BASE_GENES))
    genes = HOUSE_KEEPING
    fn = nr.deseq2_results(method=method, lab=lab, condition=condition,
                           lfc=True)
    data = convert_id_to_gene(nr, pd.read_csv(fn))
    data = data[data["gene_id"].isin(genes)]

    data["color"] = data.apply(lambda x: __assign(x), axis=1)
    data = data.sort_values(by="gene_id")

    ind = range(len(data["gene_id"].values))
    plt.barh(ind, data[DESEQ2_LOG2_CHANGE],
             color=data["color"],
             xerr=data[DESEQ2_LFCSE],
             error_kw=dict(ecolor=palette.gray(shade=70), capsize=5),
             zorder=100)
    plt.yticks(ind, data["gene_id"])
    plt.axvspan(-1, 1, color=palette.gray(shade=20), zorder=0)
    plt.axvline(1, ls="--", color=palette.black(), zorder=0)
    plt.axvline(-1, ls="--", color=palette.black(), zorder=0)
    plt.axvline(0, color=palette.black(), zorder=100)
    plt.grid(axis="both", ls=":", color=palette.gray(), zorder=0)
    plt.title(f"{condition[0]} vs {condition[1]} ({method}, {lab} lab)")
    plt.xlabel("Log$_2$ Fold Change")
    plt.ylabel("Gene")
    plt.tight_layout()
    plt.savefig("genes.png", dpi=300)
    plt.show()


def search_constant_genes(nr: NameResolver, lab, method, conditions):
    # lab = "winata"
    # method = "star"
    # conditions = [[24, 72]]
    dfs = []
    for con in conditions:
        fn = nr.deseq2_results(method=method, lab=lab, condition=con, lfc=True)
        d = pd.read_csv(fn)
        d = d[d[DESEQ2_PADJ] <= 0.05]
        d = d[d[DESEQ2_LOG2_CHANGE] > -1]
        d = d[d[DESEQ2_LOG2_CHANGE] < 1]
        d = d[["gene_id", DESEQ2_LOG2_CHANGE]]
        d = d.set_index("gene_id")
        d = d.rename(columns={DESEQ2_LOG2_CHANGE: f"{con}"})
        dfs.append(d)

    df = pd.concat(dfs, axis=1, sort=False, join="inner").reset_index()
    df = convert_id_to_gene(nr, df).set_index("gene_id")
    df = df.apply(lambda x: pow(x, 2))
    df["sum"] = df.sum(axis=1)
    df = df.sort_values(by="sum")
    needed_genes = df.index.values
    return needed_genes


def plot_fold_change(nr: NameResolver, genes: list):
    p = Palette()
    lab = "winata"
    method = "star"
    con = [24, 72]
    fn = nr.deseq2_results(method=method, lab=lab, condition=con, lfc=True)
    d = pd.read_csv(fn)
    d = convert_id_to_gene(nr, d)
    d = d[d['gene_id'].isin(genes)].fillna(1989)
    d = d[d[DESEQ2_PADJ] <= 0.05]
    d = d[d[DESEQ2_LOG2_CHANGE] > -1]
    d = d[d[DESEQ2_LOG2_CHANGE] < 1]
    plt.hist(d['log2FoldChange'].values, bins=50, color=p.ultramarine())
    plt.ylabel("Frequency")
    plt.xlabel("Log$_2$Fold Change")
    plt.title(f"Genes with Log$_2$Fold change (p<0.05)\n{lab} {method} {con}")
    plt.savefig("genes.png", dpi=300)
    plt.show()


def _extract_runs(lab: str, time: int, genotype: str) -> list:
    df = pd.read_csv("samples.csv")
    if lab not in ["winata", "yost"]:
        raise ValueError("Currently only 'winata' and 'yost' labs are "
                         "supported for this analysis")
    project = "PRJNA492280"
    if lab == "yost":
        project = "PRJNA407368"
    if genotype in ["wt", "wildtype"]:
        genotype = "wild type"

    df = df[df['BioProject'] == project]
    df = df[df['genetic_background'] == genotype]
    df = df[df['time'] == time]
    runs = df['Run'].values
    if len(runs) == 0:
        raise KeyError(f"Unable to find run for given combination of lab:"
                       f"{lab}, time:{time}, genotype:{genotype}")
    return runs


def get_average_dataframe(nr: NameResolver, lab: str, time: int,
                          genotype: str) -> pd.DataFrame:
    method = "stringtie"
    dfs = []
    for r in _extract_runs(lab, time, genotype):
        df = pd.read_csv(nr.run_output_file(r, method), sep="\t")
        df = df.sort_values(by="TPM", ascending=False)
        df = df.drop_duplicates(subset="Gene ID").reset_index(drop=True)
        df = df.set_index("Gene ID")
        dfs.append(df)

    tps = [tp[["TPM"]] for tp in dfs]
    tps = pd.concat(tps, join="inner", axis=1)
    name = "TPM_AVG"
    tps[name] = tps.mean(axis=1)
    tps = tps[[name]]
    dfs = dfs[0]
    del dfs["TPM"]
    dfs = pd.concat([dfs, tps], join="inner", axis=1)
    dfs = dfs.rename(columns={name: "TPM"}).reset_index()
    return dfs


def calculate_geometric_mean(nr: NameResolver):
    lab = "winata"
    genotype = "wt"
    cons = [24, 72]
    w = search_constant_genes(nr, lab, "star", [cons])
    dfs = {}
    for c in cons:
        df = get_average_dataframe(nr, lab, c, genotype)
        df = df[df["Gene Name"].isin(w)].reset_index(drop=True)
        df = df.sort_values(by="TPM", ascending=False).drop_duplicates(
            subset="Gene Name")
        df = df[["Gene Name", "TPM"]].set_index("Gene Name")
        df = df[df["TPM"] > 0]
        print(f"{c} Genes : {len(df)}")
        dfs[c] = df

    for c in dfs:
        values = dfs[c]["TPM"].values
        mean = round(gmean(values), 3)
        print(f"{c} G-Mean: {mean}")


def check_gene_expression(nr: NameResolver, genes: list):
    lab = "winata"
    genotype = "wt"
    time = 72
    # g_mean = 9.216
    g_mean = 6.547
    df = get_average_dataframe(nr, lab, time, genotype)
    df = df[df["Gene Name"].isin(genes)][["Gene Name", "TPM"]].reset_index(
        drop=True)
    df["NORM TMP"] = round(df["TPM"] / g_mean, 2)
    df["TPM"] = round(df["TPM"], 2)
    print(df)


def run():
    nr = NameResolver("config.json")
    genes = set(INTERESTED_GENES) - set(BASE_GENES)
    check_gene_expression(nr, genes)
    # print(palette.red(shade=30))
