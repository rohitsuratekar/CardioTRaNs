#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Cross Validation of the TPM values
import matplotlib.pyplot as plt
import numpy as np
from SecretColors import Palette
from sklearn.model_selection import LeavePOut

from constants.boolean import *
from helpers.ngs_parser import *
from pprint import pprint as pp

RANDOM_STATE = 1989
pd.options.mode.chained_assignment = None  # This needed to suppress warning


def get_data(hour, output_method):
    if output_method == OUTPUT_STRING_TIE:
        data = get_string_tie_data(GENOTYPE_WT, hour, BIO_PROJECT_WINATA_LAB)
    elif output_method == OUTPUT_SALMON:
        data = get_salmon_data(GENOTYPE_WT, hour, BIO_PROJECT_WINATA_LAB)
    else:
        raise Exception("Invalid output method")
    return average_data(data, output_method)


def get_positive_control():
    genes = POSITIVE_SET
    genes.extend(POSITIVE_CONTROL)
    return np.array(list(set(genes)))


def get_negative_control():
    genes = NEGATIVE_SET
    genes.extend(NEGATIVE_CONTROL)
    return np.array(list(set(genes)))


class TPMAnalysis:
    def __init__(self, output_method):
        self.positive_control = get_positive_control()
        self.negative_control = get_negative_control()
        self._all_genes = None
        self.output_method = output_method
        self.min_tpm = 0
        self.max_tpm = 10
        self.divisions = 100

    @property
    def all_genes(self):
        if self._all_genes is None:
            self._all_genes = np.concatenate(
                (self.positive_control, self.negative_control))
        return self._all_genes

    @property
    def col_gene(self):
        if self.output_method == OUTPUT_STRING_TIE:
            return STRING_GENE_NAME
        elif self.output_method == OUTPUT_SALMON:
            return SALMON_GENE_NAME

    @property
    def col_tpm(self):
        if self.output_method == OUTPUT_STRING_TIE:
            return STRING_TPM
        elif self.output_method == OUTPUT_SALMON:
            return SALMON_TPM

    @staticmethod
    def calculate_rates(genes, pos, neg):
        tp = len([x for x in pos if genes[x] == 1])
        fn = len([x for x in pos if genes[x] == 0])

        fp = len([x for x in neg if genes[x] == 1])
        tn = len([x for x in neg if genes[x] == 0])

        tpr = tp / (tp + fn)
        fpr = fp / (fp + tn)

        return fpr, tpr

    def split_genes(self, genes):
        pos = [x for x in genes if x in self.positive_control]
        neg = [x for x in genes if x in self.negative_control]
        return pos, neg

    def analyse(self, data, tpm_threshold, pos, neg):
        d = data[[self.col_gene, self.col_tpm]]
        d.loc[d[self.col_tpm] >= tpm_threshold, self.col_tpm] = 1
        d.loc[d[self.col_tpm] < tpm_threshold, self.col_tpm] = 0
        d = pd.Series(d[self.col_tpm].values, index=d[self.col_gene]).to_dict()
        return self.calculate_rates(d, pos, neg)

    def roc(self, hour, pos=None, neg=None):
        if pos is None:
            pos = self.positive_control
        if neg is None:
            neg = self.negative_control

        data = get_data(hour, self.output_method)
        fpr = []
        tpr = []
        thresholds = np.linspace(self.min_tpm, self.max_tpm, self.divisions)

        for tpm in thresholds:
            f, t = self.analyse(data, tpm, pos, neg)
            fpr.append(f)
            tpr.append(t)
        return np.array(fpr), np.array(tpr), thresholds


def plot_simple_roc(hour):
    t = TPMAnalysis(OUTPUT_STRING_TIE)
    p = Palette()
    fnr, tpr, thresholds = t.roc(hour)
    threshold_index = np.argmax(tpr - fnr)
    plt.plot([0, 1], [0, 1], ls="--", color=p.gray())

    plt.plot(fnr, tpr, marker="o", markerfacecolor=p.white(), color=p.red())
    plt.plot(fnr[threshold_index], tpr[threshold_index], marker="o",
             color=p.blue())

    plt.annotate(f"TPM={round(thresholds[threshold_index], 2)}",
                 xy=(fnr[threshold_index], tpr[threshold_index]),
                 xytext=(
                     fnr[threshold_index] + 0.1, tpr[threshold_index] - 0.1),
                 arrowprops=dict(arrowstyle="->"),
                 bbox=dict(boxstyle="round", fc=p.gray(shade=10)))

    plt.ylabel("True Positive Rate")
    plt.xlabel("False Positive Rate")
    plt.title(f"ROC curve at {hour} hpf")
    plt.grid(axis="both", zorder=10, alpha=0.6)
    plt.tight_layout()
    plt.savefig("roc.png", type="png", dpi=300)
    plt.show()


def plot_cross_validate_roc(hour):
    t = TPMAnalysis(OUTPUT_STRING_TIE)
    p = Palette()
    out_group = 3

    samples = LeavePOut(p=out_group)
    all_fpr = []
    all_tpr = []
    for training, _ in samples.split(t.all_genes):
        pos, neg = t.split_genes(t.all_genes[training])
        fpr, tpr, thresholds = t.roc(hour, pos, neg)
        all_fpr.append(fpr)
        all_tpr.append(tpr)

    mean_tpr = np.mean(all_tpr, axis=0)
    mean_fpr = np.mean(all_fpr, axis=0)
    std_tpr = np.std(all_tpr, axis=0)

    plt.fill_between(mean_fpr, mean_tpr + std_tpr, mean_tpr - std_tpr,
                     color=p.gray(), alpha=0.5)
    plt.plot(mean_fpr, mean_tpr, color=p.red())
    plt.ylabel("True Positive Rate")
    plt.xlabel("False Positive Rate")
    plt.title(f"ROC curve at {hour} hpf (out_group={out_group})")
    plt.grid(axis="both", zorder=10, alpha=0.6)
    plt.tight_layout()
    plt.savefig("roc.png", type="png", dpi=300)
    plt.show()


def run():
    plot_cross_validate_roc(72)
