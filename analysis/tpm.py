#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
#  TPM analysis related functions

from pprint import pprint as pp

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SecretColors import Palette

from constants.boolean import *
from constants.ngs import *
from constants.system import GENOTYPE_WT
from helpers.ngs_parser import average_data, get_string_tie_data

# RANDOM_SEED = 1610
# np.random.seed(RANDOM_SEED)
pd.options.mode.chained_assignment = None  # This needed to suppress warning


def get_positive_control():
    genes = POSITIVE_SET
    genes.extend(POSITIVE_CONTROL)
    genes = list(set(genes))
    np.random.shuffle(genes)
    return genes


def get_negative_control():
    genes = NEGATIVE_SET
    genes.extend(NEGATIVE_CONTROL)
    genes = list(set(genes))
    np.random.shuffle(genes)
    return genes


def get_data(hour):
    data = get_string_tie_data(GENOTYPE_WT, hour, BIO_PROJECT_WINATA_LAB)
    return average_data(data, OUTPUT_STRING_TIE)


class TPMAnalysis:
    def __init__(self):
        self._raw_positive = get_positive_control()
        self._raw_negative = get_negative_control()

        self.test_percentage = 25
        self.iterations = 100

        self._positive_control = None
        self._positive_test = None
        self._negative_control = None
        self._negative_test = None

    def _distribute_genes(self):
        self._positive_test = self._raw_positive[:int(len(
            self._raw_positive) * self.test_percentage / 100)]
        self._positive_control = self._raw_positive[len(
            self._positive_test):]

        self._negative_test = self._raw_negative[:int(len(
            self._raw_negative) * self.test_percentage / 100)]
        self._negative_control = self._raw_negative[len(self._negative_test):]

    @property
    def positive_control(self):
        if self._positive_control is None:
            self._distribute_genes()
        return self._positive_control

    @property
    def positive_test(self):
        if self._positive_test is None:
            self._distribute_genes()
        return self._positive_test

    @property
    def negative_control(self):
        if self._negative_control is None:
            self._distribute_genes()
        return self._negative_control

    @property
    def negative_test(self):
        if self._negative_test is None:
            self._distribute_genes()
        return self._negative_test

    @staticmethod
    def calculate_error(data, tpm_threshold, pos_genes, neg_genes):
        pos = data[data[STRING_GENE_NAME].isin(pos_genes)]
        neg = data[data[STRING_GENE_NAME].isin(neg_genes)]

        pos[pos[STRING_TPM] >= tpm_threshold] = 1
        neg[neg[STRING_TPM] >= tpm_threshold] = 1
        pos[pos[STRING_TPM] < tpm_threshold] = 0
        neg[neg[STRING_TPM] < tpm_threshold] = 0

        pos_err = len(pos) - pos[STRING_TPM].sum()
        neg_err = neg[STRING_TPM].sum()
        return tpm_threshold, pos_err + neg_err, pos_err, neg_err

    def _calculate_rates(self, pos_err, neg_err):
        true_positive = len(self.positive_test) - pos_err
        true_negative = len(self.negative_test) - neg_err
        false_positive = neg_err
        false_negative = pos_err

        tpr = true_positive / (true_positive + false_negative)
        fpr = 1 - (true_negative / (true_negative + false_positive))
        return tpr, fpr

    def analyze(self, data):
        d = data[[STRING_GENE_NAME, STRING_TPM]]
        errs = []
        for t in np.linspace(0, 10, 500):
            errs.append(self.calculate_error(d, t,
                                             self.positive_control,
                                             self.negative_control))

        best_tpm = sorted(errs, key=lambda x: x[1])[0][0]

        e = self.calculate_error(d, best_tpm,
                                 self.positive_test,
                                 self.negative_test)

        rates = self._calculate_rates(e[2], e[3])
        distance = np.sqrt(pow(1 - rates[0], 2) + pow(rates[1], 2))

        return best_tpm, distance, rates


def plot_roc_curve(hour, iterations=100):
    p = Palette()
    data = get_data(hour)
    values = []
    for _ in range(iterations):
        t = TPMAnalysis()
        values.append(t.analyze(data))

    values = sorted(values, key=lambda m: m[1])
    pp(values)

    y = [x1[2][0] for x1 in values]
    x = [x1[2][1] for x1 in values]

    plt.scatter(x, y, marker="o", color=p.red())
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.plot([0, 1], [0, 1], ls="--", color=p.gray())
    plt.annotate(f"TPM={round(values[0][0], 2)}",
                 xy=(x[0], y[0]),
                 xytext=(x[0] + 0.1, y[0] - 0.1),
                 arrowprops=dict(arrowstyle="->"))
    plt.ylabel("True Positive Rate")
    plt.xlabel("False Positive Rate")
    plt.title(f"ROC curve at {hour} hpf")
    plt.tight_layout()
    plt.savefig("roc.png", type="png", dpi=300)
    plt.show()


def run():
    plot_roc_curve(24)
