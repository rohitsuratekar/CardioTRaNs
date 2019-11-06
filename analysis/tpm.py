#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
#  TPM analysis related functions

from pprint import pprint as pp

import matplotlib.pyplot as plt
import numpy as np
from SecretColors import Palette
from matplotlib.lines import Line2D

from constants.boolean import *
from helpers.ngs_parser import *

RANDOM_SEED = 1610
np.random.seed(RANDOM_SEED)
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


def get_data(hour, output_method):
    if output_method == OUTPUT_STRING_TIE:
        data = get_string_tie_data(GENOTYPE_WT, hour, BIO_PROJECT_WINATA_LAB)
    elif output_method == OUTPUT_SALMON:
        data = get_salmon_data(GENOTYPE_WT, hour, BIO_PROJECT_WINATA_LAB)
    else:
        raise Exception("Invalid output method")
    return average_data(data, output_method)


class TPMAnalysis:
    def __init__(self):
        self._raw_positive = get_positive_control()
        self._raw_negative = get_negative_control()

        self.test_percentage = 25
        self.divisions = 100
        self.output_method = OUTPUT_STRING_TIE

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

    @property
    def gene_column(self):
        if self.output_method == OUTPUT_STRING_TIE:
            return STRING_GENE_NAME
        elif self.output_method == OUTPUT_SALMON:
            return SALMON_GENE_NAME
        else:
            raise Exception(f"Unknown output method: {self.output_method}")

    @property
    def tpm_column(self):
        if self.output_method == OUTPUT_STRING_TIE:
            return STRING_TPM
        elif self.output_method == OUTPUT_SALMON:
            return SALMON_TPM
        else:
            raise Exception(f"Unknown output method: {self.output_method}")

    def calculate_error(self, data, tpm_threshold, pos_genes, neg_genes):
        pos = data[data[self.gene_column].isin(pos_genes)]
        neg = data[data[self.gene_column].isin(neg_genes)]

        pos[pos[self.tpm_column] >= tpm_threshold] = 1
        neg[neg[self.tpm_column] >= tpm_threshold] = 1
        pos[pos[self.tpm_column] < tpm_threshold] = 0
        neg[neg[self.tpm_column] < tpm_threshold] = 0

        pos_err = len(pos) - pos[self.tpm_column].sum()
        neg_err = neg[self.tpm_column].sum()
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
        d = data[[self.gene_column, self.tpm_column]]
        errs = []
        for t in np.linspace(1, 10, self.divisions):
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


def get_plot_values(data, iterations, divisions, output_method):
    values_string = []
    for _ in range(iterations):
        t = TPMAnalysis()
        t.divisions = divisions
        t.output_method = output_method
        values_string.append(t.analyze(data))

    values_string = sorted(values_string, key=lambda m: m[1])
    y = [x1[2][0] for x1 in values_string]
    x = [x1[2][1] for x1 in values_string]
    y.insert(0, 0)
    x.insert(0, 0)
    y.append(1)
    x.append(1)

    return x, y


def plot_roc_curve(hour, iterations=10, divisions=50):
    p = Palette()
    colors = [p.red(), p.blue()]
    outputs = [OUTPUT_STRING_TIE, OUTPUT_SALMON]
    for out, color in zip(outputs, colors):
        data = get_data(hour, out)
        x, y = get_plot_values(data, iterations, divisions, out)
        plt.plot(x, y, marker="o", markerfacecolor=p.white(),
                 color=color)

    plt.xlim(-0.1, 1.1)
    plt.ylim(0, 1.1)
    plt.plot([0, 1], [0, 1], ls="--", color=p.gray())
    # plt.annotate(f"TPM={round(values_string[0][0], 2)}",
    #              xy=(x[1], y[1]),
    #              xytext=(x[1] + 0.1, y[1] - 0.1),
    #              arrowprops=dict(arrowstyle="->"),
    #              bbox=dict(boxstyle="round", fc=p.gray(shade=10)))
    plt.ylabel("True Positive Rate")
    plt.xlabel("False Positive Rate")
    plt.title(f"ROC curve at {hour} hpf")
    plt.grid(axis="both", zorder=10, alpha=0.6)

    leg = [Line2D([0], [0], ls="--", label="Random", color=p.gray())]
    for label, color in zip(outputs, colors):
        leg.append(Line2D([0], [0], label=label, color=color))
    plt.legend(handles=leg, loc="lower right")
    plt.tight_layout()
    plt.savefig("roc.png", type="png", dpi=300)
    plt.show()


def run():
    plot_roc_curve(72)
