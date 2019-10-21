#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
#  TPM analysis related functions

import matplotlib.pylab as plt
import numpy as np
from SecretColors import Palette

from constants.boolean import *
from constants.ngs import *
from constants.system import *
from helpers.ngs_parser import get_string_tie_data

RANDOM_SEED = 1610
np.random.seed(RANDOM_SEED)


class TPMControl:
    def __init__(self, positive, negative):
        self._raw_positive = positive
        self._raw_negative = negative
        self._positive = None
        self._negative = None
        self._positive_test = None
        self._negative_test = None
        self._all_genes = None
        self.no_of_reserved_pos = 10
        self.no_of_reserved_neg = 6

    def _assign(self):

        self._positive = (np.random.choice(self._raw_positive,
                                           len(self._raw_positive) -
                                           self.no_of_reserved_pos,
                                           replace=False))
        self._negative = (np.random.choice(self._raw_negative,
                                           len(self._raw_negative) -
                                           self.no_of_reserved_neg,
                                           replace=False))

        self._positive_test = [x for x in self._raw_positive if x not in
                               self._positive]

        self._negative_test = [x for x in self._raw_negative if x not in
                               self._negative]

        x = [x for x in self._positive]
        y = [x for x in self._negative]
        x.extend(y)
        self._all_genes = x

    @property
    def positive(self):
        if self._positive is None:
            self._assign()
        return self._positive

    @property
    def negative(self):
        if self._negative is None:
            self._assign()
        return self._negative

    @property
    def positive_test(self):
        if self._positive_test is None:
            self._assign()
        return self._positive_test

    @property
    def negative_test(self):
        if self._negative_test is None:
            self._assign()
        return self._negative_test

    @property
    def all_genes(self):
        if self._all_genes is None:
            self._assign()
        return self._all_genes

    def error(self, pos, neg):
        e = 0
        for p in pos:
            if p not in self.positive:
                e += 1

        for n in neg:
            if n not in self.negative:
                e += 1

        return e


def get_data(time):
    return get_string_tie_data(GENOTYPE_WT, time, BIO_PROJECT_WINATA_LAB)[0]


def get_controls() -> TPMControl:
    negative1 = NEGATIVE_SET
    negative2 = NEGATIVE_CONTROL
    positive1 = POSITIVE_SET
    positive2 = POSITIVE_CONTROL

    negative2.extend(negative1)
    positive2.extend(positive1)

    negative = list(set(negative2))
    positive = list(set(positive2))

    return TPMControl(positive, negative)


def calculate_error(data, control, current_tpm):
    data = data[data[STRING_GENE_NAME].isin(control.all_genes)].reset_index(
        drop=True)
    data.loc[data[STRING_TPM] < current_tpm, STRING_TPM] = 0
    data.loc[data[STRING_TPM] > current_tpm, STRING_TPM] = 1
    positives = data[data[STRING_TPM] == 1][STRING_GENE_NAME].values
    negatives = data[data[STRING_TPM] == 0][STRING_GENE_NAME].values
    return control.error(positives, negatives)


def tpm_trace():
    p = Palette(show_warning=False)
    c = get_controls()
    tpm_range = np.linspace(0, 40, 300)
    colors = [p.red(), p.blue(), p.yellow()]
    all_error = []
    for t, col in zip([24, 48, 72], colors):
        errors = []
        data = get_data(t)
        for tpm in tpm_range:
            errors.append(calculate_error(data, c, tpm))
        all_error.append(errors)
        plt.plot(tpm_range, errors, color=col, label=f"{t} hpf", alpha=0.7)

    all_error = np.asarray(all_error)
    all_error = sum(all_error) / 3
    plt.plot(tpm_range, all_error, color=p.black(), label="Mean Error",
             ls="--", lw=1.5)
    plt.ylabel("Error")
    plt.xlabel("TPM")
    plt.axvline(min(all_error), color=p.gray())
    plt.legend(loc=0)
    plt.show()


def tpm_stat():
    p = Palette(show_warning=False)
    c = get_controls()
    tpm_range = np.linspace(0, 40, 300)
    colors = [p.red(), p.blue(), p.yellow()]
    min_tpm = []
    for _ in range(1):
        all_error = []
        for t, col in zip([24, 48, 72], colors):
            errors = []
            data = get_data(t)
            for tpm in tpm_range:
                errors.append(calculate_error(data, c, tpm))
            all_error.append(errors)

        all_error = np.asarray(all_error)
        all_error = sum(all_error) / 3
        min_tpm.append(min(all_error))

    plt.hist(min_tpm, color=p.green())
    plt.ylabel("Frequency")
    plt.xlabel("Minimum Average TPM")
    plt.show()


def run():
    tpm_stat()
