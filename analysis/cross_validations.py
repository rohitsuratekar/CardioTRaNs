#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Cross Validation of the TPM values
import numpy as np
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import KFold
from helpers.ngs_parser import *
from constants.boolean import *
import pandas as pd
import matplotlib.pyplot as plt

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
    def __init__(self):
        self.positive_control = get_positive_control()
        self.negative_control = get_negative_control()
        self._all_genes = None

    @property
    def all_genes(self):
        if self._all_genes is None:
            self._all_genes = np.concatenate(
                (self.positive_control, self.negative_control))
        return self._all_genes

    @staticmethod
    def gene(output_method):
        if output_method == OUTPUT_STRING_TIE:
            return STRING_GENE_NAME
        elif output_method == OUTPUT_SALMON:
            return SALMON_GENE_NAME

    @staticmethod
    def tpm(output_method):
        if output_method == OUTPUT_STRING_TIE:
            return STRING_TPM
        elif output_method == OUTPUT_SALMON:
            return SALMON_TPM

    def get_binary_label(self, genes):
        return [0 if x in self.negative_control else 1 for x in genes]

    def tpm_array(self, data, om, threshold):
        d = data[[self.gene(om), self.tpm(om)]]
        d.loc[d[self.tpm(om)] >= threshold, self.tpm(om)] = 1
        d.loc[d[self.tpm(om)] < threshold, self.tpm(om)] = 0
        return d

    def predictions(self, data, genes, om):
        d = pd.Series(data[self.tpm(om)].values, index=data[self.gene(
            om)]).to_dict()
        return [d[x] for x in genes]

    def analyse(self, data, om):
        kfold = KFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
        # Split value of 5 is used so that we will get 20% in our test dataset
        for training, test in kfold.split(self.all_genes):
            binary_labels = self.get_binary_label(self.all_genes[training])
            tpms = self.tpm_array(data, om, 1)
            hypothesis = self.predictions(data[[self.gene(om), self.tpm(om)]],
                                          self.all_genes[training], om)
            fpr, tpr, _ = roc_curve(binary_labels, hypothesis, pos_label=1)
            plt.plot(fpr, tpr)
            print(_[np.argmax(tpr - fpr)])

        plt.show()


def run():
    data = get_data(24, OUTPUT_STRING_TIE)
    t = TPMAnalysis()
    t.analyse(data, OUTPUT_STRING_TIE)
