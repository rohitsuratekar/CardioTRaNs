#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
# All Wrappers

import pandas as pd
from constants.outputs import *
import yaml


class OutRNA:
    GENE_NAME = "gene_name"
    TPM = "TPM"
    AVG_TPM = "AVG_TPM"

    def __init__(self, filename: str, method: str):
        self.filename = filename
        self.method = method.strip().lower()
        available = ["salmon", "stringtie", "kallisto", "deseq2"]

        if self.method not in available:
            raise KeyError(f"Method '{method}' is not available at the "
                           f"moment. Use any of these methods: {available}")

    def _id_mapping(self) -> dict:
        with open("config.yml") as f:
            data = yaml.load(f, Loader=yaml.SafeLoader)
        df = pd.read_csv(data["biomart"], sep="\t")
        if self.method == "deseq2":
            df = dict(zip(df[BIOMART_GENE_ID], df[BIOMART_GENE_NAME]))
        else:
            df = dict(zip(df[BIOMART_TRANSCRIPT_ID_VERSION],
                          df[BIOMART_GENE_NAME]))
        return df

    def _parse_salmon(self):
        df = pd.read_csv(self.filename, sep="\t")
        mapping = self._id_mapping()
        df[OutRNA.GENE_NAME] = df[SALMON_NAME].map(lambda x: mapping[x])
        df = df.sort_values(by=SALMON_TPM, ascending=False)
        df = df.drop_duplicates(subset=OutRNA.GENE_NAME).reset_index(drop=True)
        return df

    def _parse_stringtie(self):
        df = pd.read_csv(self.filename, sep="\t")
        df = df.sort_values(by=STRINGTIE_TPM, ascending=False)
        df = df.drop_duplicates(subset=STRINGTIE_NAME).reset_index(drop=True)
        df[OutRNA.GENE_NAME] = df[STRINGTIE_NAME]
        return df

    def _parse_kallisto(self):
        df = pd.read_csv(self.filename, sep="\t")
        mapping = self._id_mapping()
        df[OutRNA.GENE_NAME] = df[KALLISTO_TARGET_ID].map(lambda x: mapping[x])
        df = df.sort_values(by=KALLISTO_TPM, ascending=False)
        df = df.drop_duplicates(subset=OutRNA.GENE_NAME).reset_index(drop=True)
        df = df.rename(columns={"tpm": "TPM"})
        return df

    def _parse_deseq2(self):
        df = pd.read_csv(self.filename)
        mapping = self._id_mapping()
        df[OutRNA.GENE_NAME] = df[DESEQ2_GENE_ID].map(lambda x: mapping[x])
        return df

    @property
    def dataframe(self) -> pd.DataFrame:
        if self.method == "salmon":
            return self._parse_salmon()
        elif self.method == "stringtie":
            return self._parse_stringtie()
        elif self.method == "kallisto":
            return self._parse_kallisto()
        elif self.method == "deseq2":
            return self._parse_deseq2()
        else:
            raise ValueError(f"Method {self.method} is not implemented yet")
