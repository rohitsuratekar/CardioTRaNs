#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
# All Wrappers

import os
import warnings

import pandas as pd
import yaml

from constants.outputs import *
from constants.zfin import *


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


class ZData:
    WHOLE_ORGANISM = "ZFA:0001094"

    def __init__(self, expression: str, xpat: str, ontology: str):
        self._expression = expression
        self._xpat = xpat
        self._ontology = ontology
        for filename in [expression, xpat, ontology]:
            if not os.path.isfile(filename):
                raise FileNotFoundError(f"ZFIN data not found at {filename}")

        self._exp_df = None
        self._xpat_df = None
        self._onto_df = None
        self.start_stage = None
        self.end_stage = None
        self.structures = []
        self.expressed = True

        self._wt_lines = None
        self._fish_xpat = None

    @property
    def expression(self):
        if self._exp_df is None:
            self._exp_df = pd.read_csv(self._expression, sep="\t")
        return self._exp_df

    @property
    def xpat(self):
        if self._xpat_df is None:
            self._xpat_df = pd.read_csv(self._xpat, sep="\t")
        return self._xpat_df

    @property
    def ontology(self):
        if self._onto_df is None:
            self._onto_df = pd.read_csv(self._ontology, sep="\t")
        return self._onto_df

    def _onto_mapping(self, key, value) -> dict:
        return dict(
            zip(self.ontology[key], self.ontology[value]))

    def clear_filters(self):
        self.start_stage = None
        self.end_stage = None
        self.structures = []

    def _apply_filters(self):
        tmp = "Temp"
        df = self.expression.copy()
        df = df[[ZFIN_EXP_GENE_SYMBOL,
                 ZFIN_EXP_START_STAGE,
                 ZFIN_EXP_END_STAGE,
                 ZFIN_EXP_SUPER_STR_NAME]]
        if self.start_stage:
            start_map = self._onto_mapping(ZFIN_ONT_STAGE_NAME,
                                           ZFIN_ONT_BEGIN_HOUR)
            df[tmp] = df[ZFIN_EXP_START_STAGE].map(lambda x: start_map[x])
            df = df[df[tmp] >= self.start_stage]
            del df[tmp]

        if self.end_stage:
            end_map = self._onto_mapping(ZFIN_ONT_STAGE_NAME,
                                         ZFIN_ONT_END_HOUR)
            df[tmp] = df[ZFIN_EXP_END_STAGE].map(lambda x: end_map[x])
            df = df[df[tmp] <= self.end_stage]
            del df[tmp]

        if self.structures:
            search_str = "|".join(self.structures).lower().strip()
            df = df[df[ZFIN_EXP_SUPER_STR_NAME].str.contains(search_str)]

        df = df.reset_index(drop=True)
        return df

    def _apply_not_expressed(self):
        tmp = "Temp"
        df = self.xpat.copy()

        # For strict implementation, check expression is not found in whole
        # organism. This will also override all structure filters.
        df = df[df[ZFIN_XPAT_ANATOMY_SUPER] == ZData.WHOLE_ORGANISM]
        df = df[df[ZFIN_XPAT_EXPRESSION_FOUND] == "f"]

        if self.start_stage:
            start_map = self._onto_mapping(ZFIN_ONT_STAGE_ID,
                                           ZFIN_ONT_BEGIN_HOUR)
            df[tmp] = df[ZFIN_XPAT_START_STAGE].map(lambda x: start_map[x])

            df = df[df[tmp] >= self.start_stage]
            del df[tmp]
        if self.end_stage:
            end_map = self._onto_mapping(ZFIN_ONT_STAGE_ID, ZFIN_ONT_END_HOUR)
            df[tmp] = df[ZFIN_XPAT_END_STAGE].map(lambda x: end_map[x])
            df = df[df[tmp] <= self.end_stage]
            del df[tmp]

        if self.structures:
            warnings.warn(F"Structure filter not applied. For strict "
                          F"implementation of Negative genes, "
                          F"whole organism is considered instead structures")

        # Extract the expression IDs and now take a file which will have
        # both expression ID and gene IDs. Filter the appropriate rows
        df = df[ZFIN_XPAT_EXPRESSION_ID].values
        xpat = pd.read_csv(self._fish_xpat, sep="\t")
        xpat = xpat[xpat[ZFIN_ASSAY_EXPRESSION_ID].isin(df)]

        # Only takes expressions IDs with wild types fish
        fish = pd.read_csv(self._wt_lines, sep="\t")
        xpat = xpat[xpat[ZFIN_ASSAY_FISH_ID].isin(fish[ZFIN_FISH_ID].values)]
        xpat = xpat.reset_index(drop=True)

        # Genes which have shown not expressing at given points
        genes = xpat[ZFIN_ASSAY_GENE_SYMBOL].values
        print(genes)

    @property
    def dataframe(self):
        if self.expressed:
            return self._apply_filters()
        else:
            return self._apply_not_expressed()

    def add_start_stage(self, time: float):
        self.start_stage = time

    def add_end_stage(self, time: float):
        self.end_stage = time

    def add_structures(self, struct):
        if isinstance(struct, str):
            self.structures.append(struct)
        else:
            self.structures.extend(struct)

    def search_non_expressed(self, fish_xpat: str, wt_lines: str):
        self._fish_xpat = fish_xpat
        self._wt_lines = wt_lines
        self.expressed = False
