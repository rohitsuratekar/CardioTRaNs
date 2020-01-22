#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
# File and Data manager

import json

import pandas as pd

from constants.extra import *
from helper.parsers import ConfigParser
from helper.utils import exists_path


class FileManager:
    def __init__(self, config: ConfigParser):
        self.config = config
        self.log = config.log

    def _validate_path(self, path):
        if not exists_path(path):
            self.log.error(f"File not found at {path}. Please check your "
                           f"configuration file {self.config.filename} ")

    def _get_sra_data(self) -> pd.DataFrame:
        self._validate_path(self.config.sra_meta)
        return pd.read_csv(self.config.sra_meta)

    def _sra_out_file(self, sra: str, srr: str, method: str):
        method = method.strip().lower()
        if method == "stringtie":
            return f"{sra}/stringtie/{srr}_gene_expression.tsv"
        elif method == "salmon":
            return f"{sra}/salmon/quant.sf"
        elif method == "kallisto":
            return f"{sra}/kallisto/abundance.tsv"
        elif method == "star":
            return f"{sra}/star/{srr}_Aligned.sortedByCoord.out.bam"
        else:
            self.log.error(f"Method {method} is not supported by this "
                           f"pipeline yet.")

    def _get_runs(self, sra):
        meta = f"{self.config.sra_folder}/{sra}/{sra}_meta.json"
        if not exists_path(meta):
            self.log.error(f"Automatic retrieval of metadata file of {sra} "
                           f"failed. Make sure if you have used latest "
                           f"version of CardioPipeLine to generate this "
                           f"metadata. If you want to use your own input "
                           f"file please pass full path to 'filename' "
                           f"argument. E.g fm.sra_output(id, salmon, "
                           f"filename=path/to/file  )")
        with open(meta) as f:
            data = json.load(f)
            return data["Runs"]

    def salmon_output(self, sra_id: str, *, filename: str = None,
                      srr: str = None) -> pd.DataFrame:
        return self.method_output(sra_id, method="salmon", filename=filename,
                                  srr=srr)

    def kallisto_output(self, sra_id: str, *, filename: str = None,
                        srr: str = None) -> pd.DataFrame:
        return self.method_output(sra_id, method="kallisto", filename=filename,
                                  srr=srr)

    def stringtie_output(self, sra_id: str, *, filename: str = None,
                         srr: str = None) -> pd.DataFrame:
        return self.method_output(sra_id, method="stringtie",
                                  filename=filename,
                                  srr=srr)

    def method_output(self, sra_id: str,
                      *,
                      method: str,
                      filename: str = None,
                      srr: str = None) -> pd.DataFrame:

        # Sanity check
        mt = method.strip().lower()

        if filename is None:
            runs = self._get_runs(sra_id)
            if len(runs) > 1:
                self.log.error(
                    f"More than 1 runs found for {sra_id}. Please pass run "
                    f"id with argument 'srr' to continue. E.G. "
                    f"fm.sra_output(id, salmon, srr=run_id)")
            elif len(runs) == 0:
                self.log.error(f"No run information found for the sra id "
                               f"{sra_id}. Please use 'filename' parameter "
                               f"to point program towards output file.")

            if srr is None:
                srr = runs[0]

            filename = f"{self.config.sra_folder}" \
                       f"/{self._sra_out_file(sra_id, srr, method)}"

        if mt == "star":
            self.log.error("STAR outputs BAM file which we can not use in "
                           "this analysis. Do you mean StringTie?")

        self.log.info(
            f"Reading data processed with {method} from the file {filename}")

        return pd.read_table(filename)

    def _get_deseq2_file(self, method, lab, conditions, lfc):
        con = sorted(conditions)
        bf = f"{self.config.deseq2_folder}/{lab}/Analysis{con[0]}v{con[1]}"
        lfc_add = ""
        if lfc:
            lfc_add = "_lfc"
            self.log.info(f"Trying get data processed with LFC")
        return f"{bf}/{method.strip().lower()}.result{lfc_add}.csv"

    def _get_deseq2_count_file(self, method, lab, conditions):
        con = sorted(conditions)
        bf = f"{self.config.deseq2_folder}/{lab}/Analysis{con[0]}v{con[1]}"
        return f"{bf}/{method.strip().lower()}.vst.csv"

    def _sanity_check(self, method: str, *,
                      lab: str = None,
                      conditions: list = None,
                      filename: str = None):
        # Sanity check
        if method.strip().lower() not in ["star", "salmon", "kallisto",
                                          "stringtie"]:
            self.log.error(f"Method '{method}' is not supported by this "
                           f"pipeline yet")

        if filename is None:
            if lab is None or conditions is None:
                self.log.error("If 'filename' is not provided, 'lab' and "
                               "'conditions' argument needed to proceed.")

            if len(conditions) != 2:
                self.log.error("Currently this pipeline handles only 2 "
                               "conditions at a time")

    def deseq2_counts(self, method: str, *, lab: str = None,
                      conditions: list = None,
                      filename: str = None) -> pd.DataFrame:

        # Sanity check
        self._sanity_check(method, lab=lab,
                           conditions=conditions,
                           filename=filename)

        if filename is None:
            filename = self._get_deseq2_count_file(method, lab, conditions)

        if not exists_path(filename):
            self.log.error(
                f"File does not exist {filename}. "
                f"Please check if you have provided correct details.",
                exception=FileNotFoundError)

        self.log.info(f"Reading DESeq2 count data from the file {filename}")

        return pd.read_csv(filename)

    def deseq2_results(self, method: str, *,
                       lab: str = None,
                       lfc: bool = False,
                       conditions: list = None,
                       filename: str = None) -> pd.DataFrame:

        # Sanity check
        self._sanity_check(method, lab=lab,
                           conditions=conditions,
                           filename=filename)

        if filename is None:
            filename = self._get_deseq2_file(method, lab, conditions, lfc)

        if not exists_path(filename):
            self.log.error(
                f"File does not exist {filename}. "
                f"Please check if you have provided correct details.",
                exception=FileNotFoundError)

        self.log.info(f"Reading DESeq2 data from the file {filename}")

        return pd.read_csv(filename)

    def _get_deseq2_data(self) -> pd.DataFrame:
        return pd.read_csv(self.config.deseq2_meta)

    def extract_sra_ids(self, time: int, *,
                        genotype: str = "wt",
                        lab: str = None,
                        bio_project: str = None) -> list:

        if lab is None and bio_project is None:
            self.log.error("Please provide either lab name or BioProject ID")

        lab_keys = {
            "winata": "PRJNA492280",
            "hills": "PRJNA407368"
        }

        if bio_project is None:
            if not lab.strip().lower() in lab_keys.keys():
                self.log.error(f"Unknown lab {lab}. Please provide "
                               f"BioProject ID instead. Currently only "
                               f"following labs are available for direct "
                               f"use: {list(lab_keys.keys())}")

            bio_project = lab_keys[lab.strip().lower()]

        data = self._get_sra_data()
        value = data[data[SRA_PROJECT] == bio_project]
        value = value[value[SRA_GENOTYPE] == genotype]
        value = value[value[SRA_TIME] == time]
        return list(value[SRA_ID].values)
