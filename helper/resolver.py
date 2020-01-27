#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
# Name resolver class

import json
import logging

from helper.logging import Log
from helper.utils import exists_path


class NameResolver:
    """
    A Simple wrapper around the configuration files which can automatically
    generates the file names for downstream analysis
    """

    def __init__(self, config_file: str, log: Log = None):
        self.config_file = config_file
        self._log = log
        self._data = None

    @property
    def data(self):
        if self._data is None:
            with open(self.config_file) as f:
                self._data = json.load(f)
        return self._data

    @property
    def log(self) -> Log:
        debug = -1
        log_format = self.data["setup"]["log_format"]
        if not self.data["setup"]["show_debug_log"]:
            debug = logging.DEBUG
        if len(str(log_format).strip()) == 0:
            log_format = "%(message)s"
        if self._log is None:
            log = Log(show_log=self.data["setup"]["show_log"],
                      add_to_file=self.data["setup"]["save_log"],
                      filename=self.data["setup"]["log_file"],
                      logging_format=log_format,
                      min_log_level=debug)
            self._log = log
            self._log.info("Logging is initialized")
        return self._log

    def _validate(self, *args):
        data = self.data
        for key in args:
            try:
                data = data[key]
            except KeyError:
                self.log.error(
                    f"Configuration file {self.config_file} does not have "
                    f"data at node '{'/'.join(args)}'", exception=KeyError)

        return data

    @property
    def sra_folder(self) -> str:
        return self._validate("rna_seq", "data_folder")

    @property
    def sra_meta(self) -> str:
        return self._validate("rna_seq", "meta")

    @property
    def deseq2_meta(self) -> str:
        return self._validate("deseq2", "meta")

    @property
    def deseq2_folder(self) -> str:
        return self._validate("deseq2", "data_folder")

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

    def get_runs(self, sra):
        meta = f"{self.sra_folder}/{sra}/{sra}_meta.json"
        if not exists_path(meta):
            self.log.error(f"Automatic retrieval of metadata file {meta} "
                           f"failed. Make sure if you have used latest "
                           f"version of CardioPipeLine to generate this "
                           f"metadata.")
        with open(meta) as f:
            data = json.load(f)
            return data["Runs"]

    def mapping_output(self, method: str, sra_id: str, *, srr: str = None):
        mt = method.strip().lower()  # Sanity check
        runs = self.get_runs(sra_id)
        if len(runs) > 1:
            self.log.error(
                f"More than 1 runs found for '{sra_id}'. Please pass run "
                f"id with argument 'srr' to continue. E.G. "
                f"fm.mapping_output(id, salmon,SRA_ID, srr=run_id)")
        elif len(runs) == 0:
            self.log.error(f"No run information found for the sra id "
                           f"'{sra_id}' in automatically generated metadata.")

        if srr is None:
            srr = runs[0]

        return f"{self.sra_folder}/{self._sra_out_file(sra_id, srr, mt)}"

    def _get_deseq2_folder(self, lab: str, condition: list):
        if len(condition) != 2:
            self.log.error("There should be exactly two conditions. e.g. ["
                           "24,30]")

        if lab.strip().lower() not in ["winata", "hills"]:
            self.log.error("Currently, only following labs are available : "
                           "winata, hills")

        condition = sorted(condition)
        folder = f"Analysis{condition[0]}v{condition[1]}"
        return f"{self.deseq2_folder}/{lab.strip().lower()}/{folder}"

    def deseq2_counts(self, method, *, lab: str, condition: list):
        prefix = self._get_deseq2_folder(lab, condition)
        return f"{prefix}/{method}.counts.csv"

    def deseq2_results(self, method, *, lab: str, condition: list,
                       lfc: bool = False):
        prefix = self._get_deseq2_folder(lab, condition)
        extra = ""
        if lfc:
            extra = "_lfc"
        return f"{prefix}/{method}.result{extra}.csv"

    def deseq2_vst(self, method, *, lab: str, condition: list):
        prefix = self._get_deseq2_folder(lab, condition)
        return f"{prefix}/{method}.vst.csv"

    def deseq2_rlog(self, method, *, lab: str, condition: list):
        prefix = self._get_deseq2_folder(lab, condition)
        return f"{prefix}/{method}.rlog.csv"

    @property
    def id_to_gene(self):
        return self._validate("biomart", "id_gene_rel")
