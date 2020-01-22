#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
# Various parsers will go here

import json
import logging

from helper.logging import Log
from helper.utils import exists_path


class ConfigParser:
    def __init__(self, filename: str, log: Log = None):
        self.filename = filename
        self._log = log
        if not exists_path(filename):
            raise FileNotFoundError(f"Configuration file {filename} not "
                                    f"found. Please provide full path.")
        with open(filename) as f:
            self.data = json.load(f)

        self.log.info("Configuration data successfully loaded")

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

    def _return_with_validation(self, *args):
        data = self.data
        for key in args:
            try:
                data = data[key]
            except KeyError:
                self.log.error(
                    f"Configuration file {self.filename} does not have "
                    f"data at node '{'/'.join(args)}'", exception=KeyError)

        return data

    @property
    def sra_meta(self) -> str:
        return self._return_with_validation("rna_seq", "meta")

    @property
    def sra_folder(self) -> str:
        return self._return_with_validation("rna_seq", "data_folder")

    @property
    def deseq2_meta(self) -> str:
        return self._return_with_validation("deseq2", "meta")

    @property
    def deseq2_folder(self) -> str:
        return self._return_with_validation("deseq2", "data_folder")

    @property
    def biomart_id_gene(self) -> str:
        return self._return_with_validation("biomart", "id_gene_rel")
