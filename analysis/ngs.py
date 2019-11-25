#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# NGS Analysis

import pandas as pd
from SecretPlots import BooleanPlot
from SecretColors import Palette
from constants.boolean import *
from constants.ngs import *
from constants.system import *
from helpers.ngs_parser import get_string_tie_data, average_data


def get_data(hour) -> pd.DataFrame:
    d = get_string_tie_data(GENOTYPE_WT, hour, BIO_PROJECT_WINATA_LAB)
    d = average_data(d, OUTPUT_STRING_TIE)
    return d


def run():
    hours = [24, 48, 72]
    tpm_cutoffs = [0.84, 0.78, 0.82]
    genes = INTERESTED_GENES
    data = None
    for hour, tpm in zip(hours, tpm_cutoffs):
        d = get_data(hour)
        d = d[d[STRING_GENE_NAME].isin(genes)]
        d = d[[STRING_GENE_NAME, STRING_TPM]]
        d[str(hour)] = (d[STRING_TPM]
                        .apply(lambda x: 1 if x >= tpm else 0)
                        )
        d = d.set_index(STRING_GENE_NAME)
        d = d[[str(hour)]]
        if data is None:
            data = d
        else:
            data = pd.concat([data, d], axis=1)

    values = data.values
    p = Palette(show_warning=False)
    b = BooleanPlot(values, 1)
    (b.change_orientation("y")
     .add_y_ticklabels([x.capitalize() for x in data.index.values])
     .add_x_ticklabels(data.columns.values)
     .add_off_color(p.amber())
     .add_on_color(p.ultramarine())
     .add_x_padding(0, 0)
     .add_y_padding(0, 0)
     .add_x_label("Time (hpf)")
     .change_aspect_ratio(1 / 2)
     .add_x_midlines(color="w", alpha=0.5)
     .add_y_midlines(color="w", alpha=0.5)
     .add_legends(bbox_to_anchor=(1.05, 0.15))
     .show(tight=True))
