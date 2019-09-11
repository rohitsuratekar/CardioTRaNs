#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# RNA-SEQ related visualizations

from collections import defaultdict

import matplotlib.pylab as plt
from SecretColors import Palette

from helpers.parsers.ngs import *
from helpers.utils import log
from matplotlib.patches import Patch

p = Palette()


def _get_mapping_key(output_method):
    if output_method == OUTPUT_STAR:
        return STAR_UNIQUELY_MAPPED_READ_PERCENTAGE
    elif output_method == OUTPUT_SALMON_LOG:
        return SAL_LOG_MAPPING_PERCENT
    else:
        log("Mapping statistics only supports for {} and {}".format(
            OUTPUT_SALMON_LOG, OUTPUT_STAR), raise_error=True)


def plot_mapping_stats(output_method):
    meta = get_all_meta_data()
    names = []
    values = []
    for key in meta:
        names.append(meta[key][R_SEQ_FILE_REF_DETAILS])
        temp = []
        for files in get_rna_seq_files(output_method, key):
            v = parse(files, output_method)[_get_mapping_key(output_method)]
            v = str(v).replace("%", "")
            temp.append(float(v))
        values.append(temp)

    plt.boxplot(values, labels=names,
                patch_artist=True,
                boxprops=dict(facecolor=p.red(), color=p.black()),
                medianprops=dict(color=p.black())
                )
    plt.ylabel("Uniquely Mapped Reads (%)")
    plt.title(output_method)
    plt.gca().set_facecolor(p.gray(shade=15))
    plt.grid(axis="y", color=p.white())
    plt.tight_layout()
    plt.show()


def plot_test(data, width: float = 0.5,
              sample_gap: float = 1,
              internal_gap: float = 0.1,
              labels: list = None,
              sample_names: list = None,
              sample_colors: list = None,
              shade=True,
              shade_color=p.gray(shade=20),
              shade_alpha: float = 0.5):
    no_of_sample = 1
    only_data = False
    # Check if it has multiple samples
    try:
        if type(data[0][0]) == list:
            no_of_sample = max([len(x) for x in data])
        else:
            only_data = True
    except TypeError:
        # If it is only single sample with all data
        data = [data]
        only_data = True

    if sample_colors is None:
        if no_of_sample > len(p.get_color_list):
            sample_colors = p.random(no_of_colors=len(data))
        else:
            sample_colors = p.get_color_list[:no_of_sample]

    if labels is None:
        labels = [str(x) for x in range(len(data))]

    if sample_names is None:
        sample_names = ["G{}".format(x) for x in range(
            no_of_sample)]

    loc = []
    for l1 in range(0, len(data)):
        for k in range(no_of_sample):
            if len(loc) == 0:
                loc.append(0)
            else:
                if k == 0:
                    loc.append(loc[-1] + sample_gap + width)
                else:
                    loc.append(loc[-1] + internal_gap + width)

    for sample in range(no_of_sample):
        temp_loc = []
        temp_data = []
        for i in loc[sample::no_of_sample]:
            temp_loc.append(i)
        for d in data:
            temp_data.append(d[sample])

        if only_data:
            temp_data = data[0]

        plt.boxplot(temp_data,
                    positions=temp_loc,
                    widths=width,
                    patch_artist=True,
                    boxprops=dict(facecolor=sample_colors[sample],
                                  color=p.black()),
                    medianprops=dict(color=p.black()))

    bound_region = list(zip(loc[::no_of_sample],
                            loc[no_of_sample - 1::no_of_sample]))

    if shade:
        for s in bound_region:
            plt.axvspan(s[0] - width / 2, s[1] + width / 2,
                        alpha=shade_alpha,
                        color=shade_color)

    legend_elements = []

    for i in range(no_of_sample):
        legend_elements.append(Patch(facecolor=sample_colors[i],
                                     label=sample_names[i]))

    plt.legend(handles=legend_elements, loc=0)

    plt.xlim(0 - width, loc[-1] + width)
    plt.xticks([(x[0] + x[1]) / 2 for x in bound_region], labels)
    plt.show()


def test():
    bar_width = 0.5
    sample_loc = 1

    meta = get_all_meta_data()
    data = defaultdict(dict)
    methods = [OUTPUT_STAR, OUTPUT_SALMON_LOG]
    for key in meta:
        for m in methods:
            temp = []
            for f in get_rna_seq_files(m, key):
                v = parse(f, m)[_get_mapping_key(m)]
                v = str(v).replace("%", "")
                temp.append(float(v))

            data[key][m] = temp

    d = []
    for k in meta:
        d.append([data[k][OUTPUT_STAR], data[k][OUTPUT_SALMON_LOG]])

    plot_test(d[0])


def t_p(data):
    print(data)


def run():
    import numpy as np
    k = np.random.uniform(1, 100, 50)
    t_p(k)


