#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
#   All statistics related functions

import pandas as pd
from SecretColors import Palette
from SecretPlots import BarGroupedPlot


def star_mapping_stat(path):
    data = {}
    with open(path) as f:
        for line in f:
            if len(line.strip()) != 0:
                values = line.strip().split("|")
                values = [x.strip() for x in values]
                if len(values) == 2:
                    data[values[0]] = values[1]

    input_reads = data["Number of input reads"]
    unique_per = data["Uniquely mapped reads %"].replace("%", "")
    average_length = data["Average input read length"]
    return [int(input_reads), float(unique_per), int(average_length)]


def get_all_star_stat(bioproject: str):
    data = pd.read_csv("samples.csv")
    data = data[data["BioProject"] == bioproject]
    base = "/mnt/windows/Enigma/Zebrafish/data/rnaseq/star"
    all_times = list(sorted(set(data["time"].values)))
    values = []
    for time in all_times:
        temp = data[data["time"] == time]
        temp = temp[temp["genetic_background"] == "wild type"]
        time_temp = []
        for srr in temp["Run"].values:
            url = f"{base}/{srr}/{srr}_Log.final.out"
            reads, mapped, average = star_mapping_stat(url)
            time_temp.append(reads)
        values.append(sorted(time_temp))

    p = Palette()
    b = BarGroupedPlot(values)
    b.show_legend = False
    b.add_x_ticklabels(all_times)
    b.add_x_label("Hours Post Fertilization")
    b.add_y_label("Number of Input Reads")
    b.add_group_gap(3)
    b.colors = [p.red(shade=30), p.red(), p.red(shade=70)]
    b.draw()
    b.ax.set_title("Yost Lab")
    b.ax.grid(axis="y", ls="--", zorder=0)
    b.ax.set_ylim(0, 48889779)
    b.show()


def run():
    project = "PRJNA407368"
    # project = "PRJNA492280"
    get_all_star_stat(project)
