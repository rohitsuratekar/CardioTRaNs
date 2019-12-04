#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Clustering algorithms

import numpy as np
import pandas as pd

from constants.other import *
from constants.string import *
from helpers.filemanager import string_info, string_links


def get_temp_data():
    with open("analysis/test.txt") as f:
        return pd.read_csv(f, delimiter=";").values


def get_data():
    k = string_info(organism=ORG_ZEBRAFISH)
    k = pd.Series(k[STRING_PREFERRED_NAME].values, index=k[
        STRING_PROTEIN_EXTERNAL_ID].values)

    links = string_links("original", organism=ORG_ZEBRAFISH, use_dask=True)

    links[STRING_PROTEIN1] = links[STRING_PROTEIN1].map(lambda x: k[x],
                                                        meta=(
                                                            STRING_PROTEIN1,
                                                            str))
    links[STRING_PROTEIN2] = links[STRING_PROTEIN2].map(lambda x: k[x],
                                                        meta=(
                                                            STRING_PROTEIN2,
                                                            str))

    return links[[STRING_PROTEIN1, STRING_PROTEIN2]].values


def generate_cluster(data):
    clusters = []
    for d in data:
        values = [len(set(x).intersection(d)) for x in clusters]
        if sum(values) == 0:
            clusters.append(list(d))
        else:
            tmp = list(d)
            tmp2 = []
            for i in range(len(values)):
                if values[i] != 0:
                    tmp.extend(clusters[i])
                else:
                    tmp2.append(clusters[i])

            tmp2.append(list(set(tmp)))
            clusters = tmp2

    return clusters


class ClusterNode:
    def __init__(self, name):
        self.name = name
        self.value = np.inf
        self.closest = None
        self.alternatives = []

    def add_distance(self, value, node):
        if value < self.value:
            self.value = value
            self.closest = node
            return True
        elif value == self.value and value != np.inf:
            new_alt = []
            for a in self.alternatives:
                if a[0] < value:
                    new_alt.append(a)
            if self.closest != node:
                new_alt.append((value, node))
            self.alternatives = new_alt

        return False

    def __repr__(self):
        return f"Node( {self.value} )"


def calculate_distance(data, start_node):
    all_nodes = list(set(data.flatten()))
    all_nodes = {x: ClusterNode(x) for x in all_nodes}
    all_nodes[start_node].value = 0
    counter = 1
    while counter != len(all_nodes):
        temp_count = 0
        for d in data:
            d1 = all_nodes[d[0]]  # type:ClusterNode
            d2 = all_nodes[d[1]]  # type:ClusterNode
            if d2.add_distance(d1.value + 1, d1.name):
                temp_count += 1
            if d1.add_distance(d2.value + 1, d2.name):
                temp_count += 1

        if temp_count == 0:
            break
        else:
            counter += temp_count

    return all_nodes


def find_connection(data, gene1, gene2):
    if gene2 not in data.flatten():
        raise Exception(f"{gene2} is not present in the database")
    d = calculate_distance(data, gene1)
    path = []
    next_gene = gene2
    while True:
        if next_gene is None:
            break
        ng = d[next_gene]  # type: ClusterNode
        if ng is not None:
            if len(ng.alternatives) > 0:
                alt = [x[1] for x in ng.alternatives]
                alt.append(ng.name)
                if ng.name == gene2:
                    path.append(ng.name)
                else:
                    path.append(alt)
            else:
                path.append(ng.name)
            next_gene = ng.closest
        else:
            path.append(ng.name)
    path = list(reversed(path))
    print(path)


def run():
    d = get_data()
    print(d)
    # find_connection(d, "nkx2.5", "gata5")
