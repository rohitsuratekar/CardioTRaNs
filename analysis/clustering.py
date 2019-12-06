#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Clustering algorithms

from collections import Counter

import numpy as np
import pandas as pd
from SecretColors import Palette
from SecretPlots import NetworkPlot
from constants.other import *
from analysis.string import NetworkFinder, check_homologue
from constants.string import *
from helpers.filemanager import zfin_ortho


def get_temp_data():
    with open("analysis/test.txt") as f:
        return pd.read_csv(f, delimiter=";").values


def get_data(organism=ORG_ZEBRAFISH):
    nf = NetworkFinder()
    nf.organism = organism
    return nf.interactions[[STRING_PROTEIN1, STRING_PROTEIN2]].values


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


def find_shortest_connection(data, gene1, gene2):
    if gene2 not in data.flatten():
        raise Exception(f"{gene2} is not present in the database")
    d = calculate_distance(data, gene1)
    all_paths = []

    def _find_path(data_list, temp_ng):
        temp_path = data_list
        while True:
            extra_paths = []
            if temp_ng is None:
                break
            ng = d[temp_ng]  # type: ClusterNode
            if ng is not None:
                if len(ng.alternatives) > 0:
                    alt = [x[1] for x in ng.alternatives]
                    extra_paths.extend(alt)
                    temp_path.append(ng.name)
                else:
                    temp_path.append(ng.name)
                temp_ng = ng.closest
            else:
                temp_path.append(ng.name)

            for e in extra_paths:
                ex_path = [x for x in temp_path]
                _find_path(ex_path, e)

        temp_path = list(reversed(temp_path))
        all_paths.append(temp_path)

    _find_path([], gene2)
    return all_paths


def paths_to_connections(data, paths):
    connections = []
    for path in paths:
        for i in range(1, len(path)):
            start = path[i - 1]
            end = path[i]
            row = np.array([start, end])
            ind = np.where(np.all(data == row, axis=1))
            if len(data[ind]) == 0:
                connections.append([end, start, 1])
            else:
                connections.append([start, end, 1])

    return connections


def run():
    organism = ORG_MOUSE
    p = Palette()
    data = zfin_ortho(organism)
    gene_counter = Counter()
    zf_test_genes = ["mapk1", "mapk3", "ptch1", "actc1a", "actb1", "smarca5",
                     "myh6", "spns2", "sptb", "sptan1", "akt1", "akt2",
                     "akt3b"]
    test_genes = []
    for x in zf_test_genes:
        try:
            test_genes.append(check_homologue(data, x, organism))
        except Exception:
            print(f"No homologue found for {x} in {organism}")

    for end_gene in test_genes:
        d = get_data(organism)
        # start_gene = "chl1a"
        start_gene = "Chl1"
        try:
            paths = find_shortest_connection(d, start_gene, end_gene)
        except Exception:
            print(f"No connection found in {end_gene}")
            continue
        data = paths_to_connections(d, paths)
        for d in data:
            gene_counter.update({d[0]})
            gene_counter.update({d[1]})
        n = NetworkPlot(data)
        n.max_columns = 5
        n.node_width = 2
        n.node_gap = 2
        n.aspect_ratio = 1
        n.add_text_options(fontsize=10)
        n.node_placement = paths[0]
        n.line_decoration = False
        n.colors = p.gray(shade=30)
        # n.colors_mapping = {"ptprz1a": p.blue(), "chl1a": p.green(),
        #                     "ptprga": p.blue(), "ca16b": p.blue(),
        #                     "cntn1a": p.blue(), end_gene: p.red()}
        n.colors_mapping = {"Ptprz1": p.blue(), "Chl1": p.green(),
                            "Ptprg": p.blue(), "Ca16b": p.blue(),
                            "Cntn1": p.blue(), end_gene: p.red()}

        n.save(f"network_{end_gene}.png", tight=True, dpi=300, format="png")

    print(gene_counter)
