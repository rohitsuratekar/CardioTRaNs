#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Scores have multiplied by 1000 These functions are little bit heavy and
# you should use powerful system to run them. To test the code you can use
# "temp" parameter while getting data


from constants.string import *
from helpers.filemanager import string_info, string_links, string_actions
from pprint import pprint as pp


def get_info_mapping():
    info = string_info()
    info = info[[STRING_PROTEIN_EXTERNAL_ID, STRING_PREFERRED_NAME]]
    info = info.set_index(STRING_PROTEIN_EXTERNAL_ID)
    return info.to_dict("index")


def map_protein_names(d, info_mapping):
    data = d.copy()
    item1 = STRING_PROTEIN1
    item2 = STRING_PROTEIN2

    # Depending on input data, change the column names
    if item1 not in data.columns.values:
        item1 = STRING_ITEM1
        item2 = STRING_ITEM2

    for i in [item1, item2]:
        data[i] = data[i].map(lambda x: info_mapping[x][STRING_PREFERRED_NAME])

    return data


def find_interactors(links_data,
                     info_mapping,
                     gene: str,
                     confidence: float = 0.4):
    # p values of confidence level are multiplied with 1000 in the database

    links = links_data[links_data[STRING_EXPERIMENTS] +
                       links_data[
                           STRING_EXPERIMENTS_TRANSFER] >= confidence * 1000]

    links = links[links[STRING_COMBINED_SCORE] >= confidence * 1000]

    links = map_protein_names(links, info_mapping)
    links = links[links[STRING_PROTEIN1] == gene.strip().lower()]
    return links[STRING_PROTEIN2].values


def check_actions(actions, info_mapping, item1: str, item2: str):
    def _action_key(value):
        if value == "inhibition":
            return 0
        elif value == "activation":
            return 1
        else:
            return -1

    ac = map_protein_names(actions, info_mapping)
    ac = (ac[(ac[STRING_ITEM1] == item1) &
             (ac[STRING_ITEM2] == item2) &
             (ac[STRING_IS_DIRECTIONAL] == "t")])
    acts = []
    for index, row in ac.iterrows():
        is_acting = row[STRING_IS_ACTING]
        a = row[STRING_ACTION]
        if is_acting == "t":
            acts.append((item1, item2, _action_key(a)))
        else:
            acts.append((item2, item1, _action_key(a)))
    return list(set(acts))


def run():
    gene = "nkx2.5"
    data = string_links()
    info = get_info_mapping()
    names = find_interactors(data, info, gene)
    actions = string_actions()
    interactions = []
    for n in names:
        k = check_actions(actions, info, gene, n)
        if len(k) > 0:
            interactions.extend(k)
        else:
            interactions.append((gene, n, -1))

    pp(interactions)
