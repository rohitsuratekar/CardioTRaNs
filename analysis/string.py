#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Scores have multiplied by 1000 These functions are little bit heavy and
# you should use powerful system to run them. To test the code you can use
# "temp" parameter while getting data


import pandas as pd

from constants.other import *
from constants.string import *
from constants.zfin import *
from helpers.filemanager import string_info, string_links, zfin_ortho

from pprint import pprint as pp


class Score:
    """
    Simple class to hold the column information
    """

    def __init__(self, name):
        self.name = name
        self.use = False


class NetworkFinder:
    """
    Class used to calculate the network from STRING database
    """

    def __init__(self):
        self._info = None
        self._links = None
        self._file_status = "original"
        self.organism = ORG_ZEBRAFISH
        self.neighborhood = Score(STRING_NEIGHBORHOOD)
        self.neighborhood_trans = Score(STRING_NEIGHBORHOOD_TRANSFER)
        self.fusion = Score(STRING_FUSION)
        self.cooccurence = Score(STRING_COOCCURENCE)
        self.homology = Score(STRING_HOMOLOGY)
        self.coexpression = Score(STRING_COEXPRESSION)
        self.coexpression_trans = Score(STRING_COEXPRESSION_TRANSFER)
        self.experiments = Score(STRING_EXPERIMENTS)
        self.experiments_trans = Score(STRING_EXPERIMENTS_TRANSFER)
        self.text_mining = Score(STRING_TEXT_MINING)
        self.text_mining_trans = Score(STRING_TEXT_MINING_TRANSFER)

        self.experiments.use = True
        self.experiments_trans.use = True

        self.threshold = 0.4

    @property
    def columns(self):
        """
        :return: List of all columns
        """
        return [self.neighborhood, self.neighborhood_trans, self.fusion,
                self.cooccurence, self.homology, self.coexpression,
                self.coexpression_trans, self.experiments,
                self.experiments_trans, self.text_mining,
                self.text_mining_trans]

    @property
    def ids(self):
        """
        :return: Dictionary or Mapping of Protein IDs (as index) and
        preferred name (as value)
        """
        if self._info is None:
            k = string_info(organism=self.organism)
            k = pd.Series(k[STRING_PREFERRED_NAME].values, index=k[
                STRING_PROTEIN_EXTERNAL_ID].values)
            self._info = k.to_dict()
        return self._info

    @property
    def interactions(self) -> pd.DataFrame:
        """
        Protein-Protein interactions and their individual scores.
        """
        if self._links is None:
            # Get data
            self._links = string_links(self._file_status,
                                       organism=self.organism)
            # Convert IDs to regular names
            self._links[STRING_PROTEIN1] = self._links[STRING_PROTEIN1].map(
                lambda x: self.ids[x])
            self._links[STRING_PROTEIN2] = self._links[STRING_PROTEIN2].map(
                lambda x: self.ids[x])

            # Calculate temporary score for given traits
            temp = "temp"
            self._links[temp] = 0
            for c in self.columns:
                if c.use:
                    self._links[temp] = self._links[temp] + self._links[c.name]

            # Take only entries which have more score than the user set
            # threshold (database values are multiplied by 1000)
            self._links = (self._links[
                self._links[temp] >= self.threshold * 1000].reset_index(
                drop=True))

            # Remove the temporary column
            del self._links[temp]

        return self._links

    def find_interactors(self, gene: str):
        """
        Finds interacting proteins with given settings
        :param gene: Gene name
        :return: List of interacting partners
        """
        data = self.interactions[
            self.interactions[STRING_PROTEIN1].str.lower() == gene.strip(
            ).lower()]
        return list(data[STRING_PROTEIN2].values)

    def generate_network(self, genes, *, level: int = 1):
        # Sanity check
        if level < 1 or type(level) != int:
            raise Exception(
                "Invalid Level. Please use integers greater than 0")
        # if only single gene is given, convert it into a list
        if type(genes) == str:
            genes = [genes]
        net_stat = []
        all_genes = []
        secondary_genes = []
        # Start network expansion
        for i in range(level):
            temp = []
            if len(net_stat) == 0:
                secondary_genes = genes
            for g in secondary_genes:
                if g not in all_genes:
                    links = self.find_interactors(g)
                    net_stat.append((g, links))
                    all_genes.append(g)
                    temp.extend(links)
            secondary_genes = temp
            all_genes = list(set(all_genes))

        return net_stat


def check_homologue(data, gene_name, organism):
    if organism == ORG_ZEBRAFISH:
        return gene_name
    d = data[data[ZFIN_ORTHO_GENE_SYMBOL] == gene_name.strip().lower()]
    col = ZFIN_ORTHO_HUMAN_SYMBOL
    if organism == ORG_MOUSE:
        col = ZFIN_ORTHO_MOUSE_SYMBOL
    d = list(set(d[col].values))
    if len(d) > 1:
        raise Exception(f"Something went wrong with {gene_name} in {organism}")
    elif len(d) == 0:
        raise Exception(f"No homologue information found for {gene_name} in "
                        f"{organism}")
    return d[0]


def run():
    n = NetworkFinder()
    n.organism = ORG_HUMAN
    d = zfin_ortho(n.organism)
    gene = check_homologue(d, "nfasca", n.organism)
    k = n.generate_network(gene)
    pp(k)
