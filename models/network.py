#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Network related models

from typing import List
from models.biology import Gene


class TRN:
    def __init__(self):
        self.components = []  # type: List[Gene]

    def add(self, *component):
        self.components.extend(component)

    def print_network(self):
        for c in self.components:
            c.print()

    def print_expression_pattern(self):
        for c in self.components:
            print("{}\t= {}".format(c, c.history))

    def expression_pattern(self):
        pattern = []
        names = []
        for c in self.components:
            pattern.append(c.history)
            names.append(c.name)

        return pattern, names

    def update(self, time: int = 1):
        for t in range(time):
            for g in self.components:
                g.check()
            for g in self.components:
                g.update()

    def regulatory_links(self) -> list:
        all_links = []
        for c in self.components:
            all_links.extend(c.outputs)

        return all_links
