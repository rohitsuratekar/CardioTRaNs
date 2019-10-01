#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Network related models

from models.biology import Gene


class TRN:
    def __init__(self):
        self.components = []

    def add(self, *component):
        self.components.extend(component)

    def print_network(self):
        for c in self.components:
            c.print()

    def print_expression_pattern(self):
        for c in self.components:
            print("{} = {}".format(c, c.history))

    def update(self, time: int = 1):
        for t in range(time):
            for c2 in self.components:
                c = c2  # type:Gene
                c.output()
                c.update_expression()
