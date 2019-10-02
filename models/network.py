#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Network related models


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
            print("{}\t= {}".format(c, c.history))

    def update(self, time: int = 1):
        for t in range(time):
            for g in self.components:
                g.check()
            for g in self.components:
                g.update()
