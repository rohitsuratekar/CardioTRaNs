#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 12/08/19, 1:44 PM
#
#  Copyright (c) 2019.
#
#  All network related models will go here


class TRN:
    def __init__(self):
        self.components = []

    def add(self, *component):
        self.components.extend(component)

    def print(self):
        for c in self.components:
            c.print()

    def print_expression_pattern(self):
        for c in self.components:
            print("{} = {}".format(c, c.is_expressed))
