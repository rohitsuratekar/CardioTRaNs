#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 12/08/19, 1:44 PM
#
#  Copyright (c) 2019.
#
#  All network related models will go here

from models.biology import Gene


class TRN:
    def __init__(self):
        self.components = []

    def add(self, component: Gene):
        self.components.append(component)

    def print(self):
        for c in self.components:
            c.print()
