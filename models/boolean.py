#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 12/08/19, 3:05 PM
#
#  Copyright (c) 2019.
#
#  All Boolean object models will go here


class Input:
    def __init__(self, *args):
        self.args = args


class OR(Input):
    def __init__(self, *args):
        super().__init__(*args)
        self.inputs = args


class AND(Input):
    def __init__(self, *args):
        super().__init__(*args)
        self.inputs = args
    #     for i in args:
    #         self.inputs.append(i.__bool__())
    #
    # def __bool__(self):
    #     return not (False in self.inputs)
    #
    # __nonzero__ = __bool__


class NOT(Input):
    def __init__(self, *args):
        super().__init__(*args)
        self.inputs = args

        if len(self.inputs) > 1:
            raise Exception("NOT gate takes only 1 input")
