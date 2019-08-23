#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 12/08/19, 3:05 PM
#
#  Copyright (c) 2019.
#
#  All Boolean object models will go here

from constants.boolean import SYMBOL_AND, SYMBOL_NOT, SYMBOL_OR


class Input:
    """
    Base class for all the logic gates.
    Do not use this class directly, use individual classes.
    """

    def __init__(self, *args):
        self._input = args
        self._final_inputs = []

    def solve(self):
        raise NotImplementedError("Subclass must implement solve")

    def _resolve_input(self):
        for x in self._input:
            if type(x) == bool:
                self._final_inputs.append(x)
            elif issubclass(type(x), Input):
                self._final_inputs.append(x.solve())
            elif type(x) == int:
                if x not in [0, 1]:
                    raise Exception("If you are using integer for Boolean "
                                    "function, then only 1 and 0 are allowed")
                else:
                    self._final_inputs.append(x.__bool__())
            else:
                try:
                    self._final_inputs.append(x.__bool__())
                except AttributeError:
                    raise Exception("{} is not a valid input for logic "
                                    "classes".format(type(x)))

    @property
    def inputs(self):
        return self._input


class OR(Input):
    """
    Simple OR Gate
    """

    def __init__(self, *args):
        super().__init__(*args)

        if len(args) < 2:
            raise Exception("OR gate should have at least 2 inputs")

    def solve(self):
        if len(self._final_inputs) == 0:
            self._resolve_input()

        return any(self._final_inputs)

    def __str__(self):
        return SYMBOL_OR


class AND(Input):
    """
    Simple AND Gate.
    """

    def __init__(self, *args):
        super().__init__(*args)

        if len(args) < 2:
            raise Exception("AND gate should have at least 2 inputs")

    def solve(self):
        if len(self._final_inputs) == 0:
            self._resolve_input()

        return all(self._final_inputs)

    def __str__(self):
        return SYMBOL_AND


class NOT(Input):
    """
    Simple NOT Gate
    """

    def __init__(self, *args):
        super().__init__(*args)

        if len(args) > 1:
            raise Exception("NOT gate should have only 1 input")

    def solve(self):
        if len(self._final_inputs) == 0:
            self._resolve_input()

        return not self._final_inputs[0]

    def __str__(self):
        return SYMBOL_NOT
