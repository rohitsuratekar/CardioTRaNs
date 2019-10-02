#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Some modified versions of the Boolean objects

from constants.boolean import SYMBOL_AND, SYMBOL_NOT, SYMBOL_OR, SYMBOL_COPY


class Input:
    def __init__(self, *args):
        self._values = args
        self._simplified = None

        if len(self._values) < 1:
            raise Exception("At lest one component is needed for {} "
                            "gate".format(self.name))

    @property
    def name(self):
        raise NotImplementedError

    @property
    def values(self):
        return self._values

    @property
    def types(self):
        return [type(x) == bool for x in self.values]

    def solve(self):
        raise NotImplementedError

    def _simplify(self):
        x = []
        for v in self.values:
            if type(v) == bool:
                x.append(v)
            elif issubclass(type(v), Input):
                x.append(v.solve())
            elif type(v) == int:
                if v not in [0, 1]:
                    raise Exception("If you are using integer for Boolean "
                                    "function, then only 1 and 0 are allowed")
                else:
                    x.append(v.__bool__())
            else:
                try:
                    x.append(v.__bool__())
                except AttributeError:
                    raise Exception("{} is not a valid input for logic "
                                    "classes".format(type(v)))

        self._simplified = x


class OR(Input):

    @property
    def name(self):
        return "OR"

    def solve(self):
        self._simplify()
        return any(self._simplified)

    def __str__(self):
        return SYMBOL_OR


class AND(Input):

    @property
    def name(self):
        return "AND"

    def solve(self):
        self._simplify()
        return all(self._simplified)

    def __str__(self):
        return SYMBOL_AND


class NOT(Input):

    @property
    def name(self):
        return "NOT"

    def solve(self):
        self._simplify()
        if len(self._simplified) != 1:
            raise Exception("NOT gate will require ONLY one input")
        return not self._simplified[0]

    def __str__(self):
        return SYMBOL_NOT


class COPY(Input):

    @property
    def name(self):
        return "COPY"

    def solve(self):
        self._simplify()
        if len(self._simplified) != 1:
            raise Exception("COPY gate will require ONLY one input")
        return self._simplified[0]

    def __str__(self):
        return SYMBOL_COPY


def run():
    pass
