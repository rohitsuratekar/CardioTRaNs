#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Models related to Biology

from models.boolean import *


class Gene:
    def __init__(self, name, is_expressed: bool = True):
        self.name = name
        self._is_expressed = is_expressed.__bool__()
        self._input_function = None
        self._future_state = None
        self._history = [self.is_expressed]
        self._outputs = []

    def __bool__(self):
        return self.is_expressed

    def __str__(self):
        return self.name

    @property
    def is_expressed(self):
        return self._is_expressed

    @property
    def history(self):
        return self._history

    def input(self, function: Input):
        self._input_function = function

    @property
    def outputs(self):
        return self._outputs

    def add_output(self, gene: str, link_type: int):
        self._outputs.append([self.name, gene, link_type])

    def check(self):
        if self._input_function is None:
            self._future_state = self.is_expressed
        else:
            self._future_state = self._input_function.solve()

    def update(self):
        if self._future_state is None:
            raise Exception("To update the gene expression state, "
                            "run Gene.check() first")
        self._is_expressed = self._future_state
        self._history.append(self._future_state)
        self._future_state = None

    def print(self):

        def _expand(k, class_type: str):
            if class_type == SYMBOL_NOT:
                return [class_type, k[0]]

            j = []
            for i2 in k:
                j.append(i2)
                j.append(class_type)

            j.pop()
            return j

        if self._input_function is None:
            print("{0}*\t= {0}".format(self.name))

        else:
            ip = self._input_function.values
            if len(ip) == 1:
                if type(self._input_function) == COPY:
                    print("{}*\t= {}".format(self.name, ip[0]))
                else:
                    print("{}*\t= {} ( {} ) ".format(self.name,
                                                     self._input_function,
                                                     ip[0]))
            else:
                t = _expand(ip, str(self._input_function))
                while any([issubclass(type(x), Input) for x in t]):
                    for m in t:
                        if issubclass(type(m), Input):
                            t.insert(t.index(m), "(")
                            for m2 in _expand(m.values, str(m)):
                                t.insert(t.index(m), m2)

                            t.insert(t.index(m), ")")
                            t.remove(m)

                print("{}*\t= {}".format(self.name,
                                         " ".join([str(x) for x in t])))
