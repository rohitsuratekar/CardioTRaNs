#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Models related to Biology

from models.boolean import *


class Gene:
    def __init__(self, name: str, is_expressed: bool = True):
        self.name = name.lower().strip()
        self.__expression = is_expressed
        self.__history = [is_expressed]
        self.__input_function = None
        self.__temp_expression = None

    @property
    def is_expressed(self):
        return self.__expression

    @is_expressed.setter
    def is_expressed(self, value: bool):
        self.__expression = value

    @property
    def history(self):
        return self.__history

    def update_history(self):
        self.__history.append(self.is_expressed)

    def __bool__(self):
        return self.is_expressed

    __nonzero__ = __bool__

    def __int__(self):
        return self.is_expressed

    def __add__(self, new):
        return self.__int__() + self.__int__()

    def __str__(self):
        return self.name

    def input(self, expression: Input):
        self.__input_function = expression

    def output(self):
        if self.__input_function is None:
            self.__temp_expression = self.is_expressed
            return self.is_expressed
        print(self.__input_function)
        self.__temp_expression = self.__input_function.solve()
        return self.__temp_expression

    def update_expression(self):
        if self.__temp_expression is None:
            self.__temp_expression = self.is_expressed
        self.is_expressed = self.__temp_expression
        self.update_history()

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

        if self.__input_function is None:
            print("There is no expression function registered for {}".format(
                self.name))

        else:
            ip = self.__input_function.inputs
            if len(ip) == 1:
                if type(self.__input_function) == COPY:
                    print("{}*\t= {}".format(self.name, ip[0]))
                else:
                    print("{}*\t= {} ( {} ) ".format(self.name,
                                                     self.__input_function,
                                                     ip[0]))
            else:
                t = _expand(ip, str(self.__input_function))
                while any([issubclass(type(x), Input) for x in t]):
                    for m in t:
                        if issubclass(type(m), Input):
                            t.insert(t.index(m), "(")
                            for m2 in _expand(m.inputs, str(m)):
                                t.insert(t.index(m), m2)

                            t.insert(t.index(m), ")")
                            t.remove(m)

                print("{}*\t= {}".format(self.name,
                                         " ".join([str(x) for x in t])))
