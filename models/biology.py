#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 12/08/19, 4:03 PM
#
#  Copyright (c) 2019.

from models.boolean import *


class Gene:
    def __init__(self, name: str, is_expressed: bool = True):
        self.name = name.lower().strip()
        self.__expression = is_expressed
        self.__history = [is_expressed]
        self.__input_function = None

    @property
    def is_expressed(self):
        return self.__expression

    @is_expressed.setter
    def is_expressed(self, value: bool):
        self.__expression = value

    @property
    def history(self):
        return self.__history

    def __bool__(self):
        return self.is_expressed

    __nonzero__ = __bool__

    def __int__(self):
        return self.is_expressed

    def __add__(self, new):
        return self.__int__() + self.__int__()

    def __str__(self):
        return self.name

    def input_function(self, expression):
        self.__input_function = expression

    def enhances(self, *args):
        pass

    def represses(self, *args):
        pass

    def input(self, input_function: Input):
        self.__input_function = input_function
