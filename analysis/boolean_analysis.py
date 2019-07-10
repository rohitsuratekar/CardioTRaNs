"""
Project: CardioTrans
Author: Rohit Suratekar
Year: 2019

This file should contain all functions related to Boolean analysis
"""

from helpers.parsers import get_boolean_genes


def get_initial_condition():
    data = get_boolean_genes()
    print(data)


def run():
    get_initial_condition()
