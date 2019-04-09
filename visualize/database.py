"""
CardioTRaNs 2019
Author: Rohit Suratekar

Visualization related to database
"""

from collections import defaultdict

import matplotlib.pylab as plt
from SecretColors import palette

from helpers.data_processing import get_expression_data


def check_expression_times():
    """
    General plot showing available expression data for given developmental time
    periods
    """
    p = palette.Palette()
    data, factors = get_expression_data(True)
    times = []

    normalize_time = defaultdict(list)

    for d in data.items():
        for ex in d[1]:
            if ex.time not in normalize_time[d[0]]:
                normalize_time[d[0]].append(ex.time)

    for d in normalize_time.items():
        for expression in d[1]:
            times.append(expression)

    differentiation = [x for x in times if x == 16]
    heart_tube = [x for x in times if x == 24]
    new_times = [x for x in times if x not in [16, 24]]
    plt.hist(differentiation, 50, color=p.orange(), range=(0, 50),
             label="Differentiation (16h)")
    plt.hist(heart_tube, 50, color=p.red(), range=(0, 50),
             label="Heart Tube (24h)")
    plt.hist(new_times, 50, color=p.cerulean(), range=(0, 50), label="Others")
    plt.ylabel("Unique Factors")
    plt.xlabel("Time (hpf)")
    plt.title("Expression Data Availability")
    plt.legend(loc=0)
    plt.show()


def check_factor_expression(gene: str):
    """
    Returns list of times at which given gene is expressed

    :param gene: symbol of gene (Case Sensitive)
    """
    data, factors = get_expression_data(True)
    factor = None
    for f in factors:
        if gene == factors[f].symbol:
            factor = factors[f]
            break

    if factor is None:
        raise Exception("{}not found in the database".format(gene))

    exp = []

    for key in data:
        if key == factor.uid:
            exp.extend(data[key])

    if len(exp) == 0:
        raise Exception("Expression pattern not found for {}".format(gene))

    return [x.time for x in exp]


def plot_expression_pattern(genes: list):
    """
    Shows expression of specific gene(s) available in the literature
    :param genes: list of genes to be analyze
    """
    p = palette.Palette()
    c = p.random(len(genes), grade=50)
    if type(c) is not list:
        c = [c]
    for i, gene in enumerate(genes):
        plt.hist(check_factor_expression(gene), 50, range=(0, 50),
                 label=gene, color=c[i], alpha=0.5)

    plt.ylabel("Unique Factors")
    plt.xlabel("Time (hpf)")
    plt.title("Expression Data Availability")
    plt.legend(loc=0)
    plt.show()


def check_all_factors(time: float):
    """
    Prints all the factors expressed at given time point
    :param time: Time in hours
    """
    data, factors = get_expression_data(True)
    exp_factors = []
    for d in data.values():
        for e in d:
            if e.time == time:
                exp_factors.append(e)

    exp_factors = set([x.factor.symbol for x in exp_factors])
    print("Factors at {} hpf: {}".format(time, ", ".join(exp_factors)))


def run():
    check_all_factors(30)
