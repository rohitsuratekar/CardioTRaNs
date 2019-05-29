"""
CardioTRaNs 2019
Author: Rohit Suratekar

Models related to NGS data-sets
"""


class STGeneExpression:
    def __init__(self, data):
        self.data = data
        self.gene_id = data[0]
        self.gene_name = data[1]
        self.reference = data[2]
        self.strand = data[3]
        self.start = int(data[4])
        self.end = int(data[5])
        self.coverage = float(data[6])
        self.fpkm = float(data[7])
        self.tpm = float(data[8])
