"""
Project: CardioTrans
Author: Rohit Suratekar
Year: 2019

All models related NGS output files
"""


class STGeneExpression:
    """
    Simple class to hold the NGS information generated with StringTie tool
    """

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

    @property
    def is_expressed(self) -> bool:
        """
        Checks if gene is expressed or not
        :return: True if TPM is greater than 0
        """
        return self.tpm > 0
