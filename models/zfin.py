"""
CardioTRaNs 2019
Author: Rohit Suratekar

All models related ZFIN database
Database columns can be found on https://zfin.org/downloads
"""


class ZFINAnatomy:
    """
    Class to hold anatomy data
    """

    def __init__(self, data: list):
        if len(data) != 4:
            raise Exception("Anatomy should have 4 items")
        self.id = data[0].strip()
        self.name = data[1].strip()
        self.start_stage_id = data[2].strip()
        self.end_stage_id = data[3].strip()

    @classmethod
    def empty(cls):
        return cls


class ZFINStages:
    """
    Class to hold stage data
    """

    def __init__(self, data: list):
        self.id = data[0].strip()
        self.obo_id = data[1].strip()
        self.name = data[2].strip()
        self.start = float(data[3])
        self.end = float(data[4])

    @classmethod
    def empty(cls):
        return cls
