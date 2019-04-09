"""
CardioTRaNs 2019
Author: Rohit Suratekar

All Biological Object Models
"""

from constants import REGULATION_TYPES


class Factor:
    """
    Class to hold information regarding factors
    """

    def __init__(self, data: list):
        """
        :param data: List of column values of INPUT_FILE_FACTORS
        """
        try:
            self.uid = int(data[0])  # UID of the factor
        except ValueError:
            raise Exception("UID should be integer in line {}".format(data))
        self.symbol = data[1]  # Symbol
        self.organism = data[2]  # Organism
        self.zfin = data[3]  # ZFIN ID
        self.type = data[4]  # Type
        self.homolog = int(data[5])  # Zebrafish Homolog


class RegulationLink:
    """
    Class to hold target information
    """

    def __init__(self, data: list, factor_map: dict):
        """
        :param data: List of column values of INPUT_FILE_TARGETS
        :param factor_map: Dictionary of factors
        """
        try:
            factor_uid = int(data[0])  # UID of the factor
            target_uid = int(data[2])  # UID of the target
        except ValueError:
            raise Exception("UID should be integer in line {}".format(data))

        try:
            self.factor = factor_map[factor_uid]
        except KeyError:
            raise Exception("Factor '{}' from the line {} is not present in "
                            "the 'factor_map' provided. Please check if it is"
                            "present in INPUT_FILE_FACTORS file".
                            format(data[1], data))
        try:
            self.target = factor_map[target_uid]
        except KeyError:
            raise Exception("Target '{}' from the line {} is not present in "
                            "the 'factor_map' provided. Please check if it is"
                            "present in INPUT_FILE_FACTORS file".
                            format(data[3], data))

        self.regulation = str(data[4]).lower()  # Type of regulation
        if self.regulation not in REGULATION_TYPES:
            raise Exception(
                "Regulation type '{}' is not allowed in the line {}\nAllowed "
                "types are {}".format(data[4], data, REGULATION_TYPES))

        self.place = data[5]  # Place of regulation
        self.reported_stage = data[6]  # Reported Developmental Stage
        self.stage = data[7]  # Corrected developmental stage
        try:
            self.time = float(data[8])
        except ValueError:
            raise Exception("Time should be a number in the line \n{}".format(
                data))
        self.experiment = data[9]  # Type of experiment
        self.reference = data[10]  # Reference
        self.notes = data[11]  # Notes


class Expression:
    def __init__(self, data: list, factor_map):
        try:
            uid = int(data[0])  # UID of the factor
        except ValueError:
            raise Exception("UID should be integer in line {}".format(data))

        try:
            self.factor = factor_map[uid]
        except KeyError:
            raise Exception("Factor '{}' from the line {} is not present in "
                            "the 'factor_map' provided. Please check if it is"
                            "present in INPUT_FILE_FACTORS file".
                            format(data[1], data))
        self.place = data[2]  # Place of expression
        self.reported_stage = data[3]  # Reported/Calculated developmental stage
        self.stage = data[4]  # Corrected developmental stage
        try:
            self.time = float(data[5])
        except ValueError:
            raise Exception("Time should be a number in the line \n{}".format(
                data))
        self.experiment = data[6]  # Type of experiment
        self.reference = data[7]  # Reference
        self.notes = data[8]  # Notes
