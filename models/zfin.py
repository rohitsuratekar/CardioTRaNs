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

    def __init__(self, data: list = None, relationship: dict = None):
        """
        :param data: Columns from FILE_ANATOMY_ITEMS file
        """
        self.id = None
        self.name = None
        self.start_stage_id = None
        self.end_stage_id = None
        self.parent_id = None

        if data is not None and relationship is not None:
            if len(data) != 4:
                raise Exception("Anatomy should have 4 items")
            self.id = data[0].strip()
            self.name = data[1].strip()
            self.start_stage_id = data[2].strip()
            self.end_stage_id = data[3].strip()
            try:
                self.parent_id = relationship[self.id][0]
            except KeyError:
                self.parent_id = self.id


class ZFINStages:
    """
    Class to hold stage data
    """

    def __init__(self, data: list = None):
        """
        :param data: Columns from FILE_STAGE_ONTOLOGY file
        """

        self.id = None
        self.obo_id = None
        self.name = None
        self.start = None
        self.end = None

        if data is not None:
            if len(data) != 5:
                raise Exception("Stage list should have 5 columns")

            self.id = data[0].strip()
            self.obo_id = data[1].strip()
            self.name = data[2].strip()
            self.start = float(data[3])
            self.end = float(data[4])


class ZFINExpression:
    """
    Class to hold wild type expression data
    """

    def __init__(self, data: list = None, stages: dict = None):
        """
        :param data: Columns from FILE_WILD_TYPE_EXPRESSION file
        :param stages: Dictionary of stages from FILE_STAGE_ONTOLOGY file
        with their names as a key
        """

        self.id = None
        self.gene_symbol = None
        self.fish_name = None
        self.super_structure_id = None
        self.sub_structure_id = None
        self.start_stage_id = None
        self.end_stage_id = None
        self.assay = None
        self.assay_mmo_id = None
        self.publication_id = None
        self.probe_id = None
        self.antibody_id = None
        self.fish_id = None

        if data is not None and stages is not None:

            if len(data) != 15:
                raise Exception("Expression list should have 15 columns")

            self.id = data[0]
            self.gene_symbol = data[1]
            self.fish_name = data[2]
            self.super_structure_id = data[3]
            self.sub_structure_id = data[5]
            self.start_stage_id = stages[data[7]]
            self.end_stage_id = stages[data[8]]
            self.assay = data[9]
            self.assay_mmo_id = data[10]
            self.publication_id = data[11]
            self.probe_id = data[12]
            self.antibody_id = data[13]
            self.fish_id = data[14]
