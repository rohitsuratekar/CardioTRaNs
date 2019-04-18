"""
CardioTRaNs 2019
Author: Rohit Suratekar

Local database operations
"""

import time

import h5py

from constants import *
from models.zfin import *


class LocalDatabase:
    """
    All database related operations and functions
    """

    # VERY IMPORTANT, Always open with "a" for writing else old data will
    #  be erased
    __mode_writing = "a"
    __mode_reading = "r"

    def __init__(self):
        try:
            with h5py.File(DATABASE, self.__mode_writing) as db:
                main = db.create_group(DB_DATA)
                main.create_group(DB_ANATOMY)
                main.create_group(DB_ANATOMY_RELATIONSHIPS)
                main.create_group(DB_STAGE_ONTOLOGY)
        except ValueError:
            pass

        self.dt = h5py.special_dtype(vlen=str)  # Need to store variable
        # length strings

    @property
    def get_group_names(self):
        """
        :return: Names of tables available in the main database
        """
        with h5py.File(DATABASE, self.__mode_reading) as db:
            return list(db[DB_DATA].keys())

    def add_anatomy_items(self, items: list) -> None:
        """
        Adds anatomy item to the database
        """

        with h5py.File(DATABASE, self.__mode_writing) as db:
            ai = db[DB_DATA][DB_ANATOMY]
            for item in items:
                anatomy = ai.create_dataset(
                    item.id, (1, 4), maxshape=(None, None), dtype=self.dt)
                anatomy[0, DB_ANATOMY_ID] = item.id
                anatomy[0, DB_ANATOMY_NAME] = item.name
                anatomy[0, DB_ANATOMY_START_STAGE] = item.start_stage_id
                anatomy[0, DB_ANATOMY_END_STAGE] = item.end_stage_id

                # Add meta-data
                anatomy.attrs[DB_ATTRS_AUTHOR] = DATABASE_AUTHOR
                anatomy.attrs[DB_ATTRS_VERSION] = DATABASE_VERSION
                anatomy.attrs[DB_ATTRS_CREATION] = int(time.time())

    def add_stages(self, items: list) -> None:
        """
        Adds Stage information into the database
        """
        with h5py.File(DATABASE, self.__mode_writing) as db:
            st = db[DB_DATA][DB_STAGE_ONTOLOGY]
            for item in items:
                stage = st.create_dataset(
                    item.id, (1, 5), maxshape=(None, None), dtype=self.dt)

                stage[0, DB_STAGE_ID] = item.id
                stage[0, DB_STAGE_NAME] = item.name
                stage[0, DB_STAGE_OBO_ID] = item.obo_id
                stage[0, DB_STAGE_START] = item.start
                stage[0, DB_STAGE_END] = item.end

                # Add meta-data
                stage.attrs[DB_ATTRS_AUTHOR] = DATABASE_AUTHOR
                stage.attrs[DB_ATTRS_VERSION] = DATABASE_VERSION
                stage.attrs[DB_ATTRS_CREATION] = int(time.time())

    def get_anatomy_item(self, anatomy_id: str) -> ZFINAnatomy:
        """
        :param anatomy_id: ZFIN anatomy ID
        :return: ZFINAnatomy Object
        :raises: KeyError if anatomy ID not found
        """
        z = ZFINAnatomy()
        with h5py.File(DATABASE, self.__mode_reading) as db:
            z.id = anatomy_id
            ref = db[DB_DATA][DB_ANATOMY][z.id][0]
            z.name = ref[DB_ANATOMY_NAME]
            z.start_stage_id = ref[DB_ANATOMY_START_STAGE]
            z.end_stage_id = ref[DB_ANATOMY_END_STAGE]
        return z

    def get_stage(self, stage_id: str) -> ZFINStages:
        """
        :param stage_id: ZFIN stage ID
        :return: ZFINStages object
        :raises: KeyError if stage ID not found
        """
        with h5py.File(DATABASE, self.__mode_reading) as db:
            z = ZFINStages()
            z.id = stage_id
            ref = db[DB_DATA][DB_STAGE_ONTOLOGY][z.id][0]
            z.obo_id = ref[DB_STAGE_OBO_ID]
            z.name = ref[DB_STAGE_NAME]
            z.start = ref[DB_STAGE_START]
            z.end = ref[DB_STAGE_END]
            return z

    def get_all_anatomy_items(self) -> list:
        """
        :return: All Anatomy items of the database
        """
        anatomy_items = []
        with h5py.File(DATABASE, self.__mode_reading) as db:
            item_ids = []
            item_ids.extend(db[DB_DATA][DB_ANATOMY].keys())
            for i in item_ids:
                k = self.get_anatomy_item(i)
                anatomy_items.append(k)
        return anatomy_items

    def get_all_stages(self) -> list:
        """
        :return: All stages of the database
        """
        items = []
        with h5py.File(DATABASE, self.__mode_reading) as db:
            item_ids = []
            item_ids.extend(db[DB_DATA][DB_STAGE_ONTOLOGY].keys())
            for i in item_ids:
                items.append(self.get_stage(i))
        return items
