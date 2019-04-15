"""
CardioTRaNs 2019
Author: Rohit Suratekar

All constants of this project (including database constants)
"""

DATABASE = "CardioTrans.hdf5"  # Name of the main HDF5 database
DATABASE_VERSION = 1
DATABASE_AUTHOR = "Rohit Suratekar"

DATA_FOLDER = "data/"  # Database Folder

# All database files
FILE_ANATOMY_ITEMS = "anatomy_item_2019.04.08.txt"
FILE_ANATOMY_RELATIONSHIP = "anatomy_relationship_2019.04.08.txt"
FILE_STAGE_ONTOLOGY = "stage_ontology_2019.04.08.txt"
FILE_WILD_TYPE_EXPRESSION = "wildtype-expression_fish_2019.04.08.txt"

# Database Groups
DB_DATA = "zfin_data"
DB_ANATOMY = "anatomy"
DB_ANATOMY_RELATIONSHIPS = "anatomy_relationships"
DB_STAGE_ONTOLOGY = "stage_ontology"
DB_WILD_TYPE_EXP = "wild_type_expression"

# Anatomy constants
DB_ANATOMY_ID = 0
DB_ANATOMY_NAME = 1
DB_ANATOMY_START_STAGE = 2
DB_ANATOMY_END_STAGE = 3

# Stage Constants
DB_STAGE_ID = 0
DB_STAGE_OBO_ID = 1
DB_STAGE_NAME = 2
DB_STAGE_START = 3
DB_STAGE_END = 4

# Meta-data attributes
DB_ATTRS_VERSION = "version"
DB_ATTRS_CREATION = "creation"
DB_ATTRS_AUTHOR = "author"
