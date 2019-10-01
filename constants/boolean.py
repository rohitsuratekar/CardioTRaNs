#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
#  All constants related to the Boolean analysis

from constants.colors import *

SYMBOL_OR = COLOR_BLUE + "OR" + COLOR_END
SYMBOL_AND = COLOR_BLUE + "AND" + COLOR_END
SYMBOL_NOT = COLOR_BLUE + "NOT" + COLOR_END
SYMBOL_COPY = "COPY"

INTERESTED_GENES = ["nkx2.5", "gata4", "hand2", "hey2", "gata2a", "zfpm1",
                    "gata3", "smarca4a", "tbx5a", "mef2ca", "nkx2.3", "foxh1",
                    "gata2b", "sqstm1", "tbx5b", "zfpm2a", "gata5"]

BASE_GENES = ["nkx2.5", "gata4", "tbx5a", "mef2ca", "hand2"]

CONTROL_GENES = ["neurog1", "ins", "arr3b"]
