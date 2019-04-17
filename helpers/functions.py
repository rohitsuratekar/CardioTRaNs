"""
CardioTRaNs 2019
Author: Rohit Suratekar

All Helper functions
"""

import tkinter as tk

from helpers.parsers import *


def print_stages(start_time: int = None, end_time: int = None,
                 print_end: bool = False):
    """
    Prints developmental stages
    :param start_time: Starting Time (in hours)
    :param end_time: Ending Time (in hours)
    :param print_end: If True, prints end time as well
    :return: ZFIN ID (start time in hours) --ZFIN name
    """

    for s in parse_stage_ontology():
        if start_time is None:
            start_time = 0
        if end_time is None:
            end_time = 3000
        if end_time >= s.start >= start_time:
            if print_end:
                print("{}\t({}) -- ({}) {} ".format(s.id, s.start, s.end,
                                                    s.name))
            else:
                print("{}\t({}) --{} ".format(s.id, s.start, s.name))


def print_structure(search_item: str = None):
    items = []
    for s1 in parse_zfin_anatomy_items():
        s = s1  # type:ZFINAnatomy
        if search_item is not None:
            if search_item.lower().strip() in s.name.strip().lower():
                items.append(s)
        else:
            items.append(s)
    db = LocalDatabase()
    for x in items:
        print(db.get_anatomy_item(x.parent_id).name)


def run():
    root = tk.Tk()
    w = tk.Message(root, text="Hello", padx=30, pady=30)
    w.pack()
    root.mainloop()
