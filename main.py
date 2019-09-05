import os

from constants.system import PLOT_FOLDER
from helpers.parsers.ngs import run

# Create folder for saving plots
if not os.path.exists(PLOT_FOLDER):
    os.makedirs(PLOT_FOLDER)

run()
