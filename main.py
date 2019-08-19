import os

from constants.system import PLOT_FOLDER
from visualization.atlas import run

# Create folder for saving plots
if not os.path.exists(PLOT_FOLDER):
    os.makedirs(PLOT_FOLDER)

run()
