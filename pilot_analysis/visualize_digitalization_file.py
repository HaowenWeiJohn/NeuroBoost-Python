import numpy
import matplotlib.pyplot as plt
import mne
from pathlib import Path
from dotenv import load_dotenv
import os
import pandas as pd
import pandas as pd
import pyvista as pv
from pyvistaqt import BackgroundPlotter


load_dotenv()  # Loads from default .env file



TMS_EEG_ROOT_DIR = os.getenv('TMS_EEG_ROOT_DIR')
EXPERIMENT_NAME = os.getenv('EXPERIMENT_NAME')
PARTICIPANT_ID = os.getenv('PARTICIPANT_ID')

# check if TMS_EEG_ROOT_DIR is set and if this directory exists
assert TMS_EEG_ROOT_DIR is not None, "TMS_EEG_ROOT_DIR is not set in the environment variables."
assert os.path.exists(TMS_EEG_ROOT_DIR), f"TMS_EEG_ROOT_DIR {TMS_EEG_ROOT_DIR} does not exist."

digitalization_file_path = Path(TMS_EEG_ROOT_DIR)   / EXPERIMENT_NAME / PARTICIPANT_ID / 'digitization'

assert digitalization_file_path.exists(), f"Digitalization file {digitalization_file_path} does not exist."


df = pd.read_csv(digitalization_file_path)      # skip the header comment if present

ch_names = df['#'].astype(str).tolist()

# for channel name from 33 to 64, we add 8
for i in range(33, 65):
    ch_names[i-1] =  str(i+8)

coords_m = (df[['x', 'y', 'z']] / 1000).values  # convert mm→m; remove /1000 if already metres

ch_pos   = {name: coord for name, coord in zip(ch_names, coords_m)}



# Use blocking Plotter
plotter = pv.Plotter()
plotter.set_background("white")

# Plot spheres
for name, coord in ch_pos.items():
    sphere = pv.Sphere(radius=0.005, center=coord)
    plotter.add_mesh(sphere, color="blue")
    plotter.add_point_labels([coord], [name], point_size=10, font_size=20, text_color='black', always_visible=True)

# Top-down view
plotter.view_vector((0, 0, 1), viewup=(0, 1, 0))

# Blocking call to keep window open
plotter.show()



# montage = mne.channels.make_dig_montage(ch_pos=ch_pos,
#                                         coord_frame='head')  # 'head' is typical for Polhemus/Xensor
#
# # Visualize in 2D
# montage.plot(kind = "topomap", show_names=True)
#
#
# fig = montage.plot(kind="3d", show=False)  # 3D
# fig = fig.gca().view_init(azim=-90, elev=0)  # set view angle for tutorial
# plt.show()



