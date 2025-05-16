import matplotlib.pyplot as plt
import mne

mne.viz.set_3d_backend('pyvistaqt')  # IMPORTANT!


available_montages = mne.channels.get_builtin_montages()
print(available_montages)


# Load the 64-channel EasyCap M10 montage
montage = mne.channels.make_standard_montage('biosemi256')

# Optionally inspect channel names
print(montage.ch_names)

# print how many channels
print("Number of channels:", len(montage.ch_names))

# Visualize in 2D
montage.plot(kind = "topomap", show_names=True)


# fig = montage.plot(kind="3d", show=False)  # 3D
# fig = fig.gca().view_init(azim=-90, elev=90)  # set view angle for tutorial
# plt.show()


