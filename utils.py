import pandas as pd
import pyvista as pv
import mne


def visualize_digitalization_file_3d(
    filepath,
    always_visible=True,
    point_size=10,
    font_size=20,
    sphere_radius=0.005,
    distance_scaler = 0.001,
    background_color="white"
):
    """
    Plot digitized EEG electrode positions from a CSV file using PyVista.

    Parameters:
    - filepath (str or Path): Path to the digitization CSV file.
    - always_visible (bool): Whether to always show labels.
    - point_size (int): Point size for labels.
    - font_size (int): Font size for labels.
    - sphere_radius (float): Radius of electrode marker spheres.
    - background_color (str): Background color for the plotter.
    """
    df = pd.read_csv(filepath)

    # Extract channel names, replacing 33–64 with offset +8
    ch_names = df['#'].astype(str).tolist()


    # for i in range(33, 65):
    #     ch_names[i - 1] = str(i + 8)

    # Coordinates in meters
    coords_m = (df[['x', 'y', 'z']] / distance_scaler).values
    ch_pos = {name: coord for name, coord in zip(ch_names, coords_m)}

    # PyVista plot
    plotter = pv.Plotter()
    plotter.set_background(background_color)

    for name, coord in ch_pos.items():
        sphere = pv.Sphere(radius=sphere_radius, center=coord)
        plotter.add_mesh(sphere, color="blue")
        plotter.add_point_labels(
            [coord],
            [name],
            point_size=point_size,
            font_size=font_size,
            text_color='black',
            always_visible=always_visible
        )

    # Top-down view
    plotter.view_vector((0, 0, 1), viewup=(0, 1, 0))
    plotter.show()




def visualize_digitalization_file_2d(filepath, distance_scaler = 0.001, show_names=True):
    """
    Plot 2D EEG topomap from digitization file using MNE.

    Parameters:
    - filepath (str or Path): Path to the digitization CSV file.
    - show_names (bool): Whether to show channel names on the plot.
    """
    df = pd.read_csv(filepath)

    # Extract and correct channel names (add 8 for index 33–64)
    ch_names = df['#'].astype(str).tolist()

    # for i in range(33, 65):
    #     ch_names[i - 1] = str(i + 8)

    coords_m = (df[['x', 'y', 'z']] * distance_scaler).values  # mm → m
    ch_pos = {name: coord for name, coord in zip(ch_names, coords_m)}

    # Create montage and dummy info (sfreq can be any number, not used)
    montage = mne.channels.make_dig_montage(ch_pos=ch_pos, coord_frame='head')

    montage.plot(kind="topomap", show_names=show_names)









