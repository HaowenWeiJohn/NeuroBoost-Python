import pandas as pd
import pyvista as pv
import mne
import os

import matplotlib.pyplot as plt
import mne
import numpy as np
import scipy
import os
import json
from types import SimpleNamespace
import os
import re

from mne.viz import plot_topomap

import torch



def visualize_digitalization_file_3d(
    filepath,
    channel_names=None,
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
    - channel_names (list of str or None): Optional list of channel names to use as labels.
                                           If None, uses the '#' column in the file.
    - always_visible (bool): Whether to always show labels.
    - point_size (int): Point size for labels.
    - font_size (int): Font size for labels.
    - sphere_radius (float): Radius of electrode marker spheres.
    - background_color (str): Background color for the plotter.
    """
    df = pd.read_csv(filepath)

    # Extract channel names, replacing 33–64 with offset +8
    if channel_names is None:
        ch_names = df['#'].astype(str).tolist()
    else:
        if len(channel_names) != len(df):
            raise ValueError("Provided channel_names list must match the number of coordinate rows.")
        ch_names = channel_names

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




def visualize_digitalization_file_2d(filepath, channel_names=None, distance_scaler = 0.001, show_names=True):
    """
    Plot 2D EEG topomap from digitization file using MNE.

    Parameters:
    - filepath (str or Path): Path to the digitization CSV file.
    - channel_names (list of str or None): Optional list of channel names to use as labels.
                                           If None, uses the '#' column in the file.
    - show_names (bool): Whether to show channel names on the plot.
    """
    df = pd.read_csv(filepath)

    # Extract and correct channel names (add 8 for index 33–64)
    if channel_names is None:
        ch_names = df['#'].astype(str).tolist()
    else:
        if len(channel_names) != len(df):
            raise ValueError("Provided channel_names list must match the number of coordinate rows.")
        ch_names = channel_names

    coords_m = (df[['x', 'y', 'z']] * distance_scaler).values  # mm → m
    ch_pos = {name: coord for name, coord in zip(ch_names, coords_m)}

    # Create montage and dummy info (sfreq can be any number, not used)
    montage = mne.channels.make_dig_montage(ch_pos=ch_pos, coord_frame='head')

    montage.plot(kind="topomap", show_names=show_names)

def create_montage_from_digitization_file(filepath, channel_names=None, distance_scaler=0.001):
    """
    Create an MNE DigMontage object from a digitization file.

    Parameters:
    - filepath (str or Path): Path to the digitization CSV file.
    - channel_names (list of str or None): Optional list of channel names.
                                           If None, uses the '#' column in the file.
    - distance_scaler (float): Scaling factor for coordinates (default assumes input is in mm).

    Returns:
    - montage (mne.channels.DigMontage): The constructed montage object.
    """
    df = pd.read_csv(filepath)

    if channel_names is None:
        ch_names = df['#'].astype(str).tolist()
    else:
        if len(channel_names) != len(df):
            raise ValueError("Provided channel_names list must match the number of coordinate rows.")
        ch_names = channel_names

    coords = (df[['x', 'y', 'z']] * distance_scaler).values
    ch_pos = {name: coord for name, coord in zip(ch_names, coords)}

    montage = mne.channels.make_dig_montage(ch_pos=ch_pos, coord_frame='head')
    return montage


def reference_eeg_raw(raw, ref_channels='average', projection=False):
    """
    This function references the EEG data in a raw MNE object.

    Parameters:
    raw (mne.io.Raw): The raw MNE object containing the EEG data.
    ref_channels (str): The reference channel(s) to use for re-referencing. Default is 'average'.
    projection (bool): Whether to use the EEG projections for re-referencing. Default is False.

    Returns:
    mne.io.Raw: The re-referenced raw MNE object.
    """
    raw, _ = mne.set_eeg_reference(raw, ref_channels=ref_channels, projection=projection)
    return raw


def bandpass_filter_raw(raw, low_cutoff, high_cutoff, n_jobs=1, picks='eeg'):
    """
    This function applies a bandpass filter to the EEG data in a raw MNE object.

    Parameters:
    raw (mne.io.Raw): The raw MNE object containing the EEG data.
    low_cutoff (float): The low cutoff frequency for the bandpass filter.
    high_cutoff (float): The high cutoff frequency for the bandpass filter.
    n_jobs (int): The number of parallel jobs to run. Default is 1.

    Returns:
    mne.io.Raw: The raw MNE object with the bandpass filter applied.
    """
    raw.filter(low_cutoff, high_cutoff, n_jobs=n_jobs, picks=picks)
    return raw

def notch_filter_raw(raw, freqs, n_jobs=1, picks='eeg'):
    """
    This function applies a notch filter to the EEG data in a raw MNE object.

    Parameters:
    raw (mne.io.Raw): The raw MNE object containing the EEG data.
    freqs (list): The frequencies to notch filter.
    n_jobs (int): The number of parallel jobs to run. Default is 1.

    Returns:
    mne.io.Raw: The raw MNE object with the notch filter applied.
    """
    raw.notch_filter(freqs, picks=picks, n_jobs=n_jobs)
    return raw

def resample_raw(raw, sfreq, n_jobs=1):
    """
    This function resamples the EEG data in a raw MNE object.

    Parameters:
    raw (mne.io.Raw): The raw MNE object containing the EEG data.
    sfreq (int): The new sampling frequency to resample to.
    n_jobs (int): The number of parallel jobs to run. Default is 1.

    Returns:
    mne.io.Raw: The raw MNE object resampled to the new sampling frequency.
    """
    raw.resample(sfreq, n_jobs=n_jobs)
    return raw

def identify_bad_channels_and_interpolate(raw, thresh1=None, thresh2=None, proportion=0.3, picks='eeg'):

    """
    This function identifies bad channels in the EEG data and interpolates them.

    Parameters:
    raw (mne.io.Raw): The raw MNE object containing the EEG data.
    thresh1 (float): The threshold for identifying bad channels based on the median. Default is None.
    thresh2 (float): The threshold for identifying bad channels based on the absolute value. Default is None.
    proportion (float): The proportion of bad channels to interpolate. Default is 0.3.

    Returns:
    mne.io.Raw: The raw MNE object with bad channels interpolated
    """

    data = raw.get_data(picks=picks)
    # We found that the data shape of epochs is 3 dims
    # print(data.shape)
    if len(data.shape) > 2:
        data = np.squeeze(data)
    Bad_chns = []
    value = 0
    # Delete the much larger point
    if thresh1 is not None:
        md = np.median(np.abs(data))
        value = np.where(np.abs(data) > (thresh1 * md), 0, 1)
    if thresh2 is not None:
        value = np.where(np.abs(data) > thresh2, 0, 1)
    # Use the standard to pick out the bad channels
    Bad_chns = np.argwhere(np.mean(1 - value, axis=1) > proportion).flatten()
    if Bad_chns.size > 0:
        raw.info['bads'].extend([raw.ch_names[bad] for bad in Bad_chns])
        print('Bad channels: ', raw.info['bads'])
        raw.interpolate_bads()
    else:
        print('No bad channel currently')
    return raw





def visualize_epochs_channels(epochs, event_groups, colors, picks, tmin_vis, tmax_vis, title='', out_dir=None, verbose='INFO', fig_size=(12.8, 7.2),
                              is_plot_timeseries=True):
    """
    Visualize EEG epochs for different event types and channels.

    Args:
        epochs (mne.Epochs): The EEG epochs to visualize.
        event_groups (dict): A dictionary mapping event names to lists of event IDs. Only events in these groups will be plotted.
        colors (dict): A dictionary mapping event names to colors to use for plotting.
        picks (list): A list of EEG channels to plot.
        title (str, optional): The title to use for the plot. Default is an empty string.
        out_dir (str, optional): The directory to save the plot to. If None, the plot will be displayed on screen. Default is None.
        verbose (str, optional): The verbosity level for MNE. Default is 'INFO'.
        fig_size (tuple, optional): The size of the figure in inches. Default is (12.8, 7.2).
        is_plot_timeseries (bool, optional): Whether to plot the EEG data as a timeseries. Default is True.

    Returns:
        None

    Raises:
        None

    """

    # Set the verbosity level for MNE
    mne.set_log_level(verbose=verbose)

    # Set the figure size for the plot
    plt.rcParams["figure.figsize"] = fig_size

    # Plot each EEG channel for each event type
    if is_plot_timeseries:
        for ch in picks:
            for event_name, events in event_groups.items():
                try:
                    # Get the EEG data for the specified event type and channel
                    y = epochs.crop(tmin_vis, tmax_vis)[event_name].pick([ch]).get_data().squeeze(1)
                except KeyError:  # meaning this event does not exist in these epochs
                    continue
                y_mean = np.mean(y, axis=0)
                y1 = y_mean + scipy.stats.sem(y, axis=0)  # this is the upper envelope
                y2 = y_mean - scipy.stats.sem(y, axis=0)

                time_vector = np.linspace(tmin_vis, tmax_vis, y.shape[-1])

                # Plot the EEG data as a shaded area
                plt.fill_between(time_vector, y1, y2, where=y2 <= y1, facecolor=colors[event_name], interpolate=True,
                                 alpha=0.5)
                plt.plot(time_vector, y_mean, c=colors[event_name], label='{0}, N={1}'.format(event_name, y.shape[0]))

            # Set the labels and title for the plot
            plt.xlabel('Time (sec)')
            plt.ylabel('BioSemi Channel {0} (μV), shades are SEM'.format(ch))
            plt.legend()
            plt.title('{0} - Channel {1}'.format(title, ch))

            # Save or show the plot
            if out_dir:
                plt.savefig(os.path.join(out_dir, '{0} - Channel {1}.png'.format(title, ch)))
                plt.clf()
            else:
                plt.show()



def visualize_topomap_windows(t_min, t_max, split_window, data, time_vector, info, fig_title=None):
    """
    Plots the topomap for specified time windows using a data array, time vector, and MNE info object.

    Parameters:
    - t_min: float, Start time in seconds.
    - t_max: float, End time in seconds.
    - split_window: float, Size of each time window in seconds.
    - data: numpy.ndarray, The data array of shape (n_channels, n_times).
    - time_vector: numpy.ndarray, The time vector corresponding to the data points.
    - info: mne.Info, The MNE info object containing channel information.
    - fig_title: str, Optional title for the entire figure.

    Raises:
    - ValueError: If t_min or t_max is out of the time range of the time_vector.
    """

    # Check if t_min and t_max are within the range of the time_vector
    if not (time_vector[0] <= t_min <= time_vector[-1]):
        raise ValueError(f"t_min ({t_min}) is out of range. The available time range is "
                         f"{time_vector[0]} to {time_vector[-1]} seconds.")

    if not (time_vector[0] <= t_max <= time_vector[-1]):
        raise ValueError(f"t_max ({t_max}) is out of range. The available time range is "
                         f"{time_vector[0]} to {time_vector[-1]} seconds.")

    # Calculate the number of windows
    num_windows = int((t_max - t_min) / split_window)

    # Create subplots
    fig, axs = plt.subplots(1, num_windows, figsize=(22, 5), sharey=True)

    for i in range(num_windows):
        start_time = t_min + i * split_window
        end_time = start_time + split_window

        # Get the indices for the time window
        time_mask = (time_vector >= start_time) & (time_vector < end_time)

        # Average the data over the time window
        data_mean = data[:, time_mask].mean(axis=1)

        # Plot the topomap
        im, cn = plot_topomap(
            data_mean,
            info,
            axes=axs[i],
            res=512,
            show=False
        )

        # Set the title for each subplot
        axs[i].set_title(f"{int(start_time * 1e3)}-{int(end_time * 1e3)} ms")

    # Add the figure title if provided
    if fig_title:
        fig.suptitle(fig_title, fontsize=16)

    plt.tight_layout()
    plt.show()





def count_dirs_starting_with(path, prefix):
    """
    This function counts the number of directories or files in the specified path
    that start with the given prefix.

    Parameters:
    path (str): The path to the directory containing the files or subdirectories.
    prefix (str): The prefix to match at the beginning of the directories or files.

    Returns:
    int: The number of directories or files that start with the prefix.
    list: A list of directories or files that match the prefix.
    """
    # Get a list of directories or files that start with the given prefix
    matching_items = [d for d in os.listdir(path) if d.startswith(prefix)]

    # Count the number of matching directories or files
    num_matching = len(matching_items)

    return num_matching, matching_items



def load_json_file(file_path, simple_name_space=False):
    """
    This function loads a JSON file and returns the data as a dictionary.

    Parameters:
    file_path (str): The path to the JSON file.
    simple_name_space (bool): Whether to use a simple namespace for the JSON data.

    Returns:
    dict: The JSON data as a dictionary.
    """
    with open(file_path, 'r') as file:
        if simple_name_space:
            data = json.load(file, object_hook=lambda d: SimpleNamespace(**d))
        else:
            data = json.load(file)

    return data


def count_bids_runs(subject_dir):
    """
    This function counts the number of runs in the EEG data for a given subject.

    Parameters:
    subject_dir (str): The path to the subject's EEG directory.

    Returns:
    int: The number of runs for this subject.
    """
    # Get all the files in the subject's EEG directory
    eeg_files = os.listdir(subject_dir)

    # Use a set to store the unique run numbers found
    runs = set()

    # Regex pattern to match 'run-X' in the filenames
    run_pattern = re.compile(r'run-(\d+)')

    # Check each file and find matches for the 'run-X' pattern
    for file in eeg_files:
        match = run_pattern.search(file)
        if match:
            runs.add(match.group(1))  # Add the run number to the set

    # The number of unique runs found
    num_runs = len(runs)

    return num_runs, sorted(runs)



def display_directory_structure(start_path, target_file=None, indent_level=0):
    try:
        if not os.path.exists(start_path):
            raise FileNotFoundError(f"Directory not found: {start_path}")

        for item in os.listdir(start_path):
            item_path = os.path.join(start_path, item)
            indent = "    " * indent_level
            if os.path.isdir(item_path):
                print(f"{indent}- {item}/")
                display_directory_structure(item_path, target_file, indent_level + 1)
            else:
                if item == target_file:
                    # Highlight the target file with asterisks
                    print(f"{indent}* {item} *")
                else:
                    print(f"{indent}- {item}")
    except FileNotFoundError as fnf_error:
        print(fnf_error)
    except PermissionError:
        print(f"Permission Denied: {start_path}")


def interaug(timg, label, batch_size):
    aug_data = []
    aug_label = []
    for cls4aug in range(4):
        cls_idx = np.where(label == cls4aug)
        tmp_data = timg[cls_idx]
        tmp_label = label[cls_idx]

        tmp_aug_data = np.zeros((int(batch_size / 4), 1, 22, 1000))
        for ri in range(int(batch_size / 4)):
            for rj in range(8):
                rand_idx = np.random.randint(0, tmp_data.shape[0], 8)
                tmp_aug_data[ri, :, :, rj * 125:(rj + 1) * 125] = tmp_data[rand_idx[rj], :, :,
                                                                  rj * 125:(rj + 1) * 125]

        aug_data.append(tmp_aug_data)
        aug_label.append(tmp_label[:int(batch_size / 4)])
    aug_data = np.concatenate(aug_data)
    aug_label = np.concatenate(aug_label)
    aug_shuffle = np.random.permutation(len(aug_data))
    aug_data = aug_data[aug_shuffle, :, :]
    aug_label = aug_label[aug_shuffle]

    aug_data = torch.from_numpy(aug_data).cuda()
    aug_data = aug_data.float()
    aug_label = torch.from_numpy(aug_label).cuda()
    aug_label = aug_label.long()
    return aug_data, aug_label



