import mne
import os
import numpy as np
import matplotlib.pyplot as plt

def plot_butterfly_topomap_with_gfp(epochs, plot_range_ms=(-10, 200), times_ms=None, plot_title="Joint Plot", fig_dpi=100):
    """
    Plot the averaged evoked response as a joint plot (butterfly + topomaps)
    using MNE's plot_joint(), and let MNE plot the Global Field Power (GFP)
    automatically (which is essentially the same as GMFA).

    Parameters
    ----------
    epochs : mne.Epochs
        The epoched EEG data in MNE format.
    plot_range_ms : tuple of (float, float)
        The time range (in ms) to crop and plot.
    times_ms : list or None
        Time points (in ms) for the topomaps in the joint plot.
    plot_title : str
        Title for the joint plot.
    fig_dpi : int
        The resolution (dots per inch) of the output figure.
    """
    # Convert millisecond values to seconds
    tmin_sec = plot_range_ms[0] / 1000.0
    tmax_sec = plot_range_ms[1] / 1000.0

    # Crop the epochs to the desired time range (does not alter original data)
    epochs_cropped = epochs.copy().crop(tmin=tmin_sec, tmax=tmax_sec)

    # Average the epochs to get an Evoked object
    evoked = epochs_cropped.average()

    # Create the joint plot
    if times_ms:
        times_sec = [ms / 1000.0 for ms in times_ms]
        fig = evoked.plot_joint(
            times=times_sec,
            title=plot_title,
            show=False,
            ts_args=dict(gfp=True)  # Plot MNE's built-in GFP
        )
    else:
        fig = evoked.plot_joint(
            title=plot_title,
            show=False,
            ts_args=dict(gfp=True)  # Plot MNE's built-in GFP
        )

    # Set figure resolution
    fig.set_dpi(fig_dpi)

    # Show the final figure
    plt.show()

# -------------------------------------------------------------------------
# Usage example
# -------------------------------------------------------------------------
# Define paths and load the dataset
output_path = "../data/"
output_filename = "sub-01_ses-eegT1_final.set"
epochs = mne.read_epochs_eeglab(os.path.join(output_path, output_filename))

# Plot with custom figure resolution (e.g., 300 DPI for high-quality export)
plot_butterfly_topomap_with_gfp(
    epochs,
    plot_range_ms=(-10, 200),
    times_ms=[15, 30, 60, 100],
    plot_title="Evoked Response with Topomaps at 15, 30, 60, 100 ms",
    fig_dpi=500
)
