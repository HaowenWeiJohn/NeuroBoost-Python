import mne
import os
import numpy as np
import matplotlib.pyplot as plt

def plot_butterfly_topomap_with_gfp(epochs, plot_range_ms=(-10, 200), times_ms=None, plot_title="Joint Plot", fig_dpi=100):
    """
    Plot the averaged evoked response as a joint plot (butterfly + topomaps)
    using MNE's plot_joint(), and let MNE plot the Global Field Power (GFP)
    automatically.

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

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object with the joint plot.
    """
    # Convert millisecond values to seconds
    tmin_sec = plot_range_ms[0] / 1000.0
    tmax_sec = plot_range_ms[1] / 1000.0

    # Crop the epochs to the desired time range (does not alter original data)
    epochs_cropped = epochs.copy().crop(tmin=tmin_sec, tmax=tmax_sec)

    # Average the epochs to get an Evoked object
    evoked = epochs_cropped.average()

    # Create the joint plot (with topomap time points if provided)
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

    # Set the figure resolution
    fig.set_dpi(fig_dpi)
    return fig

# Define parameters
output_path = "../data/"
output_filename_base = "sub-01_ses-eegT1"  # Base name for files

# User-defined loop parameters
start_epoch_count = 10    # Start from this number of epochs
increment_epoch_count = 5 # Step size (increment)
stop_epoch_count = 120    # Stop at this number of epochs

# Generate list of epoch counts based on user-defined parameters
epoch_counts = list(range(start_epoch_count, stop_epoch_count + 1, increment_epoch_count))

# Define the time range and topomap time points
plot_range_ms = (-10, 200)  # Time window for the plot (in ms)
times_ms = [15, 30, 60, 100]  # Time points for topomap visualization (in ms)

# Loop through different epoch counts and generate the plots
for n_epochs in epoch_counts:
    # Construct the filename dynamically
    filename = f"{output_filename_base}_{n_epochs}epochs.set"
    filepath = os.path.join(output_path, filename)

    # Check if the file exists before loading
    if not os.path.exists(filepath):
        print(f"Skipping {n_epochs} epochs - File not found: {filename}")
        continue

    # Load epochs from the EEGLAB .set file
    epochs = mne.read_epochs_eeglab(filepath, verbose=False)

    # Generate the butterfly + topomap plot and obtain the figure object
    fig = plot_butterfly_topomap_with_gfp(
        epochs,
        plot_range_ms=plot_range_ms,
        times_ms=times_ms,
        plot_title=f"Evoked Response ({n_epochs} Epochs) with Topomaps at {times_ms} ms",
        fig_dpi=500
    )

    # Save the figure with a descriptive filename
    save_path = os.path.join(output_path, f"evoked_{n_epochs}epochs.png")
    fig.savefig(save_path, dpi=500)
    plt.close(fig)  # Close the figure to free memory

    print(f"Finished plotting for {n_epochs} epochs. Saved at: {save_path}")

print("All epoch plots generated.")
