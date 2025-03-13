import os
import numpy as np
import mne
import matplotlib.pyplot as plt

###############################################################################
#                       USER-PARAMETERS AT THE TOP
###############################################################################
output_path = "../data/"
file_pattern_base = "sub-01_ses-eegT1"  # e.g. "sub-01_ses-eegT1_XXepochs.set"

start_epoch_count = 10
increment_epoch_count = 5
stop_epoch_count = 120

# The reference file to load (the "gold standard")
# Must exist in the format: sub-01_ses-eegT1_{reference_count}epochs.set
reference_count = 120

# Time window (ms) for comparing subsets to the reference
compare_time_window_ms = (-1, 100)

# Figure DPI
fig_dpi = 120

###############################################################################
#                      HELPER FUNCTIONS
###############################################################################
def crop_evoked(evoked_obj, tmin_ms, tmax_ms):
    """
    Return a copy of `evoked_obj` cropped to the specified window in milliseconds.
    """
    tmin_s = tmin_ms / 1000.0
    tmax_s = tmax_ms / 1000.0
    return evoked_obj.copy().crop(tmin=tmin_s, tmax=tmax_s)


def rms_diff(evoked_a, evoked_b):
    """
    Compute RMS difference between two evoked objects (assuming same channels & times).
    We'll flatten channels x time, then do RMS of differences.
    """
    data_a = evoked_a.data  # shape: (n_channels, n_times)
    data_b = evoked_b.data
    diff = data_a - data_b
    return np.sqrt(np.mean(diff**2))


def compute_correlation(evoked_a, evoked_b):
    """
    Compute Pearson correlation between flattened data of two evokeds.
    """
    arr_a = evoked_a.data.ravel()
    arr_b = evoked_b.data.ravel()
    r = np.corrcoef(arr_a, arr_b)[0, 1]
    return r

###############################################################################
#                           MAIN SCRIPT
###############################################################################
def main():
    # ------------------------------------------------------------------
    # 1) Load reference file once, create reference Evoked
    # ------------------------------------------------------------------
    ref_fname = f"{file_pattern_base}_{reference_count}epochs.set"
    ref_fullpath = os.path.join(output_path, ref_fname)
    if not os.path.isfile(ref_fullpath):
        raise FileNotFoundError(f"Reference file not found: {ref_fullpath}")

    print(f"Loading reference file: {ref_fullpath}")
    ref_epochs = mne.read_epochs_eeglab(ref_fullpath, verbose=False)
    ref_evoked_full = ref_epochs.average()
    ref_evoked_cropped = crop_evoked(ref_evoked_full,
                                     compare_time_window_ms[0],
                                     compare_time_window_ms[1])

    # Prepare data structures for results
    epoch_counts = range(start_epoch_count, stop_epoch_count+1, increment_epoch_count)
    rms_results = []
    corr_results = []

    # ------------------------------------------------------------------
    # 2) Loop over desired epoch counts, load them one by one, compare
    # ------------------------------------------------------------------
    for nEpochs in epoch_counts:
        fname = f"{file_pattern_base}_{nEpochs}epochs.set"
        fullpath = os.path.join(output_path, fname)

        if not os.path.isfile(fullpath):
            print(f"File not found: {fullpath} (skipping).")
            rms_results.append(np.nan)
            corr_results.append(np.nan)
            continue

        print(f"Loading subset file: {fullpath}")
        # load just this subset
        epochs_subset = mne.read_epochs_eeglab(fullpath, verbose=False)

        # Average to get Evoked
        evoked_subset_full = epochs_subset.average()
        # Crop to the same time window
        evoked_subset_crop = crop_evoked(evoked_subset_full,
                                         compare_time_window_ms[0],
                                         compare_time_window_ms[1])

        # If it's the reference file (nEpochs == reference_count),
        # RMS=0, correlation=1
        if nEpochs == reference_count:
            this_rms = 0.0
            this_corr = 1.0
        else:
            this_rms = rms_diff(evoked_subset_crop, ref_evoked_cropped)
            this_corr = compute_correlation(evoked_subset_crop, ref_evoked_cropped)

        rms_results.append(this_rms)
        corr_results.append(this_corr)

    # ------------------------------------------------------------------
    # 3) Plot results (RMS difference vs nEpochs)
    # ------------------------------------------------------------------
    plt.figure(figsize=(6,4), dpi=fig_dpi)
    plt.plot(epoch_counts, rms_results, 'o-', label="RMS Difference")
    plt.xlabel("Number of Epochs")
    plt.ylabel("RMS Difference (vs reference)")
    plt.title(f"Evoked Stability\nComparison Window: {compare_time_window_ms} ms")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # ------------------------------------------------------------------
    # 4) Plot correlation vs nEpochs
    # ------------------------------------------------------------------
    plt.figure(figsize=(6,4), dpi=fig_dpi)
    plt.plot(epoch_counts, corr_results, 'o-', color='green', label="Correlation")
    plt.xlabel("Number of Epochs")
    plt.ylabel("Correlation (vs reference)")
    plt.title(f"Evoked Stability\nComparison Window: {compare_time_window_ms} ms")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    print("Done! Inspect the two plots to see how stable the evoked becomes as the number of epochs increases.")

# Run it
if __name__ == "__main__":
    main()
