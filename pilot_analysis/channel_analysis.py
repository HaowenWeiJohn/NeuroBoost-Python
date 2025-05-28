import numpy
import matplotlib.pyplot as plt
import mne
from pathlib import Path
from dotenv import load_dotenv
import os

from mne import channel_type

from utils import visualize_digitalization_file_3d, visualize_digitalization_file_2d, \
    create_montage_from_digitization_file, visualize_epochs_channels

load_dotenv()  # Loads from default .env file

channel_names = ['EEG 001', 'EEG 002', 'EEG 003', 'EEG 004', 'EEG 005', 'EEG 006', 'EEG 007', 'EEG 008',
                'EEG 009', 'EEG 010', 'EEG 011', 'EEG 012', 'EEG 013', 'EEG 014', 'EEG 015', 'EEG 016',
                'EEG 017', 'EEG 018', 'EEG 019', 'EEG 020', 'EEG 021', 'EEG 022', 'EEG 023', 'EEG 024',
                'EEG 025', 'EEG 026', 'EEG 027', 'EEG 028', 'EEG 029', 'EEG 030', 'EEG 031', 'EEG 032',
                'EMG 001',
                'EEG 033', 'EEG 034', 'EEG 035', 'EEG 036', 'EEG 037', 'EEG 038', 'EEG 039', 'EEG 040',
                'EEG 041', 'EEG 042', 'EEG 043', 'EEG 044', 'EEG 045', 'EEG 046', 'EEG 047', 'EEG 048',
                'EEG 049', 'EEG 050', 'EEG 051', 'EEG 052', 'EEG 053', 'EEG 054', 'EEG 055', 'EEG 056',
                'EEG 057', 'EEG 058', 'EEG 059', 'EEG 060', 'EEG 061', 'EEG 062', 'EEG 063', 'EEG 064',
                ]

electrodes = [
    'Fp1', 'Fpz', 'Fp2',
    'AF7', 'AF3', 'AFZ', 'AF4', 'AF8',
    'F7', 'F3', 'F1', 'FZ', 'F2', 'F4', 'F8',
    'FT7', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'FT8',
    'C5', 'C3', 'C1', 'Cz', 'C2', 'C4', 'C6',
    'TP9', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6', 'TP10',
    'P9', 'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10',
    'PO9', 'PO7', 'PO3', 'POz', 'PO4', 'PO8', 'PO10',
    'O9', 'O1', 'Oz', 'O2', 'O10',
    'Iz'
]

# 32 eeg, 1 emg, 32 eeg, list
channel_types = ['eeg'] * 32 + ['emg'] + ['eeg'] * 32

stimulus_event_id = {'Stimulus/A': 10001}

baseline = (-0.2, 0)  # Baseline period (200 ms before the event)

TMS_EEG_ROOT_DIR = os.getenv('TMS_EEG_ROOT_DIR')
EXPERIMENT_NAME = os.getenv('EXPERIMENT_NAME')
PARTICIPANT_ID = os.getenv('PARTICIPANT_ID')

# check if TMS_EEG_ROOT_DIR is set and if this directory exists
assert TMS_EEG_ROOT_DIR is not None, "TMS_EEG_ROOT_DIR is not set in the environment variables."
assert os.path.exists(TMS_EEG_ROOT_DIR), f"TMS_EEG_ROOT_DIR {TMS_EEG_ROOT_DIR} does not exist."

digitalization_file_path = Path(TMS_EEG_ROOT_DIR)   / EXPERIMENT_NAME / PARTICIPANT_ID / 'digitization.csv'

assert digitalization_file_path.exists(), f"Digitalization file {digitalization_file_path} does not exist."

session = 'Pos10_80'

# session_path = Path(TMS_EEG_ROOT_DIR) / session / Path(session+'.vhdr')

session_path = Path(TMS_EEG_ROOT_DIR) / EXPERIMENT_NAME / PARTICIPANT_ID / Path(session+'.vhdr')


assert session_path.exists(), f"Session file {session_path} does not exist."

raw = mne.io.read_raw_brainvision(session_path, preload=True)

# Print the sampling frequency in Hz
print("Sampling frequency (Hz):", raw.info['sfreq'])

# Print the list of channel names as specified in the header
print("Channel names:", raw.info['ch_names'])

montage = create_montage_from_digitization_file(digitalization_file_path, channel_names= channel_names[0:32] + channel_names[33:65])

print("Montage channel names:", montage.ch_names)

old_channel_names = raw.ch_names

mapping = {old_channel_names[i]: channel_names[i] for i in range(len(old_channel_names))}

raw.rename_channels(mapping)

# set channel types
raw.set_channel_types({name: ctype for name, ctype in zip(channel_names, channel_types)})

# set montage
eeg_channel_names = [name for name, ctype in zip(channel_names, channel_types) if ctype == 'eeg']
raw.pick(eeg_channel_names).set_montage(montage)



# create events
events, event_id = mne.events_from_annotations(raw)


# create epochs
epochs = mne.Epochs(raw, events, event_id=stimulus_event_id,
                    tmin=-0.02, tmax=0.0,  # 1 second before and after the event
                    baseline=None,
                    preload=True
                    )



# evoked response
evoked = epochs.average()

# visualize the evoked response
evoked.plot(ylim=dict(eeg=[-300, 300]))





# # visualize epochs
# evoked.plot(ylim=dict(eeg=[-300, 300]))
#
browser = epochs.plot(
    n_epochs   = 3,      # whatever you like
    n_channels = 5,
    scalings   = 'auto',
    block      = True     # blocks script until you close the window
)












# print("Number of epochs:", len(epochs))



browser = raw.plot(
    n_channels = len(raw.ch_names),  # show every channel
    duration   = 20,                 # seconds visible in one screenful
    scalings   = 'auto',             # per‑channel autoscaling
    block      = True                # pause script until window closes
)