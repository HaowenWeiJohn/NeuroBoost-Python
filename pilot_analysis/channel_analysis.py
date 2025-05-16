import numpy
import matplotlib.pyplot as plt
import mne
from pathlib import Path
from dotenv import load_dotenv
import os

load_dotenv()  # Loads from default .env file

stimulus_event_id = {'Stimulus/A': 10001}

TMS_EEG_ROOT_DIR = os.getenv('TMS_EEG_ROOT_DIR')
EXPERIMENT_NAME = os.getenv('EXPERIMENT_NAME')
PARTICIPANT_ID = os.getenv('PARTICIPANT_ID')

# check if TMS_EEG_ROOT_DIR is set and if this directory exists
assert TMS_EEG_ROOT_DIR is not None, "TMS_EEG_ROOT_DIR is not set in the environment variables."
assert os.path.exists(TMS_EEG_ROOT_DIR), f"TMS_EEG_ROOT_DIR {TMS_EEG_ROOT_DIR} does not exist."

digitalization_file_path = Path(TMS_EEG_ROOT_DIR)   / EXPERIMENT_NAME / PARTICIPANT_ID / 'digitization'

assert digitalization_file_path.exists(), f"Digitalization file {digitalization_file_path} does not exist."

session = 'Pos10_80'

# session_path = Path(TMS_EEG_ROOT_DIR) / session / Path(session+'.vhdr')

session_path = Path(TMS_EEG_ROOT_DIR) / EXPERIMENT_NAME / PARTICIPANT_ID / session / Path(session+'.vhdr')


assert session_path.exists(), f"Session file {session_path} does not exist."

raw = mne.io.read_raw_brainvision(session_path, preload=True)

# Print the sampling frequency in Hz
print("Sampling frequency (Hz):", raw.info['sfreq'])

# Print the list of channel names as specified in the header
print("Channel names:", raw.info['ch_names'])





























