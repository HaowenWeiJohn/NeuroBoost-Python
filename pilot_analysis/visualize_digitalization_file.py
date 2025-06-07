from utils import visualize_digitalization_file_3d, visualize_digitalization_file_2d
import os
from dotenv import load_dotenv
from pathlib import Path

load_dotenv()  # Loads from default .env file

stimulus_event_id = {'Stimulus/A': 10001}

TMS_EEG_ROOT_DIR = os.getenv('TMS_EEG_ROOT_DIR')
EXPERIMENT_NAME = os.getenv('EXPERIMENT_NAME')
PARTICIPANT_ID = os.getenv('PARTICIPANT_ID')

channel_names = [
    'Fp1', 'Fpz', 'Fp2',
    'AF7', 'AF3', 'AFZ', 'AF4', 'AF8',
    'F7', 'F3', 'F1', 'FZ', 'F2', 'F4', 'F8',
    'FT7', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'FT8',
    'C5', 'C3', 'C1', 'Cz', 'C2', 'C4', 'C6',
    'TP9',  'CP5', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6', 'TP10',
    'P9', 'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10',
    'PO9', 'PO7', 'PO3', 'POz', 'PO4', 'PO8', 'PO10',
    'O9', 'O1', 'Oz', 'O2', 'O10',
    'Iz'
]

digitalization_file_path = Path(TMS_EEG_ROOT_DIR) / EXPERIMENT_NAME / PARTICIPANT_ID / 'digitization.csv'

assert digitalization_file_path.exists(), f"Digitalization file {digitalization_file_path} does not exist."

visualize_digitalization_file_2d(digitalization_file_path, channel_names=channel_names)

visualize_digitalization_file_3d(digitalization_file_path, channel_names=channel_names)




