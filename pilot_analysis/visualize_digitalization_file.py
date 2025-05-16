from utils import visualize_digitalization_file_3d, visualize_digitalization_file_2d
import os
from dotenv import load_dotenv
from pathlib import Path

load_dotenv()  # Loads from default .env file

stimulus_event_id = {'Stimulus/A': 10001}

TMS_EEG_ROOT_DIR = os.getenv('TMS_EEG_ROOT_DIR')
EXPERIMENT_NAME = os.getenv('EXPERIMENT_NAME')
PARTICIPANT_ID = os.getenv('PARTICIPANT_ID')


digitalization_file_path = Path(TMS_EEG_ROOT_DIR)   / EXPERIMENT_NAME / PARTICIPANT_ID / 'digitization'

assert digitalization_file_path.exists(), f"Digitalization file {digitalization_file_path} does not exist."


visualize_digitalization_file_2d(digitalization_file_path)


visualize_digitalization_file_3d(digitalization_file_path)




