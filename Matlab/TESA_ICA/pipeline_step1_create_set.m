clear; close all; clc;

% Define the pipeline name for tracking/logging purposes
PIPELINE_NAME = 'TESA_ICA';   
STEP_NAME = 'step1_create_set';

%% 0. Load environment variables from .env file ----------------------------
%    • Use fullfile for cross-platform compatibility (avoids slash issues)
%    • ENV_PATH can be changed to point to any specific location
ENV_PATH = fullfile('..', '.env');   % ← parametric
env = loadenv(ENV_PATH);             % Load environment only if the file exists

%% 1. Read experiment parameters from environment variables ----------------
%    These variables are essential for locating data
ROOT_DIR        = getenv('TMS_EEG_ROOT_DIR');
EXPERIMENT_NAME = getenv('EXPERIMENT_NAME');
PARTICIPANT_ID  = getenv('PARTICIPANT_ID');

% Sanity checks to ensure required paths are available
assert(~isempty(ROOT_DIR),       'TMS_EEG_ROOT_DIR is not set.');
assert(isfolder(ROOT_DIR),       'Folder "%s" does not exist.', ROOT_DIR);

%% 2. Construct data paths -------------------------------------------------
%    Build the full path to the participant's data folder
data_root = fullfile(ROOT_DIR, EXPERIMENT_NAME, PARTICIPANT_ID, 'data');

% %    Path to digitization file (e.g., electrode positions)
% digitization_file = fullfile(data_root, 'digitization.csv');
% assert(isfile(digitization_file), ...
%        'Digitalization file "%s" does not exist.', digitization_file);

%    Define the session name (e.g., stimulation condition)
session = 'Pos10_80';  % ← adjust if you change session names

%    Construct the full path to the BrainVision header file (.vhdr)
vhdr_file = fullfile(data_root, [session '.vhdr']);
assert(isfile(vhdr_file), ...
       'Session file "%s" does not exist.', vhdr_file);

%% 3. Load EEG data into EEGLAB --------------------------------------------
%    Launch EEGLAB and import BrainVision data
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
EEG = pop_loadbv(data_root, [session '.vhdr']);
EEG = eeg_checkset(EEG);  % Ensure internal consistency and update GUI
EEG.setname = sprintf('%s_%s_%s', EXPERIMENT_NAME, PARTICIPANT_ID, session);

%% 4. Preprocessing: Clean up channel list and apply standard names --------
%    Remove EMG or extra input channels (e.g., 'Input 33' from Biosemi)
EEG = pop_select(EEG, 'nochannel', {'Input 33'});

%    Rename the remaining EEG channels to match standard 10-5 layout
channel_names = {
    'Fp1', 'Fpz', 'Fp2', ...
    'AF7', 'AF3', 'AFZ', 'AF4', 'AF8', ...
    'F7', 'F3', 'F1', 'FZ', 'F2', 'F4', 'F8', ...
    'FT7', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'FT8', ...
    'C5', 'C3', 'C1', 'Cz', 'C2', 'C4', 'C6', ...
    'TP9', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6', 'TP10', ...
    'P9', 'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10', ...
    'PO9', 'PO7', 'PO3', 'POz', 'PO4', 'PO8', 'PO10', ...
    'O9', 'O1', 'Oz', 'O2', 'O10', ...
    'Iz'};

%    Ensure the number of provided labels matches the number of channels
assert(numel(EEG.chanlocs) == numel(channel_names), ...
       'Dataset has %d channels but you supplied %d names.', ...
       numel(EEG.chanlocs), numel(channel_names));

%    Assign new channel labels
for k = 1:numel(channel_names)
    EEG.chanlocs(k).labels = channel_names{k};
end

%    Perform final consistency checks and lookup standard 10-5 positions
EEG = eeg_checkset(EEG);
EEG = pop_chanedit(EEG, 'lookup','standard-10-5-cap385.elp');

%% 5. Display basic EEG dataset information --------------------------------
fprintf('\n--- EEG Dataset Overview ---\n');
fprintf('Set name         : %s\n', EEG.setname);
fprintf('Number of channels: %d\n', EEG.nbchan);
fprintf('Number of data points per channel: %d\n', EEG.pnts);
fprintf('Sampling rate    : %.2f Hz\n', EEG.srate);
fprintf('Duration         : %.2f seconds\n', EEG.pnts / EEG.srate);
fprintf('Number of epochs : %d\n', EEG.trials);
fprintf('Data dimensions  : [%d channels x %d time points x %d epochs]\n', ...
        size(EEG.data, 1), size(EEG.data, 2), size(EEG.data, 3));
fprintf('Channel labels   : %s ... %s\n', EEG.chanlocs(1).labels, EEG.chanlocs(end).labels);


output = fullfile(ROOT_DIR, EXPERIMENT_NAME, PARTICIPANT_ID, 'output');


%% 6. Save the dataset for this pipeline step --------------------------------
out_dir = char( fullfile(ROOT_DIR, EXPERIMENT_NAME, PARTICIPANT_ID, ...
                         'output', PIPELINE_NAME) );   % <-- char(...)

if ~exist(out_dir, 'dir'); mkdir(out_dir); end

save_fname = char( sprintf('%s-%s-%s-%s.set', ...
                           EXPERIMENT_NAME, PARTICIPANT_ID, session, STEP_NAME) );

EEG.setname  = save_fname;
EEG.filename = save_fname;
EEG.filepath = out_dir;

pop_saveset(EEG, 'filename', save_fname, 'filepath', out_dir);

fprintf('\nData saved to  : %s\n', fullfile(out_dir, save_fname));
fprintf('Pipeline stage : %s\n', STEP_NAME);