%% ---------------------------------------------------------------
clear; close all; clc;

%% 0. Choose where to look for the .env file  ---------------------
%    • Pass an explicit path, e.g.  '../.env'
%    • Leave empty ("") to fall back to current folder
ENV_PATH = fullfile('..', '.env');   % ← parametric

env = loadenv(ENV_PATH);             % loads only if file exists

%% 1. Read the variables (will be empty strings if .env missing)
ROOT_DIR        = getenv('TMS_EEG_ROOT_DIR');
EXPERIMENT_NAME = getenv('EXPERIMENT_NAME');
PARTICIPANT_ID  = getenv('PARTICIPANT_ID');

% 1a.  Basic sanity checks
assert(~isempty(ROOT_DIR),       'TMS_EEG_ROOT_DIR is not set.');
assert(isfolder(ROOT_DIR),       'Folder "%s" does not exist.', ROOT_DIR);

%% ---------------------------------------------------------------
%  2.  Build paths (MATLAB’s fullfile avoids slash/back-slash issues)
data_root = fullfile(ROOT_DIR, EXPERIMENT_NAME, PARTICIPANT_ID);

digitization_file = fullfile(data_root, 'digitization.csv');
assert(isfile(digitization_file), ...
       'Digitalization file "%s" does not exist.', digitization_file);

session = 'Pos10_80';              % ← adjust if you change session names


vhdr_file   = fullfile(data_root, [session '.vhdr']);
assert(isfile(vhdr_file), ...
       'Session file "%s" does not exist.', vhdr_file);



%% ---------------------------------------------------------------
%  3.  Launch EEGLAB and import the BrainVision recording
[ALLEEG,EEG,CURRENTSET,ALLCOM] = eeglab;

EEG = pop_loadbv(data_root, [session '.vhdr']);

EEG = eeg_checkset(EEG);   % housekeeping + GUI refresh
EEG.setname = sprintf('%s_%s_%s', EXPERIMENT_NAME, PARTICIPANT_ID, session);

% OPTIONAL: save an immediate .set copy so the rest of the pipeline
% can reload quickly if you crash or want to inspect intermediate stages
%EEG = pop_saveset(EEG, ...
%         'filename', [EEG.setname '_raw.set'], ...
%         'filepath',  data_root);

%PIPELINE_NAME = "TESA_ICA";


%% -----------------------------------------------------------------


% replace the channel names

%channel_names = { ...
%     'EEG 001','EEG 002','EEG 003','EEG 004','EEG 005','EEG 006','EEG 007','EEG 008', ...
%     'EEG 009','EEG 010','EEG 011','EEG 012','EEG 013','EEG 014','EEG 015','EEG 016', ...
%     'EEG 017','EEG 018','EEG 019','EEG 020','EEG 021','EEG 022','EEG 023','EEG 024', ...
%     'EEG 025','EEG 026','EEG 027','EEG 028','EEG 029','EEG 030','EEG 031','EEG 032', ...
%     'EMG 001',  ...                            % 33rd channel
%     'EEG 033','EEG 034','EEG 035','EEG 036','EEG 037','EEG 038','EEG 039','EEG 040', ...
%     'EEG 041','EEG 042','EEG 043','EEG 044','EEG 045','EEG 046','EEG 047','EEG 048', ...
%     'EEG 049','EEG 050','EEG 051','EEG 052','EEG 053','EEG 054','EEG 055','EEG 056', ...
%     'EEG 057','EEG 058','EEG 059','EEG 060','EEG 061','EEG 062','EEG 063','EEG 064' };

% remove EMG channels
EEG = pop_select(EEG, 'nochannel', {'Input 33'});

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


assert(numel(EEG.chanlocs) == numel(channel_names), ...
       'Dataset has %d channels but you supplied %d names.', ...
       numel(EEG.chanlocs), numel(channel_names));

for k = 1:numel(channel_names)
    EEG.chanlocs(k).labels = channel_names{k};
end

EEG = eeg_checkset(EEG);   % housekeeping + GUI refresh
EEG = pop_chanedit(EEG, 'lookup','standard-10-5-cap385.elp');












