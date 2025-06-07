clear; close all; clc;

% ===========================================================================
% DATA INSPECTION AND INITIAL PROCESSING SCRIPT
% ===========================================================================
% This script performs initial processing on multiple TMS-EEG files:
% 1. Load EEG data from .vhdr files
% 2. Remove bad channels using joint probability
% 3. Find TMS pulses
% 4. Epoch data around TMS pulses
% 5. Save datasets and generate visualization plots
% ===========================================================================

% Define the analysis name for output organization
ANALYSIS_NAME = "DataInspection";

% Add TESA toolbox to path
addpath('C:\Program Files\MATLAB\R2024b\toolbox\eeglab2024.2\plugins\TESA1.1.1');

%% =======================================================================
%  PROCESSING PARAMETERS
%  =======================================================================

% Channel rejection parameters
activity_prob = 5;             % Threshold (std dev) for automatic bad channel detection

% Bad epoch detection parameters
tr_rej_stnddev_first = 5;      % Threshold (std dev) for joint probability bad epoch detection

% TMS pulse detection parameters
trigger_label = 'TMS';         % Event marker for TMS pulses in EEG data

% Epoching parameters
epoching_long = [-1 1];        % Epoch window around TMS pulse (seconds)
demeaning_interval = [-1000 999]; % Time window for epoch demeaning (ms)

%% =======================================================================
%  ENVIRONMENT SETUP AND PATH CONFIGURATION
%  =======================================================================

% Load environment variables from .env file
ENV_PATH = fullfile('..', '.env');
env = loadenv(ENV_PATH);

% Read experiment parameters from environment variables
ROOT_DIR        = getenv('TMS_EEG_ROOT_DIR');    % Base directory for all data
EXPERIMENT_NAME = getenv('EXPERIMENT_NAME');     % Name of current experiment
PARTICIPANT_ID  = getenv('PARTICIPANT_ID');      % Current participant identifier

% Validate that required environment variables are set
assert(~isempty(ROOT_DIR),       'TMS_EEG_ROOT_DIR is not set.');
assert(isfolder(ROOT_DIR),       'Folder "%s" does not exist.', ROOT_DIR);
assert(~isempty(EXPERIMENT_NAME), 'EXPERIMENT_NAME is not set.');
assert(~isempty(PARTICIPANT_ID), 'PARTICIPANT_ID is not set.');

% Construct data paths
data_root = fullfile(ROOT_DIR, EXPERIMENT_NAME, PARTICIPANT_ID, 'data');
assert(isfolder(data_root), 'Data directory "%s" does not exist.', data_root);

% Create output directory for this analysis
output_root = fullfile(ROOT_DIR, EXPERIMENT_NAME, PARTICIPANT_ID, 'output', ANALYSIS_NAME);
if ~exist(output_root, 'dir'); mkdir(output_root); end

%% =======================================================================
%  FIND ALL .VHDR FILES FOR PROCESSING
%  =======================================================================

% Search for all BrainVision header files in the data directory
vhdr_files = dir(fullfile(data_root, '*.vhdr'));
num_files = length(vhdr_files);

if num_files == 0
    error('No .vhdr files found in directory: %s', data_root);
end

fprintf('\n=== DATA INSPECTION PROCESSING ===\n');
fprintf('Experiment: %s\n', EXPERIMENT_NAME);
fprintf('Analysis: %s\n', ANALYSIS_NAME);
fprintf('Participant ID: %s\n', PARTICIPANT_ID);
fprintf('Data directory: %s\n', data_root);
fprintf('Output directory: %s\n', output_root);
fprintf('Found %d .vhdr files to process\n', num_files);
fprintf('=====================================\n\n');

%% =======================================================================
%  INITIALIZE EEGLAB
%  =======================================================================
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%% =======================================================================
%  PROCESS EACH .VHDR FILE
%  =======================================================================

for file_idx = 1:num_files
    
    % Get current file information
    current_file = vhdr_files(file_idx);
    [~, current_datasets_savename, ~] = fileparts(current_file.name);
    
    fprintf('Processing file %d/%d: %s\n', file_idx, num_files, current_file.name);
    
    % Create output subdirectory for this dataset
    current_output_folder = char(fullfile(output_root, current_datasets_savename));
    if ~exist(current_output_folder, 'dir'); mkdir(current_output_folder); end
    
    try
        %% ===============================================================
        %  STEP 1: LOAD EEG DATA
        %  ===============================================================
        
        fprintf('  Step 1: Loading EEG data...\n');
        
        % Load EEG data using BrainVision format reader
        EEG = pop_loadbv(data_root, current_file.name);
        EEG = eeg_checkset(EEG);
        EEG.setname = sprintf('%s_%s_%s_%s', ANALYSIS_NAME, EXPERIMENT_NAME, PARTICIPANT_ID, current_datasets_savename);
        
        % Remove auxiliary channels (e.g., EMG)
        if any(strcmp({EEG.chanlocs.labels}, 'Input 33'))
            EEG = pop_select(EEG, 'nochannel', {'Input 33'});
        end
        
        % Apply standard 10-5 electrode system naming convention
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
        
        % Validate channel count matches expected electrode configuration
        if numel(EEG.chanlocs) == numel(channel_names)
            % Assign standardized channel labels
            for k = 1:numel(channel_names)
                EEG.chanlocs(k).labels = channel_names{k};
            end
            
            % Update channel locations
            EEG = eeg_checkset(EEG);
            EEG = pop_chanedit(EEG, 'lookup','standard-10-5-cap385.elp');
        else
            warning('Channel count mismatch for %s: expected %d, got %d. Skipping channel renaming.', ...
                current_datasets_savename, numel(channel_names), numel(EEG.chanlocs));
        end
        
        fprintf('    Loaded: %d channels, %.1f Hz, %.1f seconds\n', ...
                EEG.nbchan, EEG.srate, EEG.pnts/EEG.srate);
        
        %% ===============================================================
        %  STEP 2: FIND TMS PULSES
        %  ===============================================================
        
        fprintf('  Step 2: Detecting TMS pulses...\n');
        
        % Detect TMS pulses using Cz electrode as reference
        EEG = pop_tesa_findpulse(EEG, 'Cz', ...
                                 'refract', 4, ...
                                 'rate', 1e4, ...
                                 'tmsLabel', trigger_label, ...
                                 'plots', 'off');  % Turn off plots for batch processing
        EEG = eeg_checkset(EEG);
        
        % Report number of detected pulses
        tms_events = sum(strcmp({EEG.event.type}, trigger_label));
        fprintf('    Detected %d TMS pulses\n', tms_events);
        
        %% ===============================================================
        %  STEP 3: REMOVE BAD CHANNELS USING KURTOSIS
        %  ===============================================================
        
        fprintf('  Step 3: Detecting and removing bad channels...\n');
        
        % Store original channel information for later interpolation
        EEG.allchan = EEG.chanlocs;
        
        % Remove channels based on kurtosis (measure of signal distribution)
        % High kurtosis indicates presence of artifacts or noise
        EEG = pop_rejchan(EEG, 'elec', 1:size(EEG.data,1), ...
                          'threshold', activity_prob, ...  % Threshold in std deviations
                          'norm', 'on', ...               % Normalize measures
                          'measure', 'kurt');             % Use kurtosis as rejection criterion
        
        % Report number of removed channels
        num_removed_channels = length(EEG.allchan) - EEG.nbchan;
        fprintf('    Removed %d bad channels\n', num_removed_channels);
        
        %% ===============================================================
        %  STEP 4: SAVE DATASET (FIRST TIME - AFTER CHANNEL REMOVAL)
        %  ===============================================================
        
        fprintf('  Step 4: Saving dataset after channel removal...\n');
        
        EEG.setname = [current_datasets_savename '_ChannelRemoved'];
        EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
        
        %% ===============================================================
        %  STEP 5: EPOCH DATA AROUND TMS PULSES
        %  ===============================================================
        
        fprintf('  Step 5: Epoching data around TMS pulses...\n');
        
        % Epoch data around TMS pulses
        EEG = pop_epoch(EEG, {trigger_label}, epoching_long, 'epochinfo', 'yes');
        
        % Apply demeaning correction
        EEG = pop_rmbase(EEG, demeaning_interval);
        
        fprintf('    Created %d epochs of %.1f seconds each\n', EEG.trials, epoching_long(2)-epoching_long(1));
        
        %% ===============================================================
        %  STEP 6: REMOVE BAD EPOCHS USING JOINT PROBABILITY
        %  ===============================================================
        
        fprintf('  Step 6: Detecting and removing bad epochs...\n');
        
        % Apply joint probability method to detect bad epochs
        EEG = pop_jointprob(EEG, 1, 1:size(EEG.data,1), tr_rej_stnddev_first, tr_rej_stnddev_first, 0, 0);
        
        % Count and report rejected epochs
        num_rejected_epochs = length(find(EEG.reject.rejjp));
        fprintf('    Number of rejected epochs: %d\n', num_rejected_epochs);
        
        % Store rejected epoch indices for record keeping
        bad_epochs = find(EEG.reject.rejjp);
        save(fullfile(current_output_folder, ['bad_epochs_' current_datasets_savename]), 'bad_epochs');
        
        % Remove marked epochs from dataset
        EEG = pop_rejepoch(EEG, EEG.reject.rejjp, 0);
        
        fprintf('    Remaining epochs after rejection: %d\n', EEG.trials);
        
        %% ===============================================================
        %  STEP 7: SAVE DATASET (SECOND TIME - AFTER EPOCHING AND CLEANING)
        %  ===============================================================
        
        fprintf('  Step 7: Saving final epoched and cleaned dataset...\n');
        
        EEG.setname = [current_datasets_savename '_EpochedCleaned'];
        EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
        
        %% ===============================================================
        %  STEP 8: GENERATE VISUALIZATION AND SAVE FIGURES
        %  ===============================================================
        
        fprintf('  Step 8: Generating visualization plots...\n');
        
        % Generate TMS-evoked potential plot
        visualize_eeg_evoked(EEG, ...
            'BaselineWindow', demeaning_interval, ...
            'TimeLimits', [-5, 30], ...
            'AmplitudeLimits', [-300, 300], ...
            'PlotTitle', sprintf('TMS-evoked potentials (All channels) - %s (n=%d epochs)', current_datasets_savename, EEG.trials), ...
            'SaveFig', true, ...
            'OutputPath', current_output_folder, ...
            'Filename', [current_datasets_savename '_evoked_potentials_AllChannels'], ...
            'TitleInterpreter', 'none');


        
        % Generate four-channel TMS-evoked potential plot
        visualize_eeg_evoked(EEG, ...
            'BaselineWindow', demeaning_interval, ...
            'TimeLimits', [-5, 30], ...
            'AmplitudeLimits', [-300, 300], ...
            'PlotTitle', sprintf('TMS-evoked potentials (Left Frontal channels) - %s (n=%d epochs)', current_datasets_savename, EEG.trials), ...
            'ChannelNames', {'F3', 'F1', 'FC3', 'FC1', 'C3', 'C1'}, ...
            'SaveFig', true, ...
            'OutputPath', current_output_folder, ...
            'Filename', [current_datasets_savename '_evoked_potentials_LeftFrontal'], ...
            'TitleInterpreter', 'none');

        % right frontal channels
        visualize_eeg_evoked(EEG, ...
            'BaselineWindow', demeaning_interval, ...
            'TimeLimits', [-5, 30], ...
            'AmplitudeLimits', [-300, 300], ...
            'PlotTitle', sprintf('TMS-evoked potentials (Right Frontal channels) - %s (n=%d epochs)', current_datasets_savename, EEG.trials), ...
            'ChannelNames', {'F4', 'F2', 'FC4', 'FC2', 'C4', 'C2'}, ...
            'SaveFig', true, ...
            'OutputPath', current_output_folder, ...
            'Filename', [current_datasets_savename '_evoked_potentials_RightFrontal'], ...
            'TitleInterpreter', 'none');

        % center frontal channels
        visualize_eeg_evoked(EEG, ...
            'BaselineWindow', demeaning_interval, ...
            'TimeLimits', [-5, 30], ...
            'AmplitudeLimits', [-300, 300], ...
            'PlotTitle', sprintf('TMS-evoked potentials (Center Frontal channels) - %s (n=%d epochs)', current_datasets_savename, EEG.trials), ...
            'ChannelNames', {'FZ', 'FCz', 'Cz'}, ...
            'SaveFig', true, ...
            'OutputPath', current_output_folder, ...
            'Filename', [current_datasets_savename '_evoked_potentials_CenterFrontal'], ...
            'TitleInterpreter', 'none');

        % Generate individual channel plots for all remaining channels
        fprintf('  Step 8b: Generating individual channel plots...\n');
        for ch_idx = 1:EEG.nbchan
            current_channel = EEG.chanlocs(ch_idx).labels;
            fprintf('    Plotting channel %d/%d: %s\n', ch_idx, EEG.nbchan, current_channel);
            
            visualize_eeg_evoked(EEG, ...
                'BaselineWindow', demeaning_interval, ...
                'TimeLimits', [-5, 30], ...
                'AmplitudeLimits', [-300, 300], ...
                'PlotTitle', sprintf('TMS-evoked potential - Channel %s - %s (n=%d epochs)', current_channel, current_datasets_savename, EEG.trials), ...
                'ChannelNames', {current_channel}, ...
                'SaveFig', true, ...
                'OutputPath', current_output_folder, ...
                'Filename', [current_datasets_savename '_evoked_potential_' current_channel], ...
                'TitleInterpreter', 'none');
        end
        fprintf('    Completed individual channel plots for %d channels\n', EEG.nbchan);
        
        % Close all figure windows to free up memory and avoid clutter
        close all;
        
        fprintf('  ✓ Successfully processed %s\n\n', current_datasets_savename);
        
    catch ME
        % Handle errors for individual files without stopping the entire batch
        fprintf('  ✗ ERROR processing %s: %s\n\n', current_datasets_savename, ME.message);
        continue;
    end
end

%% =======================================================================
%  PROCESSING SUMMARY
%  =======================================================================

fprintf('\n=== PROCESSING COMPLETE ===\n');
fprintf('Processed %d files\n', num_files);
fprintf('Results saved in: %s\n', output_root);
fprintf('============================\n');



