clear; close all; clc;

% ===========================================================================
% TESA (TMS-EEG Signal Analyzer) Complete Processing Pipeline
% ===========================================================================
% This pipeline implements a comprehensive TMS-EEG preprocessing workflow
% following TESA guidelines for artifact removal and signal cleaning.
% 
% Key processing steps:
% 1. Data loading and channel setup
% 2. TMS pulse detection and marking
% 3. Bad channel detection and removal
% 4. Epoching around TMS pulses
% 5. TMS pulse artifact removal and interpolation
% 6. Downsampling and bad trial rejection
% 7. Two-stage ICA decomposition for artifact removal
% 8. Channel interpolation and re-referencing
% 9. Final baseline correction and visualization
%
% Pipeline optimized for single-pulse TMS protocols with concurrent EEG.
% ===========================================================================

% Define the pipeline name for tracking/logging purposes
PIPELINE_NAME = 'TESA_ICA';   

%% =======================================================================
%  PREPROCESSING PARAMETERS
%  =======================================================================
% These parameters control the entire preprocessing pipeline and should be
% adjusted based on your experimental protocol and data characteristics.

% Channel rejection parameters
activity_prob = 5;              % Threshold (std dev) for automatic bad channel detection
                               % Higher values = more conservative (fewer channels rejected)

% Temporal parameters  
epoching_long = [-1 1];        % Epoch window around TMS pulse (seconds)
                               % [-1 1] = 1 second before to 1 second after TMS
demeaning_interval = [-1000 999]; % Time window for epoch demeaning (ms)
                               % Removes DC drift by subtracting mean voltage
baseline_long = [-1000 -2];   % Baseline correction window (ms)
                               % Should end before TMS pulse to avoid artifacts

% TMS artifact handling
trigger_label = 'TMS';         % Event marker for TMS pulses in EEG data
interp_interval = [-1 6];     % Time window around TMS pulse to interpolate (ms)
                               % [-1 6] captures direct TMS artifact and recharge

% Signal processing parameters
downsample = 1000;             % Target sampling rate (Hz) after TMS artifact removal
                               % Reduces computational load while preserving TEP features
tr_rej_stnddev_first = 5;      % Threshold (std dev) for automatic trial rejection
                               % Based on joint probability of channel activity
low_pass_filt = 90;           % Low-pass filter cutoff (Hz) - removes high-freq noise
hp_filt = 1;                  % High-pass filter cutoff (Hz) - removes slow drifts
notch_filt = [58 62];         % Notch filter range (Hz) for line noise removal
                               % [58 62] targets 60 Hz power line interference
final_ref = [];               % Final reference: [] = average reference
                               % Alternative: specific channels like {'Cz'}

%% =======================================================================
%  ENVIRONMENT SETUP AND DATA LOADING
%  =======================================================================

%% 0. Load environment variables from .env file ----------------------------
%    Cross-platform configuration management for data paths
ENV_PATH = fullfile('..', '.env');   % Path to environment configuration file
env = loadenv(ENV_PATH);             % Load variables if file exists

%% 1. Read experiment parameters from environment variables ----------------
%    Essential paths for locating raw TMS-EEG data files
ROOT_DIR        = getenv('TMS_EEG_ROOT_DIR');    % Base directory for all data
EXPERIMENT_NAME = getenv('EXPERIMENT_NAME');     % Name of current experiment
PARTICIPANT_ID  = getenv('PARTICIPANT_ID');      % Current participant identifier

% Validate that required environment variables are set
assert(~isempty(ROOT_DIR),       'TMS_EEG_ROOT_DIR is not set.');
assert(isfolder(ROOT_DIR),       'Folder "%s" does not exist.', ROOT_DIR);

%% 2. Construct data paths -------------------------------------------------
%    Build standardized file paths for input data and output results
data_root = fullfile(ROOT_DIR, EXPERIMENT_NAME, PARTICIPANT_ID, 'data');

% -------------------------------------------------------------------------
%    Dataset identifier - change this to match your file naming convention
current_datasets_savename = 'Pos10_80';   % e.g., stimulation condition/session
% -------------------------------------------------------------------------

%    Locate BrainVision header file (.vhdr) containing EEG data
vhdr_file = fullfile(data_root, [current_datasets_savename '.vhdr']);
assert(isfile(vhdr_file), ...
       'Session file "%s" does not exist.', vhdr_file);

% Create output directory structure for processed data and intermediate files
current_output_folder = char( fullfile( ...
        ROOT_DIR, EXPERIMENT_NAME, PARTICIPANT_ID, ...
        'output', PIPELINE_NAME, current_datasets_savename) );
if ~exist(current_output_folder, 'dir'); mkdir(current_output_folder); end

%% =======================================================================
%  STAGE 1: DATA IMPORT AND INITIAL SETUP
%  =======================================================================

%% 3. Load EEG data into EEGLAB -------------------------------------------
%    Import raw TMS-EEG data using BrainVision format reader
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
EEG = pop_loadbv(data_root, [current_datasets_savename '.vhdr']);
EEG = eeg_checkset(EEG);
EEG.setname = sprintf('%s_%s_%s', ...
        EXPERIMENT_NAME, PARTICIPANT_ID, current_datasets_savename);

%% 4. Preprocessing: Clean up channel list and apply standard names --------
%    Remove non-EEG channels and standardize electrode nomenclature

% Remove auxiliary channels that don't contain EEG data
% 'Input 33' is an EMG channel that is not needed for TMS-EEG analysis
EEG = pop_select(EEG, 'nochannel', {'Input 33'});

%    Apply standard 10-5 electrode system naming convention
%    This ensures compatibility with EEGLAB functions and facilitates
%    cross-study comparisons. Order should match your cap layout.
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

%    Validate channel count matches expected electrode configuration
assert(numel(EEG.chanlocs) == numel(channel_names), ...
       'Dataset has %d channels but you supplied %d names.', ...
       numel(EEG.chanlocs), numel(channel_names));

%    Assign standardized channel labels to electrode locations
for k = 1:numel(channel_names)
    EEG.chanlocs(k).labels = channel_names{k};
end

%    Update channel locations and perform internal consistency checks
EEG = eeg_checkset(EEG);
% Load standard 10-5 electrode positions for topographic plotting
EEG = pop_chanedit(EEG, 'lookup','standard-10-5-cap385.elp');

%% 5. Display basic EEG dataset information --------------------------------
%    Provide overview of loaded dataset characteristics
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


%% Save initial dataset after channel setup
pop_saveset(...
    EEG,...
    'filename', [current_datasets_savename '_SetCreation.set'],...
    'filepath', current_output_folder);

%% =======================================================================
%  STAGE 2: TMS PULSE DETECTION AND BAD CHANNEL REMOVAL
%  =======================================================================

%% 6. TMS pulse detection and event marking --------------------------------
%    Automatically detect TMS pulses in the EEG data using the Cz electrode
%    This creates event markers that will be used for epoching
EEG = pop_tesa_findpulse(EEG,'Cz', ...         % Reference electrode for pulse detection
                         'refract',4, ...       % Refractory period (ms) between pulses
                         'rate',1e4, ...        % Pulse detection threshold rate
                         'tmsLabel', trigger_label, ... % Event label for detected pulses
                         'plots','on');         % Display detection results
EEG = eeg_checkset(EEG);



%% 7. Automatic bad channel detection and removal ---------------------------
%    Identify and remove channels with excessive noise or artifacts
%    This step is crucial for preventing bad channels from contaminating ICA

% Store original channel information for later interpolation
EEG.allchan = EEG.chanlocs;

% Remove channels based on kurtosis (measure of signal distribution)
% High kurtosis indicates presence of artifacts or noise
EEG = pop_rejchan(EEG, 'elec', 1:size(EEG.data,1) ,...
                  'threshold',activity_prob,...  % Threshold in std deviations
                  'norm','on',...                % Normalize measures
                  'measure','kurt');             % Use kurtosis as rejection criterion

% Save dataset after bad channel removal
EEG.setname = [current_datasets_savename '_ChanRem'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);


%% =======================================================================
%  STAGE 3: EPOCHING AND INITIAL ARTIFACT REMOVAL
%  =======================================================================

%% 8. Epoch data around TMS pulses ------------------------------------------
%    Create time-locked epochs centered on each TMS pulse
%    This segments continuous data into trials for further processing
EEG = pop_epoch( EEG, { trigger_label }, epoching_long , 'epochinfo', 'yes');

% Save epoched dataset
EEG.setname = [current_datasets_savename '_epochs'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

% plot and save the figures
visualize_eeg_evoked( ...
        EEG, ...                 % Epoched dataset
        demeaning_interval, ...  % Use demeaning interval as baseline window
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-50 50], ...            % Amplitude scale: ±50 µV (wider range for less processed data)
        'After Epoching - TMS-evoked potentials', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_after_epoching']); % filename base

%% 9. Demeaning correction (DC removal) -------------------------------------
%    Remove the mean voltage across each epoch to eliminate DC drifts
%    This step improves ICA decomposition by centering the data around zero
EEG = pop_rmbase( EEG, demeaning_interval);

% Save demeaned dataset
EEG.setname = [current_datasets_savename '_Demeaning'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

% plot and save the figures
visualize_eeg_evoked( ...
        EEG, ...                 % Demeaned dataset
        [], ...                  % No additional baseline correction (already demeaned)
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-50 50], ...            % Amplitude scale: ±50 µV
        'After Demeaning - TMS-evoked potentials', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_after_demeaning']); % filename base

%% 10. Remove TMS pulse artifact (first pass) ------------------------------
%     Replace data around TMS pulse with flat line to prevent contamination
%     This temporary removal prevents TMS artifacts from affecting downsampling
EEG = pop_tesa_removedata( EEG, interp_interval);

% Save dataset with TMS pulse removed
EEG.setname = [current_datasets_savename '_TMSpulseREM'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

% plot and save the figures
% Note: This will show the flat line where TMS pulse was removed
visualize_eeg_evoked( ...
        EEG, ...                 % Dataset with TMS pulse removed
        [], ...                  % No baseline correction (data already demeaned)
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-50 50], ...            % Amplitude scale: ±50 µV
        'After TMS Pulse Removal - showing flat interpolation window', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_after_tms_pulse_removal']); % filename base

%% 11. Interpolate missing data around TMS pulse (first pass) ---------------
%     Use cubic spline interpolation to fill the gap left by pulse removal
%     [1,1] means use 1 ms before and after the removed window for interpolation
EEG = pop_tesa_interpdata( EEG, 'cubic', [1,1] );

% Save dataset with interpolated TMS pulse
EEG.setname = [current_datasets_savename '_TMSpInt'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

% plot and save the figures
visualize_eeg_evoked( ...
        EEG, ...                 % Dataset with interpolated TMS pulse
        [], ...                  % No baseline correction (data already demeaned)
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-50 50], ...            % Amplitude scale: ±50 µV
        'After TMS Pulse Interpolation', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_after_tms_pulse_interpolation']); % filename base


%% =======================================================================
%  STAGE 4: DOWNSAMPLING AND TRIAL REJECTION
%  =======================================================================

%% 12. Downsample data to reduce computational load -------------------------
%     Reduce sampling rate after TMS artifact removal to prevent aliasing
%     This speeds up subsequent ICA processing while preserving TEP features
EEG = pop_resample( EEG, downsample);

EEG.setname = [current_datasets_savename '_DownSamp'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% 13. Automatic bad trial detection and removal ----------------------------
%     Remove epochs with excessive artifacts using joint probability method
%     This identifies trials where multiple channels show unusual activity

% Calculate joint probability across channels and mark bad epochs
EEG = pop_jointprob(EEG,1, 1:size(EEG.data,1),tr_rej_stnddev_first,tr_rej_stnddev_first,0,0);
disp(['number of rejected epochs: ' num2str(length(find(EEG.reject.rejjp)))]);

% Store rejected epoch indices for record keeping
bad_epochs=find(EEG.reject.rejjp);
save([current_output_folder '\bad_epochs_' current_datasets_savename], 'bad_epochs');

% Remove marked epochs from dataset
EEG = pop_rejepoch( EEG, EEG.reject.rejjp, 0);

%% =======================================================================
%  STAGE 5: FIRST ICA DECOMPOSITION (TMS-MUSCLE ARTIFACTS)
%  =======================================================================

%% 14. Prepare data for ICA by removing interpolated segments ---------------
%     Replace interpolated data with flat values because interpolated data
%     can bias ICA decomposition (adds artificial correlations)
EEG = pop_tesa_removedata( EEG, interp_interval );

EEG.setname = [current_datasets_savename '_TMSpulse0'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

% Plot and save figures before first ICA
visualize_eeg_evoked( ...
        EEG, ...                 % Data prepared for first ICA
        [], ...                  % No baseline correction (data already demeaned)
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-50 50], ...            % Amplitude scale: ±50 µV
        'Before First ICA - Data prepared for component analysis', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_before_first_ICA']); % filename base

%% 15. First ICA decomposition for TMS-evoked muscle artifacts --------------
%     Decompose EEG data into independent components to isolate artifacts
%     FastICA with symmetric approach and tanh nonlinearity is optimal for TMS-EEG
EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );

EEG.setname = [current_datasets_savename '_ICA-TMSmuscle'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% 16. First component selection (focus on TMS-muscle artifacts) ------------
%     Automatically identify and remove components containing TMS-evoked muscle activity
%     This first pass focuses specifically on muscle artifacts in the 11-30 ms window

EEG = pop_tesa_compselect( EEG, ...
    'compCheck','off',...              % Disable manual component checking
    'comps', 15, ...                   % Number of components to evaluate
    'figSize','small',...              % Component plot size
    'plotTimeX',[-100 399],...         % Time window for component visualization
    'plotFreqX',[1 100],...            % Frequency range for spectral analysis
    'tmsMuscle','on',...               % Enable TMS-muscle artifact detection
    'tmsMuscleThresh',8,...            % Threshold for muscle artifact identification
    'tmsMuscleWin',[11 30],...         % Time window where muscle artifacts occur
    'tmsMuscleFeedback','off',...      % Disable feedback plots
    'blink','off',...                  % Disable blink detection (first pass)
    'blinkThresh',2.5,...              % Blink detection threshold
    'blinkElecs',{'Fp1','Fp2'},...     % Electrodes for blink detection
    'blinkFeedback','off',...          % Disable blink feedback plots
    'move','off',...                   % Disable movement detection (first pass)
    'moveThresh',2,...                 % Movement detection threshold
    'moveElecs',{'F7','F8'},...        % Electrodes for movement detection
    'moveFeedback','off',...           % Disable movement feedback plots
    'muscle','off',...                 % Disable general muscle detection (first pass)
    'muscleThresh',0.6,...             % General muscle detection threshold
    'muscleFreqWin',[30 100],...       % Frequency window for muscle activity
    'muscleFeedback','off',...         % Disable muscle feedback plots
    'elecNoise','off',...              % Disable electrode noise detection (first pass)
    'elecNoiseThresh',4,...            % Electrode noise threshold
    'elecNoiseFeedback','off' );       % Disable electrode noise feedback

EEG.setname = [current_datasets_savename '_CompSel1'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

% Plot and save figures after first component selection
visualize_eeg_evoked( ...
        EEG, ...                 % Data after first ICA component removal
        [], ...                  % No baseline correction (data already demeaned)
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-50 50], ...            % Amplitude scale: ±50 µV
        'After First ICA - TMS-muscle artifacts removed', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_after_first_ICA']); % filename base

%% =======================================================================
%  STAGE 6: FILTERING AND SECOND ICA DECOMPOSITION
%  =======================================================================

%% 17. Load dataset after first component selection -------------------------
EEG = pop_loadset(...
    'filename', [current_datasets_savename  '_CompSel1.set'],...
    'filepath',  current_output_folder);

%% 18. Remove data around TMS pulse again for filtering ---------------------
%     This prevents filtering artifacts at the edges of the removed window
EEG = pop_tesa_removedata( EEG, interp_interval );

EEG.setname = [current_datasets_savename '_TMSp0ext'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% 19. Interpolate data before filtering ------------------------------------
%     Replace flat values with smooth interpolation to prevent filtering artifacts
%     [5,5] uses 5 ms before and after the removed window for interpolation
EEG = pop_tesa_interpdata( EEG, 'cubic', [5,5] );

EEG.setname = [current_datasets_savename '_TMSp0Int'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% 20. Apply frequency domain filtering -------------------------------------
%     Remove frequency components outside the range of interest
%     Butterworth filters provide good frequency response with minimal ringing

% Bandpass filter: remove slow drifts and high-frequency noise
EEG = pop_tesa_filtbutter( EEG, hp_filt, low_pass_filt, 4, 'bandpass' );
% Notch filter: remove power line interference (50/60 Hz)
EEG = pop_tesa_filtbutter( EEG, notch_filt(1), notch_filt(2), 4, 'bandstop' );

EEG.setname = [current_datasets_savename '_notch&90_1_BandpassFilter'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% 21. Remove interpolated data again for second ICA ------------------------
%     Prepare clean data for final ICA decomposition
EEG = pop_tesa_removedata( EEG, interp_interval );

EEG.setname = [current_datasets_savename '_TMSp0'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

% Plot and save figures before second ICA
visualize_eeg_evoked( ...
        EEG, ...                 % Filtered data prepared for second ICA
        [], ...                  % No baseline correction (data already demeaned)
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-30 30], ...            % Amplitude scale: ±30 µV (smaller range after filtering)
        'Before Second ICA - After filtering, ready for final component analysis', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_before_second_ICA']); % filename base

%% 22. Second ICA decomposition for remaining artifacts ---------------------
%     Final ICA run to remove any remaining artifacts (blinks, movements, etc.)
EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );

EEG.setname = [current_datasets_savename '_ICA2'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% 23. Second component selection (comprehensive artifact removal) ----------
%     Remove all types of artifacts in this final component selection pass
EEG = pop_tesa_compselect( EEG,...
    'compCheck','off',...              % Disable manual checking
    'comps',[],...                     % Evaluate all components
    'figSize','medium',...             % Larger plots for detailed inspection
    'plotTimeX',[-100 399],...         % Time window for visualization
    'plotFreqX',[1 100],...            % Frequency range for analysis
    'tmsMuscle','on',...               % Continue TMS-muscle detection
    'tmsMuscleThresh',8,...            % Same threshold as before
    'tmsMuscleWin',[11 30],...         % TMS-muscle time window
    'tmsMuscleFeedback','off',...      % Disable feedback
    'blink','off',...                   % Enable blink artifact detection
    'blinkThresh',2.5,...              % Blink detection sensitivity
    'blinkElecs',{'Fp1','Fp2'},...     % Frontal electrodes for blinks
    'blinkFeedback','off',...          % Disable blink feedback
    'move','on',...                    % Enable movement artifact detection
    'moveThresh',2,...                 % Movement detection sensitivity
    'moveElecs',{'F7','F8'},...        % Temporal electrodes for movements
    'moveFeedback','off',...           % Disable movement feedback
    'muscle','on',...                  % Enable general muscle detection
    'muscleThresh',0.6,...             % Muscle detection threshold
    'muscleFreqWin',[30 100],...       % High-frequency muscle activity
    'muscleFeedback','off',...         % Disable muscle feedback
    'elecNoise','on',...               % Enable electrode noise detection
    'elecNoiseThresh',4,...            % Electrode noise threshold
    'elecNoiseFeedback','off' );       % Disable noise feedback

EEG.setname = [current_datasets_savename '_CompSel2'];

EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

% Plot and save figures after second component selection
visualize_eeg_evoked( ...
        EEG, ...                 % Data after comprehensive artifact removal
        [], ...                  % No baseline correction (data already demeaned)
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-30 30], ...            % Amplitude scale: ±30 µV
        'After Second ICA - All artifacts removed (muscle, blink, movement, noise)', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_after_second_ICA']); % filename base

%% =======================================================================
%  STAGE 7: FINAL PROCESSING AND RECONSTRUCTION
%  =======================================================================

%% 24. Load dataset after second component selection ------------------------
EEG = pop_loadset(...
    'filename', [current_datasets_savename  '_CompSel2.set'],...
    'filepath', current_output_folder);

%% 25. Final interpolation of TMS pulse window ------------------------------
%     Restore smooth data in the TMS pulse region for final analysis
EEG = pop_tesa_interpdata( EEG, 'cubic', [5,5] );

EEG.setname = [current_datasets_savename '_After2ndCompSel_TMSpInt'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% 26. Interpolate bad channels back into the dataset ----------------------
%     Restore channels that were removed earlier using spherical interpolation
%     This gives a complete electrode montage for final analysis
EEG = pop_interp(EEG, EEG.allchan, 'spherical');

EEG.setname = [current_datasets_savename '_ChanInt'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% 27. Re-reference to average of all electrodes ----------------------------
%     Apply average reference to minimize the influence of reference electrode choice
%     Data were originally recorded with a single reference (e.g., linked mastoids)
EEG = pop_reref( EEG, final_ref);

EEG.setname = [current_datasets_savename '_AvgReref'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

% Plot and save figures after re-referencing
visualize_eeg_evoked( ...
        EEG, ...                 % Data after average re-referencing
        [], ...                  % No baseline correction (will be applied next)
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-30 30], ...            % Amplitude scale: ±25 µV (tighter range after re-referencing)
        'After Average Re-referencing - Ready for final baseline correction', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_after_average_reref']); % filename base

%% 28. Final baseline correction --------------------------------------------
%     Remove pre-stimulus baseline to isolate TMS-evoked responses
%     Baseline window should end before TMS pulse to avoid contamination
EEG = pop_rmbase( EEG, baseline_long);

EEG.setname = [current_datasets_savename '_baselineCorr'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% =======================================================================
%  STAGE 8: FINAL DATASET AND VISUALIZATION
%  =======================================================================

%% 29. Save final processed dataset -----------------------------------------
EEG.setname= [current_datasets_savename 'TESA_Processed'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% 30. Generate final TMS-evoked potential visualization -------------------
%     Create publication-quality plot of the cleaned TMS-evoked potentials
%     Display average response across all channels with appropriate scaling

% Generate TMS-evoked potential plot with baseline correction and automatic saving
visualize_eeg_evoked( ...
        EEG, ...                 % Final processed dataset
        baseline_long, ...       % Baseline window (e.g. [-1000 -2] ms)
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-30 30], ...            % Amplitude scale: ±25 µV (adjust as needed)
        'Final TESA-processed TMS-evoked potentials', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_Final_scaled_avg_ref']); % filename base

%% 31. Pipeline completion and cleanup --------------------------------------
diary off;           % Stop logging if active

% ===========================================================================
% Pipeline completion message
% ===========================================================================
fprintf('\n=== TESA PROCESSING COMPLETE ===\n');
fprintf('Processed dataset: %s\n', EEG.setname);
fprintf('Output directory : %s\n', current_output_folder);
fprintf('Final data shape : [%d channels × %d timepoints × %d epochs]\n', ...
        size(EEG.data,1), size(EEG.data,2), size(EEG.data,3));
fprintf('Sampling rate    : %.1f Hz\n', EEG.srate);
fprintf('====================================\n');


