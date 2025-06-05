clear; close all; clc;

% ===========================================================================
% SOUND + SSP-SIR TMS-EEG Processing Pipeline (Default TESA Implementation)
% ===========================================================================
% This pipeline implements an advanced TMS-EEG preprocessing workflow
% using the default TESA implementations of SOUND (Suppression of Unwanted 
% Noise and Disturbances) and SSP-SIR (Signal Space Separation - Source-Informed 
% Reconstruction) algorithms for state-of-the-art artifact removal and signal cleaning.
% 
% Key processing steps:
% 1. Data loading and basic preprocessing
% 2. SOUND algorithm for channel-specific noise removal (TESA default)
% 3. Manual artifact rejection with visual inspection
% 4. ICA decomposition for ocular artifact removal
% 5. SSP-SIR for muscle artifact removal (TESA default)
% 6. Final processing and signal reconstruction
%
% Pipeline optimized for research-grade TMS-EEG data analysis with maximum
% flexibility and comprehensive quality control using default TESA functions.
% ===========================================================================

% Define the pipeline name for tracking/logging purposes
PIPELINE_NAME = 'SOUND-SSPSIR-DEFAULT';   

%% =======================================================================
%  PREPROCESSING PARAMETERS
%  =======================================================================
% These parameters control the entire preprocessing pipeline and should be
% adjusted based on your experimental protocol and data characteristics.

% Temporal parameters
epoching_long = [-1 1];           % Epoch window around TMS pulse (seconds)
                                  % [-1 1] = 1 second before to 1 second after TMS
baseline_long = [-1000 -2];      % Baseline correction window (ms)
                                  % Should end before TMS pulse to avoid artifacts
demeaning_interval = [-1000 999]; % Time window for epoch demeaning (ms)

% TMS artifact handling
trigger_label = 'TMS';  % Define trigger label for TMS pulses
interp_interval = [-1 6];        % Time window around TMS pulse to interpolate (ms)
                                  % [-1 6] captures direct TMS artifact and recharge

% Trial rejection parameters
tr_rej_stnddev_first = 5;        % Standard deviation threshold for automatic trial rejection
                                 % Epochs exceeding this threshold will be marked for removal

% SOUND algorithm parameters (TESA default)
sound_lambdaValue = 0.1;          % Regularization parameter for minimum-norm estimation
                                  % 0.1 was used in the original SOUND paper
sound_iter = 5;                   % Number of iterations for noise estimation
                                  % 5 iterations found sufficient in original studies

% SSP-SIR parameters (TESA default)
sspsir_timerange = [-2 50];       % Time window for SSP-SIR muscle artifact estimation (ms)
                                  % Should start before interpolation window
max_SSP_SIR_PC_to_remove = 5;     % Maximum number of principal components to remove

% Signal processing parameters
downsample = 1000;                % Target sampling rate (Hz) after preprocessing
low_pass_filt = 90;              % Low-pass filter cutoff (Hz)
hp_filt = 1;                     % High-pass filter cutoff (Hz) - removes slow drifts
notch_filt = [58 62];            % Notch filter range (Hz) for line noise removal
final_ref = [];                  % Final reference: [] = average reference

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

% Start logging workspace for reproducibility
diary([current_output_folder '\' current_datasets_savename '_workspace_SOUND_SSP-SIR_DEFAULT']);

%% =======================================================================
%  STAGE 1: DATA IMPORT AND BASIC PREPROCESSING
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
% 
%% Save initial dataset after channel setup
pop_saveset(...
    EEG,...
    'filename', [current_datasets_savename '_SetCreation.set'],...
    'filepath', current_output_folder);

%% =======================================================================
%  STAGE 2: TMS PULSE DETECTION AND INITIAL PROCESSING
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

%% 7. High-pass filtering to remove slow drifts ---------------------------
EEG = pop_eegfiltnew(EEG, 'locutoff', hp_filt, 'plotfreqz', 1);
close; % Close filter response plot

%% 8. Epoching around TMS pulses ------------------------------------------
EEG = pop_epoch(EEG, {trigger_label}, epoching_long, 'epochinfo', 'yes');

%% 9. TMS pulse artifact removal and interpolation ------------------------
EEG_before_interpolation = EEG; % Store dataset before interpolation
EEG = pop_tesa_removedata(EEG, interp_interval); % Remove TMS pulse
EEG = pop_tesa_interpdata(EEG, 'cubic', [1,1]); % Interpolate with cubic spline

%% 10. Baseline correction ------------------------------------------------
EEG = pop_rmbase(EEG, baseline_long);

%% 11. Visualization after basic preprocessing ----------------------------
% plot and save the figures
visualize_eeg_evoked( ...
        EEG, ...                 % Basic preprocessed dataset
        [-750 0], ...            % Baseline window for visualization
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-50 50], ...            % Amplitude scale: ±50 µV (wider range for less processed data)
        'After Basic Pre-processing', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_basic_prep']); % filename base

% Save dataset
EEG.setname = [current_datasets_savename '_basic_prep'];
pop_saveset(EEG, 'filename', [current_datasets_savename '_basic_prep'], ...
    'filepath', current_output_folder);

fprintf('Basic preprocessing completed.\n');

%% =======================================================================
%  STAGE 3: SOUND ALGORITHM FOR CHANNEL-SPECIFIC NOISE REMOVAL (TESA DEFAULT)
%  =======================================================================

fprintf('\n--- Starting SOUND Algorithm (TESA Default) ---\n');

% Load dataset from previous stage
EEG = pop_loadset('filename', [current_datasets_savename '_basic_prep.set'], ...
    'filepath', current_output_folder);

% Apply SOUND algorithm using TESA default implementation
fprintf('Running TESA SOUND algorithm...\n');
EEG = pop_tesa_sound(EEG, ...
    'lambdaValue', sound_lambdaValue, ...    % Regularization parameter
    'iter', sound_iter, ...                  % Number of iterations
    'leadfieldInFile', [], ...               % Use default lead field
    'leadfieldChansFile', [], ...            % Use default channel selection
    'replaceChans', [], ...                  % Replace all channels
    'multipleConds', []);                    % Single condition

%% 12. Visualization after SOUND ------------------------------------------
% plot and save the figures
visualize_eeg_evoked( ...
        EEG, ...                 % Dataset after SOUND algorithm
        [-750 0], ...            % Baseline window for visualization
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-50 50], ...            % Amplitude scale: ±50 µV
        'After SOUND Algorithm (TESA Default)', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_afterSOUND']); % filename base

% Save dataset after SOUND
EEG.setname = [current_datasets_savename '_afterSOUND'];
pop_saveset(EEG, 'filename', [current_datasets_savename '_afterSOUND'], ...
    'filepath', current_output_folder);

fprintf('SOUND algorithm completed.\n');

%% =======================================================================
%  STAGE 4: AUTOMATIC ARTIFACT REJECTION
%  =======================================================================

fprintf('\n--- Automatic Artifact Rejection ---\n');

% Load dataset from previous stage
EEG = pop_loadset('filename', [current_datasets_savename '_afterSOUND.set'], ...
    'filepath', current_output_folder);

%% 13. Automatic bad trial detection and removal ----------------------------
%     Remove epochs with excessive artifacts using joint probability method
%     This identifies trials where multiple channels show unusual activity

fprintf('Running automatic trial rejection using joint probability method...\n');
fprintf('Standard deviation threshold: %.1f\n', tr_rej_stnddev_first);

% Calculate joint probability across channels and mark bad epochs
EEG = pop_jointprob(EEG, 1, 1:size(EEG.data,1), tr_rej_stnddev_first, tr_rej_stnddev_first, 0, 0);

% Display and save rejected epochs information
num_rejected = length(find(EEG.reject.rejjp));
fprintf('Number of rejected epochs: %d\n', num_rejected);

% Store rejected epoch indices for record keeping
bad_epochs = find(EEG.reject.rejjp);
save([current_output_folder '\bad_epochs_' current_datasets_savename], 'bad_epochs');

% Remove marked epochs from dataset
EEG = pop_rejepoch(EEG, EEG.reject.rejjp, 0);

fprintf('Automatic trial rejection completed. %d epochs removed.\n', num_rejected);

%% 14. Visualization after automatic trial rejection -----------------------
% plot and save the figures
visualize_eeg_evoked( ...
        EEG, ...                 % Dataset after automatic trial rejection
        [-750 0], ...            % Baseline window for visualization
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-50 50], ...            % Amplitude scale: ±50 µV
        'After Automatic Trial Rejection', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_trialrej']); % filename base

% Save dataset after automatic trial rejection
EEG.setname = [current_datasets_savename '_trialrej'];
pop_saveset(EEG, 'filename', [current_datasets_savename '_trialrej'], ...
    'filepath', current_output_folder);

fprintf('Automatic artifact rejection completed.\n');

%% =======================================================================
%  STAGE 5: ICA DECOMPOSITION FOR OCULAR ARTIFACTS (TESA METHOD)
%  =======================================================================

fprintf('\n--- ICA Decomposition for Ocular Artifacts (TESA Method) ---\n');
% 
% % Load dataset from previous stage
EEG = pop_loadset('filename', [current_datasets_savename '_trialrej.set'], ...
    'filepath', current_output_folder);

%% 15. Prepare data for ICA by removing interpolated segments ---------------
%     Replace interpolated data with flat values because interpolated data
%     can bias ICA decomposition (adds artificial correlations)
EEG = pop_tesa_removedata( EEG, interp_interval );

% Save dataset prepared for ICA
EEG.setname = [current_datasets_savename '_TMSpulse0_beforeICA'];
pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

% Plot and save figures before ICA
visualize_eeg_evoked( ...
        EEG, ...                 % Data prepared for ICA
        [], ...                  % No baseline correction (data already processed)
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-50 50], ...            % Amplitude scale: ±50 µV
        'Before ICA - Data prepared for component analysis', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_before_ICA']); % filename base

%% 16. ICA decomposition for ocular artifacts ------------------------------
%     Decompose EEG data into independent components to isolate artifacts
%     FastICA with symmetric approach and tanh nonlinearity is optimal for TMS-EEG
fprintf('Running ICA decomposition...\n');
EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );

% Save dataset after ICA
EEG.setname = [current_datasets_savename '_After_ICA'];
pop_saveset(EEG, 'filename', [current_datasets_savename '_After_ICA'], ...
    'filepath', current_output_folder);

fprintf('ICA decomposition completed.\n');

%% 17. Automatic component selection for ocular artifacts ------------------
%     Automatically identify and remove components containing ocular artifacts
%     Focus specifically on blinks and eye movements, not muscle artifacts

fprintf('Starting automatic component selection for ocular artifacts...\n');

EEG = pop_tesa_compselect( EEG, ...
    'compCheck','off',...              % Disable manual component checking
    'comps', 15, ...                   % Number of components to evaluate
    'figSize','small',...              % Component plot size
    'plotTimeX',[-100 399],...         % Time window for component visualization
    'plotFreqX',[1 100],...            % Frequency range for spectral analysis
    'tmsMuscle','off',...               % Enable TMS-muscle artifact detection
    'tmsMuscleThresh',8,...            % Threshold for muscle artifact identification
    'tmsMuscleWin',[11 30],...         % Time window where muscle artifacts occur
    'tmsMuscleFeedback','on',...       % Enable feedback plot
    'blink','on',...                  % Disable blink detection (first pass)
    'blinkThresh',2.5,...              % Blink detection threshold
    'blinkElecs',{'Fp1','Fp2'},...     % Electrodes for blink detection
    'blinkFeedback','on',...           % Enable blink feedback plots
    'move','on',...                   % Disable movement detection (first pass)
    'moveThresh',2,...                 % Movement detection threshold
    'moveElecs',{'F7','F8'},...        % Electrodes for movement detection
    'moveFeedback','off',...            % Enable movement feedback plots
    'muscle','off',...                 % Disable general muscle detection (first pass)
    'muscleThresh',0.6,...             % General muscle detection threshold
    'muscleFreqWin',[30 100],...       % Frequency window for muscle activity
    'muscleFeedback','on',...          % Enable muscle feedback plots
    'elecNoise','off',...              % Disable electrode noise detection (first pass)
    'elecNoiseThresh',4,...            % Electrode noise threshold
    'elecNoiseFeedback','on' );        % Enable electrode noise feedback

% Save information about removed components
if isfield(EEG, 'reject') && isfield(EEG.reject, 'gcompreject')
    bad_ICAcomp = find(EEG.reject.gcompreject==1);
    save([current_output_folder '\bad_ICAcomp_' current_datasets_savename], 'bad_ICAcomp');
    
    if ~isempty(bad_ICAcomp)
        fprintf('Removed %d ICA components: %s\n', length(bad_ICAcomp), num2str(bad_ICAcomp));
    else
        fprintf('No ICA components selected for removal.\n');
    end
else
    bad_ICAcomp = [];
    fprintf('No ICA components selected for removal.\n');
end

%% 18. Visualization after ICA component removal ---------------------------
% plot and save the figures
visualize_eeg_evoked( ...
        EEG, ...                 % Dataset after ICA component removal
        [-750 0], ...            % Baseline window for visualization
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-50 50], ...            % Amplitude scale: ±50 µV
        'After ICA Ocular Artifact Removal (TESA Method)', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_ICAremoval']); % filename base

% Save dataset after ICA removal
EEG.setname = [current_datasets_savename '_ICAremoval'];
pop_saveset(EEG, 'filename', [current_datasets_savename '_ICAremoval'], ...
    'filepath', current_output_folder);

fprintf('ICA component removal completed.\n');

%% =======================================================================
%  STAGE 6: SSP-SIR FOR MUSCLE ARTIFACT REMOVAL (TESA DEFAULT)
%  =======================================================================

fprintf('\n--- SSP-SIR Muscle Artifact Removal (TESA Default) ---\n');

% Load dataset from previous stage
EEG = pop_loadset('filename', [current_datasets_savename '_ICAremoval.set'], ...
    'filepath', current_output_folder);

% Replace interpolated time with constant 0 values for SSP-SIR
EEG = pop_tesa_removedata(EEG, interp_interval);

fprintf('Testing incremental PC removal (1 to %d components)...\n', max_SSP_SIR_PC_to_remove);

% Create datasets with incremental number of removed components using TESA default
for iii = 1:max_SSP_SIR_PC_to_remove
    fprintf('Processing with %d PC(s) removed...\n', iii);
    
    % Load fresh dataset for each iteration
    EEG_temp = pop_loadset('filename', [current_datasets_savename '_ICAremoval.set'], ...
        'filepath', current_output_folder);
    EEG_temp = pop_tesa_removedata(EEG_temp, interp_interval);
    
    % Apply SSP-SIR using TESA default implementation
    EEG_temp = pop_tesa_sspsir(EEG_temp, ...
        'artScale', 'automatic', ...           % Automatic scaling
        'PC', iii, ...                         % Number of components to remove
        'timeRange', sspsir_timerange);        % Time range for artifact estimation
    
    % Save figures generated by TESA
    if ishandle(gcf)
        saveas(gcf, [current_output_folder '\' current_datasets_savename '_PC' num2str(iii)], 'fig');
        saveas(gcf, [current_output_folder '\' current_datasets_savename '_PC' num2str(iii)], 'jpg');
        close all;
    end
    
    % Save dataset with current number of PCs removed
    EEG_temp.setname = ['After SSP-SIR PC ' num2str(iii)];
    pop_saveset(EEG_temp, 'filename', [current_datasets_savename '_After_SSP-SIR_PC_' num2str(iii) '.set'], ...
        'filepath', current_output_folder);
    
    % Visualization after SSP-SIR
    visualize_eeg_evoked( ...
            EEG_temp, ...            % Dataset after SSP-SIR with current PC removal
            [], ...                  % No baseline correction (will be done in function if needed)
            [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
            [-50 50], ...            % Amplitude scale: ±50 µV
            ['After SSP-SIR (TESA Default) - PC removed: ' num2str(iii)], ...
            true, ...                % save_fig: enable automatic figure saving
            current_output_folder, ... % output_path: directory for saved figures
            [current_datasets_savename '_afterSSP-SIR_PCs_' num2str(iii)]); % filename base
    
    % Create comprehensive visualization with GFP analysis
    EEG_vis = pop_reref(EEG_temp, []);
    EEG_vis = pop_rmbase(EEG_vis, [-100 -2]);
    TEPs = mean(EEG_vis.data, 3);
    
    % Calculate Global Field Power
    gfp = std(mean(EEG_vis.data, 3)).^2;
    gfp_baseline.mean = mean(gfp(:, 4500:4995), 2);
    gfp_baseline.std = std(gfp(:, 4500:4995));
    
    % Find peaks in GFP
    [PKS, LOCS] = findpeaks(gfp, ...
        5000, ...
        'MinPeakHeight', gfp_baseline.mean + (3*gfp_baseline.std), ...
        'MinPeakProminence', 0.75);
    
    LOCS = LOCS * 1000; % Convert to ms
    LOCS = LOCS - 1000; % Adjust zero point
    peaks_after_pulse = find(LOCS > -5 & LOCS < 50);
    LOCS = LOCS(peaks_after_pulse);
    PKS = PKS(peaks_after_pulse);
    
    % Handle edge cases for topoplots
    if length(LOCS) == 1
        final_LOCS = [LOCS LOCS];
        final_PKS = [PKS PKS];
    elseif isempty(LOCS)
        final_LOCS = [59 59];
        final_PKS = [gfp(1059*5) gfp(1059*5)];
    else
        final_LOCS = LOCS;
        final_PKS = PKS;
    end
    
    % Create comprehensive figure with topoplots, GFP, and butterfly plot
    pop_topoplot(EEG_vis, 1000, final_LOCS, 'After SSP-SIR (TESA Default)', [4 length(final_LOCS)], 0, 'electrodes', 'on');
    title('\muV');
    
    % Make window larger
    H = gcf;
    set(H, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    
    % Add GFP subplot
    subplot(3, 1, 2);
    plot(EEG_vis.times(4950:5300), gfp(4950:5300), 'r', 'linewidth', 3);
    hold on;
    title(['GFP, PC: ' num2str(iii)]);
    xlabel('time (ms)');
    ylabel('\muV^2');
    xticks(-20:1:60);
    set(gca, 'Units', 'centimeters', 'TickLength', [0.004, 0.004]);
    plot(final_LOCS, final_PKS, 'b*', 'LineWidth', 2);
    
    % Add butterfly plot
    subplot(3, 1, 3);
    for x = 1:length(EEG_temp.chanlocs)
        plot(EEG_vis.times(4950:5300), TEPs(x, (4950:5300)));
        hold on;
    end
    title(['TEPs -10 +60, PC: ' num2str(iii)]);
    xlabel('time (ms)');
    ylabel('\muV');
    xticks(-20:1:60);
    set(gca, 'Units', 'centimeters', 'TickLength', [0.004, 0.004]);
    
    % Save comprehensive figure
    saveas(gcf, [current_output_folder '\' current_datasets_savename '_merge_' num2str(iii)], 'fig');
    saveas(gcf, [current_output_folder '\' current_datasets_savename '_merge_' num2str(iii)], 'jpg');
    close;
end

fprintf('SSP-SIR processing completed.\n');

%% =======================================================================
%  STAGE 7: FINAL PROCESSING AND RECONSTRUCTION
%  =======================================================================

fprintf('\n--- Final Processing ---\n');

% Manual selection of optimal number of PCs to remove
num_PCs_to_remove = input(['How many components to remove for ' current_datasets_savename ' ? ']);
EEG = pop_loadset('filename', [current_datasets_savename '_After_SSP-SIR_PC_' num2str(num_PCs_to_remove) '.set'], ...
    'filepath', current_output_folder);

fprintf('Selected %d PC(s) for removal.\n', num_PCs_to_remove);

%% 19. Final interpolation of TMS pulse -----------------------------------
EEG = pop_tesa_removedata(EEG, interp_interval);
EEG = pop_tesa_interpdata(EEG, 'cubic', [1,1]);

%% 20. Frequency filtering ------------------------------------------------
EEG = pop_eegfiltnew(EEG, 'hicutoff', low_pass_filt, 'plotfreqz', 1);
EEG = pop_tesa_filtbutter(EEG, notch_filt(1), notch_filt(2), 4, 'bandstop');
close all;

%% 21. Re-referencing to average ------------------------------------------
EEG = pop_reref(EEG, final_ref);

%% 22. Downsampling -------------------------------------------------------
EEG = pop_resample(EEG, downsample);

%% 23. Final baseline correction ------------------------------------------
EEG = pop_rmbase(EEG, baseline_long);

%% 24. Final visualization ------------------------------------------------
% plot and save the figures
visualize_eeg_evoked( ...
        EEG, ...                 % Final processed dataset
        baseline_long, ...       % Baseline window (e.g. [-1000 -2] ms)
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-30 30], ...            % Amplitude scale: -20 to +30 µV (tighter range after processing)
        'Final SOUND + SSP-SIR (TESA Default) Processed TEPs', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_Final_scaled_avg_ref']); % filename base

%% 25. Save final processed dataset ---------------------------------------
EEG.setname = [current_datasets_savename '_Final_SOUND_SSP-SIR_DEFAULT_Processed'];
pop_saveset(EEG, 'filename', [current_datasets_savename '_Final_scaled_avg_ref'], ...
    'filepath', current_output_folder);

%% =======================================================================
%  PIPELINE COMPLETION
%  =======================================================================

diary off; % Stop logging

% Display completion summary
fprintf('\n=== SOUND + SSP-SIR (TESA DEFAULT) PROCESSING COMPLETE ===\n');
fprintf('Processed dataset: %s\n', EEG.setname);
fprintf('Output directory : %s\n', current_output_folder);
fprintf('Final data shape : [%d channels × %d timepoints × %d epochs]\n', ...
        size(EEG.data,1), size(EEG.data,2), size(EEG.data,3));
fprintf('Sampling rate    : %.1f Hz\n', EEG.srate);
fprintf('SSP-SIR PCs removed: %d\n', num_PCs_to_remove);
fprintf('============================================\n');

fprintf('\nProcessing Summary:\n');
fprintf('- Basic preprocessing with TMS pulse detection\n');
fprintf('- SOUND algorithm for channel-specific noise removal (TESA default)\n');
fprintf('- Automatic trial rejection: %d epochs removed\n', num_rejected);
fprintf('- ICA decomposition with %d components removed\n', length(bad_ICAcomp));
fprintf('- SSP-SIR muscle artifact removal with %d PCs (TESA default)\n', num_PCs_to_remove);
fprintf('- Final filtering and reconstruction\n');

fprintf('\nProcessing completed successfully!\n');



