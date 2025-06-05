clear; close all; clc;

% ===========================================================================
% SOUND + SSP-SIR TMS-EEG Processing Pipeline
% ===========================================================================
% This pipeline implements an advanced TMS-EEG preprocessing workflow
% combining SOUND (Suppression of Unwanted Noise and Disturbances) and
% SSP-SIR (Signal Space Separation - Source-Informed Reconstruction) algorithms
% for state-of-the-art artifact removal and signal cleaning.
% 
% Key processing steps:
% 1. Data loading and basic preprocessing
% 2. Lead field matrix creation for source modeling
% 3. SOUND algorithm for channel-specific noise removal
% 4. Manual artifact rejection with visual inspection
% 5. ICA decomposition for ocular artifact removal
% 6. SSP-SIR for muscle artifact removal with incremental testing
% 7. Final processing and signal reconstruction
%
% Pipeline optimized for research-grade TMS-EEG data analysis with maximum
% flexibility and comprehensive quality control.
% ===========================================================================

% Define the pipeline name for tracking/logging purposes
PIPELINE_NAME = 'SOUND-SSPSIR';   

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

% SOUND algorithm parameters
lambda_value = 0.1;               % Regularization parameter for minimum-norm estimation
                                  % 0.1 was used in the original SOUND paper
iter = 5;                         % Number of iterations for noise estimation
                                  % 5 iterations found sufficient in original studies
ref_channel = 'C4';              % Reference channel for spherical lead field creation

% SSP-SIR parameters
ssp_sir_timerange = [-1 50];      % Time window for SSP-SIR muscle artifact estimation (ms)
                                  % Should start at interpolation beginning
max_SSP_SIR_PC_to_remove = 6;     % Maximum number of principal components to test
                                  % Creates datasets with 0, 1, 2, ... 5 PCs removed

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
diary([current_output_folder '\' current_datasets_savename '_workspace_SOUND_SSP-SIR']);

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

%% 7. Create spherical lead field matrix for SOUND algorithm --------------
%    This matrix allows conversion between sensor and source space
fprintf('\n--- Creating Lead Field Matrix ---\n');

% Create dummy EEGLAB structure for constructing leadfield with correct reference
EEG_dummy = EEG;
EEG_dummy.data(63,length(EEG_dummy.data)) = 0; % Add reference channel data
EEG_dummy.nbchan = 63; % Update channel count
EEG_dummy.chanlocs(63).type = 'EEG';
EEG_dummy.chanlocs(63).labels = ref_channel;
EEG_dummy = pop_chanedit(EEG_dummy, 'lookup', 'standard-10-5-cap385.elp');

% Construct spherical lead field matrix
[LFM_sphere] = construct_spherical_lead_field(EEG_dummy, [],[],[],[],[],[],[],1000);

%% 8. High-pass filtering to remove slow drifts ---------------------------
EEG = pop_eegfiltnew(EEG, 'locutoff', hp_filt, 'plotfreqz', 1);
close; % Close filter response plot

%% 9. Epoching around TMS pulses ------------------------------------------
EEG = pop_epoch(EEG, {trigger_label}, epoching_long, 'epochinfo', 'yes');

%% 10. TMS pulse artifact removal and interpolation ------------------------
EEG_before_interpolation = EEG; % Store dataset before interpolation
EEG = pop_tesa_removedata(EEG, interp_interval); % Remove TMS pulse
EEG = pop_tesa_interpdata(EEG, 'cubic', [1,1]); % Interpolate with cubic spline

%% 11. Baseline correction ------------------------------------------------
EEG = pop_rmbase(EEG, baseline_long);

%% 12. Re-reference lead field matrix to specified reference -------------
ref_chan_num = find(strcmp({EEG_dummy.chanlocs.labels}, ref_channel));
LFM_sphere_ref = repmat(LFM_sphere(ref_chan_num,:), size(LFM_sphere, 1)-1, 1);
LFM_sphere = LFM_sphere([1:ref_chan_num-1 ref_chan_num+1:end],:) - LFM_sphere_ref;

%% 13. Visualization after basic preprocessing ----------------------------
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

% Save dataset and lead field
EEG.setname = [current_datasets_savename '_basic_prep'];
pop_saveset(EEG, 'filename', [current_datasets_savename '_basic_prep'], ...
    'filepath', current_output_folder);
save([current_output_folder '\LFM_sphere'], 'LFM_sphere');

fprintf('Basic preprocessing completed.\n');

%% =======================================================================
%  STAGE 3: SOUND ALGORITHM FOR CHANNEL-SPECIFIC NOISE REMOVAL
%  =======================================================================

fprintf('\n--- Starting SOUND Algorithm ---\n');

% Load dataset from previous stage
EEG = pop_loadset('filename', [current_datasets_savename '_basic_prep.set'], ...
    'filepath', current_output_folder);
load([current_output_folder '\LFM_sphere'], 'LFM_sphere');

EEG_clean = EEG;
chanN = size(EEG.data, 1);

% Initialize storage for noise estimates
all_sigmas = zeros(size(EEG.data, 1), size(EEG.data, 2));

fprintf('Processing %d time samples with SOUND algorithm...\n', size(EEG.data, 2));
tic

% Parallel processing for each time sample
parfor SOUND_maincounter = 1:(size(EEG.data, 2))
    % Find the number of trials
    N_tr = size(EEG.data, 3);
    
    % Inform about progress
    if mod(SOUND_maincounter, 500) == 0
        fprintf('Processing time sample %d of %d\n', SOUND_maincounter, size(EEG.data, 2));
    end
    
    % Reshape data for current time sample
    tmp_data = reshape(EEG.data(:, SOUND_maincounter, :), [EEG.nbchan, N_tr]);
    
    % Run DDWiener for initial noise estimation
    [y_solved, sigmas] = DDWiener(tmp_data);
    
    % Run SOUND_alt for refined noise estimation
    [~, ~, sigmas] = SOUND_alt(tmp_data, LFM_sphere, iter, lambda_value, sigmas, 1);
    
    all_sigmas(:, SOUND_maincounter) = sigmas;
end

% Apply temporal smoothing to noise estimates
fprintf('Applying temporal smoothing to noise estimates...\n');
scaling_w = 2;
smoothness = 1;
all_sigmas_ext = [repmat(all_sigmas(:,1), [1, scaling_w]), all_sigmas, ...
                  repmat(all_sigmas(:,end), [1, scaling_w])];
step = scaling_w + 1;

% Create Gaussian smoothing function
smoothing_func = gaussmf(1:(2*scaling_w + 1), [smoothness, scaling_w + 1]);
smoothing_func = smoothing_func / sum(smoothing_func);
figure; plot(smoothing_func); title('SOUND Temporal Smoothing Function');
saveas(gcf, [current_output_folder '\' current_datasets_savename '_smoothing_func'], 'png');
close;
smoothing_func = repmat(smoothing_func, [size(EEG.data, 1), 1]);

% Apply SOUND cleaning with smoothed noise estimates
for SOUND_maincounter = 1:size(EEG.data, 2)
    N_tr = size(EEG.data, 3);
    
    % Get smoothed noise estimate for this sample
    sigmas = sum(all_sigmas_ext(:, (step-scaling_w):(step+scaling_w)) .* smoothing_func, 2);
    
    % Apply SOUND correction
    W = diag(1./sigmas);
    WL = W * LFM_sphere;
    WLLW = WL * WL';
    
    data = reshape(EEG.data(:, SOUND_maincounter, :), [chanN, N_tr]);
    
    x = WL' * ((WLLW + lambda_value*trace(WLLW)/chanN*eye(chanN)) \ (W*data));
    corrected_data = LFM_sphere * x;
    EEG_clean.data(:, SOUND_maincounter, :) = reshape(corrected_data, [size(corrected_data,1), 1, N_tr]);
    
    step = step + 1;
end

fprintf('SOUND processing completed in %.2f minutes.\n', toc/60);

%% 12. Visualization after SOUND ------------------------------------------
% plot and save the figures
visualize_eeg_evoked( ...
        EEG_clean, ...           % Dataset after SOUND algorithm
        [-750 0], ...            % Baseline window for visualization
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-50 50], ...            % Amplitude scale: ±50 µV
        'After SOUND Algorithm', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_afterSOUND']); % filename base

% Save dataset after SOUND
EEG = EEG_clean;
EEG.setname = [current_datasets_savename '_afterSOUND'];
pop_saveset(EEG, 'filename', [current_datasets_savename '_afterSOUND'], ...
    'filepath', current_output_folder);

fprintf('SOUND algorithm completed.\n');

%% =======================================================================
%  STAGE 4: ICA DECOMPOSITION FOR OCULAR ARTIFACTS
%  =======================================================================

fprintf('\n--- ICA Decomposition ---\n');

% Load dataset from previous stage
EEG = pop_loadset('filename', [current_datasets_savename '_trialrej.set'], ...
    'filepath', current_output_folder);

% ICA decomposition
fprintf('Running ICA decomposition...\n');
EEG = pop_runica(EEG, 'chanind', 1:length(EEG.chanlocs), 'interupt', 'on');

% Save dataset after ICA
EEG.setname = [current_datasets_savename '_After_ICA'];
pop_saveset(EEG, 'filename', [current_datasets_savename '_After_ICA'], ...
    'filepath', current_output_folder);

fprintf('ICA decomposition completed.\n');

%% 14. Manual ICA component selection for ocular artifacts ----------------
fprintf('\n--- Manual ICA Component Selection ---\n');

% Load dataset with ICA weights
EEG = pop_loadset('filename', [current_datasets_savename '_After_ICA.set'], ...
    'filepath', current_output_folder);

eeglab redraw
close all

% Interactive component selection
fprintf('Instructions: Select ocular components to remove, then press OK\n');
EEG = pop_selectcomps(EEG, 1:size(EEG.icaweights, 1));
waitfor(findobj('parent', figure(2), 'string', 'OK'), 'userdata');
close all;

% Save bad ICA components
EEG.BadIC = find(EEG.reject.gcompreject==1);
bad_ICAcomp = EEG.BadIC;
save([current_output_folder '\bad_ICAcomp_' current_datasets_savename], 'bad_ICAcomp');

% Remove selected components
if ~isempty(EEG.BadIC)
    EEG = pop_subcomp(EEG, EEG.BadIC, 0);
    fprintf('Removed %d ICA components: %s\n', length(EEG.BadIC), num2str(EEG.BadIC));
else
    fprintf('No ICA components selected for removal.\n');
end

%% 15. Visualization after ICA removal ------------------------------------
% plot and save the figures
visualize_eeg_evoked( ...
        EEG, ...                 % Dataset after ICA component removal
        [-750 0], ...            % Baseline window for visualization
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-50 50], ...            % Amplitude scale: ±50 µV
        'After ICA Removal', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_ICAremoval']); % filename base

% Save dataset after ICA removal
EEG.setname = [current_datasets_savename '_ICAremoval'];
pop_saveset(EEG, 'filename', [current_datasets_savename '_ICAremoval'], ...
    'filepath', current_output_folder);

fprintf('ICA component removal completed.\n');

%% =======================================================================
%  STAGE 5: SSP-SIR FOR MUSCLE ARTIFACT REMOVAL
%  =======================================================================

fprintf('\n--- SSP-SIR Muscle Artifact Removal ---\n');

% Load dataset from previous stage
EEG = pop_loadset('filename', [current_datasets_savename '_ICAremoval.set'], ...
    'filepath', current_output_folder);

% Replace interpolated time with constant 0 values
EEG = pop_tesa_removedata(EEG, interp_interval);

% Load lead field matrix
load([current_output_folder '\LFM_sphere'], 'LFM_sphere');

EEG_clean_muscle = EEG;
EEG_clean_ica = EEG;

% Calculate rank for SSP-SIR
dimN = rank(mean(EEG_clean_ica.data(1:end,:,:), 3));

% Storage for PC variance tracking
list_PCs_variance = zeros(1, max_SSP_SIR_PC_to_remove);

fprintf('Testing incremental PC removal (0 to %d components)...\n', max_SSP_SIR_PC_to_remove-1);

% Create datasets with incremental number of removed components
for iii = 1:max_SSP_SIR_PC_to_remove
    fprintf('Processing with %d PC(s) removed...\n', iii-1);
    
    % Run SSP_SIR_GB with current number of components
    [~, artifact_topographies, ~, ~, filt_ker, SSP_SIR_operator, PCs_variance] = SSP_SIR_GB(...
        mean(EEG_clean_ica.data(1:end,:,:), 3), ...
        iii, ...
        LFM_sphere, ...
        [], ...
        [], ...
        0, ...
        EEG_clean_ica.srate, ...
        EEG.times, ...
        ssp_sir_timerange);
    
    list_PCs_variance(iii) = PCs_variance;
    
    % Save SSP-SIR component figures
    saveas(gcf, [current_output_folder '\' current_datasets_savename '_PC' num2str(iii-1)], 'fig');
    saveas(gcf, [current_output_folder '\' current_datasets_savename '_PC' num2str(iii-1)], 'jpg');
    close all;
    
    % Apply SSP-SIR correction to all trials
    for i = 1:size(EEG_clean_ica.data, 3)
        EEG_clean_muscle.data(1:end,:,i) = filt_ker.*(SSP_SIR_operator*EEG_clean_ica.data(1:end,:,i)) + ...
            EEG_clean_ica.data(1:end,:,i) - filt_ker.*EEG_clean_ica.data(1:end,:,i);
    end
    
    % Save topographies of removed components
    EEG_clean_muscle.SSPSIR_topographies = artifact_topographies;
    
    % Save dataset with current number of PCs removed
    EEG = EEG_clean_muscle;
    EEG.setname = 'After SSP-SIR';
    EEG = pop_saveset(EEG, 'filename', [current_datasets_savename '_After_SSP-SIR_PC_' num2str(iii-1) '.set'], ...
        'filepath', current_output_folder);
    
    % Visualization after SSP-SIR
    visualize_eeg_evoked( ...
            EEG, ...                 % Dataset after SSP-SIR with current PC removal
            [], ...                  % No baseline correction (will be done in function if needed)
            [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
            [-50 50], ...            % Amplitude scale: ±50 µV
            ['After SSP-SIR - PC removed: ' num2str(iii-1)], ...
            true, ...                % save_fig: enable automatic figure saving
            current_output_folder, ... % output_path: directory for saved figures
            [current_datasets_savename '_afterSSP-SIR_PCs_' num2str(iii-1)]); % filename base
    
    % Create comprehensive visualization with GFP analysis
    EEG_vis = pop_reref(EEG, []);
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
    pop_topoplot(EEG_vis, 1000, final_LOCS, 'After SSP-SIR', [4 length(final_LOCS)], 0, 'electrodes', 'on');
    title('\muV');
    
    % Make window larger
    H = gcf;
    set(H, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    
    % Add GFP subplot
    subplot(3, 1, 2);
    plot(EEG_vis.times(4950:5300), gfp(4950:5300), 'r', 'linewidth', 3);
    hold on;
    title(['GFP, PC: ' num2str(iii-1)]);
    xlabel('time (ms)');
    ylabel('\muV^2');
    xticks(-20:1:60);
    set(gca, 'Units', 'centimeters', 'TickLength', [0.004, 0.004]);
    plot(final_LOCS, final_PKS, 'b*', 'LineWidth', 2);
    
    % Add butterfly plot
    subplot(3, 1, 3);
    for x = 1:length(EEG.chanlocs)
        plot(EEG_vis.times(4950:5300), TEPs(x, (4950:5300)));
        hold on;
    end
    title(['TEPs -10 +60, PC: ' num2str(iii-1)]);
    xlabel('time (ms)');
    ylabel('\muV');
    xticks(-20:1:60);
    set(gca, 'Units', 'centimeters', 'TickLength', [0.004, 0.004]);
    
    % Save comprehensive figure
    saveas(gcf, [current_output_folder '\' current_datasets_savename '_merge_' num2str(iii-1)], 'fig');
    saveas(gcf, [current_output_folder '\' current_datasets_savename '_merge_' num2str(iii-1)], 'jpg');
    close;
end

% Display PC variance summary
fprintf('\nSSP-SIR PC Variance Summary:\n');
for i = 1:max_SSP_SIR_PC_to_remove
    fprintf('PC %d removed: %.2f variance\n', i-1, list_PCs_variance(i));
end

fprintf('SSP-SIR processing completed.\n');

%% =======================================================================
%  STAGE 6: FINAL PROCESSING AND RECONSTRUCTION
%  =======================================================================

fprintf('\n--- Final Processing ---\n');

% Manual selection of optimal number of PCs to remove
num_PCs_to_remove = input(['How many components to remove for ' current_datasets_savename ' ? ']);
EEG = pop_loadset('filename', [current_datasets_savename '_After_SSP-SIR_PC_' num2str(num_PCs_to_remove) '.set'], ...
    'filepath', current_output_folder);

fprintf('Selected %d PC(s) for removal.\n', num_PCs_to_remove);

%% 16. Final interpolation of TMS pulse -----------------------------------
EEG = pop_tesa_removedata(EEG, interp_interval);
EEG = pop_tesa_interpdata(EEG, 'cubic', [1,1]);

%% 17. Frequency filtering ------------------------------------------------
EEG = pop_eegfiltnew(EEG, 'hicutoff', low_pass_filt, 'plotfreqz', 1);
EEG = pop_tesa_filtbutter(EEG, notch_filt(1), notch_filt(2), 4, 'bandstop');
close all;

%% 18. Re-referencing to average ------------------------------------------
EEG = pop_reref(EEG, final_ref);

%% 19. Downsampling -------------------------------------------------------
EEG = pop_resample(EEG, downsample);

%% 20. Final baseline correction ------------------------------------------
EEG = pop_rmbase(EEG, baseline_long);

%% 21. Final visualization ------------------------------------------------
% plot and save the figures
visualize_eeg_evoked( ...
        EEG, ...                 % Final processed dataset
        baseline_long, ...       % Baseline window (e.g. [-1000 -2] ms)
        [-100 400], ...          % Display window: 100 ms pre to 400 ms post-TMS
        [-20 30], ...            % Amplitude scale: -20 to +30 µV (tighter range after processing)
        'Final SOUND + SSP-SIR Processed TEPs', ...
        true, ...                % save_fig: enable automatic figure saving
        current_output_folder, ... % output_path: directory for saved figures
        [current_datasets_savename '_Final_scaled_avg_ref']); % filename base

%% 22. Save final processed dataset ---------------------------------------
EEG.setname = [current_datasets_savename '_Final_SOUND_SSP-SIR_Processed'];
pop_saveset(EEG, 'filename', [current_datasets_savename '_Final_scaled_avg_ref'], ...
    'filepath', current_output_folder);

%% =======================================================================
%  PIPELINE COMPLETION
%  =======================================================================

diary off; % Stop logging

% Display completion summary
fprintf('\n=== SOUND + SSP-SIR PROCESSING COMPLETE ===\n');
fprintf('Processed dataset: %s\n', EEG.setname);
fprintf('Output directory : %s\n', current_output_folder);
fprintf('Final data shape : [%d channels × %d timepoints × %d epochs]\n', ...
        size(EEG.data,1), size(EEG.data,2), size(EEG.data,3));
fprintf('Sampling rate    : %.1f Hz\n', EEG.srate);
fprintf('SSP-SIR PCs removed: %d\n', num_PCs_to_remove);
fprintf('============================================\n');

fprintf('\nProcessing Summary:\n');
fprintf('- Basic preprocessing with lead field creation\n');
fprintf('- SOUND algorithm for channel-specific noise removal\n');
fprintf('- Manual artifact rejection: %d epochs removed\n', length(rejected_epochs));
fprintf('- ICA decomposition with %d components removed\n', length(bad_ICAcomp));
fprintf('- SSP-SIR muscle artifact removal with %d PCs\n', num_PCs_to_remove);
fprintf('- Final filtering and reconstruction\n');

fprintf('\nProcessing completed successfully!\n');
