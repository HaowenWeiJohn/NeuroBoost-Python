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
PIPELINE_NAME = 'SOUND_SSP-SIR';   

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
trigger_original = 16;            % TMS trigger code in EEG data
interp_interval = [-1 6];        % Time window around TMS pulse to interpolate (ms)
                                  % [-1 6] captures direct TMS artifact and recharge

% SOUND algorithm parameters
lambda_value = 0.1;               % Regularization parameter for minimum-norm estimation
                                  % 0.1 was used in the original SOUND paper
iter = 5;                         % Number of iterations for noise estimation
                                  % 5 iterations found sufficient in original studies
ref_channel = 'TP9';              % Reference channel for spherical lead field creation

% SSP-SIR parameters
ssp_sir_timerange = [-1 50];      % Time window for SSP-SIR muscle artifact estimation (ms)
                                  % Should start at interpolation beginning
max_SSP_SIR_PC_to_remove = 6;     % Maximum number of principal components to test
                                  % Creates datasets with 0, 1, 2, ... 5 PCs removed

% Signal processing parameters
downsample = 1000;                % Target sampling rate (Hz) after preprocessing
low_pass_filt = 90;              % Low-pass filter cutoff (Hz)
hp_filt = 1;                     % High-pass filter cutoff (Hz) - removes slow drifts
notch_filt = [48 52];            % Notch filter range (Hz) for line noise removal
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
current_datasets_savename = 'Center_110_part1';   % e.g., stimulation condition/session
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

%% 4. Setup channel locations and remove auxiliary channels ---------------
%    Standardize electrode nomenclature and remove non-EEG channels
EEG = pop_chanedit(EEG, 'lookup','standard-10-5-cap385.elp');
EEG = pop_select(EEG, 'nochannel', {'HEOG', 'VEOG'}); % Remove EOG channels

% Define EEG electrodes and reference
[EEG.chanlocs(1:end).type] = deal('EEG');

%% 5. Create spherical lead field matrix for SOUND algorithm --------------
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

% Save original dataset
EEG.setname = current_datasets_savename;
EEG = pop_saveset(EEG, 'filename', [EEG.setname '_orig'], ...
    'filepath', current_output_folder);

%% 6. High-pass filtering to remove slow drifts ---------------------------
EEG = pop_eegfiltnew(EEG, 'locutoff', hp_filt, 'plotfreqz', 1);
close; % Close filter response plot

%% 7. Epoching around TMS pulses ------------------------------------------
EEG = pop_epoch(EEG, {['S ', num2str(trigger_original)]}, epoching_long, 'epochinfo', 'yes');

%% 8. TMS pulse artifact removal and interpolation ------------------------
EEG_before_interpolation = EEG; % Store dataset before interpolation
EEG = pop_tesa_removedata(EEG, interp_interval); % Remove TMS pulse
EEG = pop_tesa_interpdata(EEG, 'cubic', [1,1]); % Interpolate with cubic spline

%% 9. Baseline correction ------------------------------------------------
EEG = pop_rmbase(EEG, baseline_long);

%% 10. Re-reference lead field matrix to specified reference -------------
ref_chan_num = find(strcmp({EEG_dummy.chanlocs.labels}, ref_channel));
LFM_sphere_ref = repmat(LFM_sphere(ref_chan_num,:), size(LFM_sphere, 1)-1, 1);
LFM_sphere = LFM_sphere([1:ref_chan_num-1 ref_chan_num+1:end],:) - LFM_sphere_ref;

%% 11. Visualization after basic preprocessing ----------------------------
EEG_epoch = EEG;
EEG_epoch = pop_rmbase(EEG_epoch, [-750 0]);
TEP_dummy = mean(EEG_epoch.data, 3);

figure
plot(EEG_epoch.times, TEP_dummy(:,:));
xlim([-100 400]);
ylim([-50 50]);
title('After Basic Pre-processing');
saveas(gcf, [current_output_folder '\' current_datasets_savename '_basic_prep']);
saveas(gcf, [current_output_folder '\' current_datasets_savename '_basic_prep'], 'png');
close;

% Save dataset and lead field
EEG.setname = [current_datasets_savename '_basic_prep'];
pop_saveset(EEG, 'filename', [current_datasets_savename '_basic_prep'], ...
    'filepath', current_output_folder);
save([current_output_folder '\LFM_sphere'], 'LFM_sphere');

fprintf('Basic preprocessing completed.\n');

%% =======================================================================
%  STAGE 2: SOUND ALGORITHM FOR CHANNEL-SPECIFIC NOISE REMOVAL
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
EEG_epoch = EEG_clean;
EEG_epoch = pop_rmbase(EEG_epoch, [-750 0]);
TEP_dummy = mean(EEG_epoch.data, 3);

figure
plot(EEG_epoch.times, TEP_dummy(:,:));
xlim([-100 400]);
ylim([-50 50]);
title('After SOUND Algorithm');
saveas(gcf, [current_output_folder '\' current_datasets_savename '_afterSOUND']);
saveas(gcf, [current_output_folder '\' current_datasets_savename '_afterSOUND'], 'png');
close all;

% Save dataset after SOUND
EEG = EEG_clean;
EEG.setname = [current_datasets_savename '_afterSOUND'];
pop_saveset(EEG, 'filename', [current_datasets_savename '_afterSOUND'], ...
    'filepath', current_output_folder);

fprintf('SOUND algorithm completed.\n');

%% =======================================================================
%  STAGE 3: MANUAL ARTIFACT REJECTION
%  =======================================================================

fprintf('\n--- Manual Artifact Rejection ---\n');

% Load dataset from previous stage
EEG = pop_loadset('filename', [current_datasets_savename '_afterSOUND.set'], ...
    'filepath', current_output_folder);

% Manual epoch rejection with visual inspection
fprintf('Opening manual artifact rejection interface...\n');
fprintf('Instructions: Click on bad epochs (they turn yellow), then press UPDATE MARKS\n');

eeglab redraw
close all
pop_eegplot(EEG, 1, 1, 0); % Plot data for manual inspection
waitfor(findobj('parent', gcf, 'string', 'UPDATE MARKS'), 'userdata');
close all

% Save indices of rejected epochs
rejected_epochs = find(EEG.reject.rejmanual);
fprintf('Rejected epochs: %s\n', num2str(rejected_epochs));
save([current_output_folder '\rejected_epochs'], 'rejected_epochs');

% Remove marked epochs
EEG = pop_rejepoch(EEG, EEG.reject.rejmanual, 0);

%% 13. Visualization after trial rejection -------------v------------------
EEG_epoch = EEG;
EEG_epoch = pop_rmbase(EEG_epoch, [-750 0]);
TEP_dummy = mean(EEG_epoch.data, 3);

figure
plot(EEG_epoch.times, TEP_dummy(:,:));
xlim([-100 400]);
ylim([-50 50]);
title('After Trial Rejection');
saveas(gcf, [current_output_folder '\' current_datasets_savename '_trialrej']);
saveas(gcf, [current_output_folder '\' current_datasets_savename '_trialrej'], 'png');
close all;

% Save dataset after trial rejection
EEG.setname = [current_datasets_savename '_trialrej'];
pop_saveset(EEG, 'filename', [current_datasets_savename '_trialrej'], ...
    'filepath', current_output_folder);

fprintf('Manual artifact rejection completed.\n');

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
EEG_epoch = EEG;
EEG_epoch = pop_rmbase(EEG_epoch, [-750 0]);
TEP_dummy = mean(EEG_epoch.data, 3);

figure
plot(EEG_epoch.times, TEP_dummy(:,:));
xlim([-100 400]);
ylim([-50 50]);
title('After ICA Removal');
saveas(gcf, [current_output_folder '\' current_datasets_savename '_ICAremoval']);
saveas(gcf, [current_output_folder '\' current_datasets_savename '_ICAremoval'], 'png');
close all;

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
    EEG_epoch = EEG;
    TEP_dummy = mean(EEG_epoch.data, 3);
    
    figure
    plot(EEG_epoch.times, TEP_dummy(:,:));
    xlim([-100 400]);
    ylim([-50 50]);
    title(['After SSP-SIR - PC removed: ' num2str(iii-1)]);
    saveas(gcf, [current_output_folder '\' current_datasets_savename '_afterSSP-SIR_PCs_' num2str(iii-1)]);
    saveas(gcf, [current_output_folder '\' current_datasets_savename '_afterSSP-SIR_PCs_' num2str(iii-1)], 'png');
    close all;
    
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
figure
plot(EEG.times, mean(EEG.data, 3));
xlim([-100 400]);
ylim([-20 30]);
title('Final SOUND + SSP-SIR Processed TEPs');
xlabel('Time (ms)');
ylabel('Amplitude (\muV)');
saveas(gcf, [current_output_folder '\' current_datasets_savename '_Final_scaled_avg_ref']);
saveas(gcf, [current_output_folder '\' current_datasets_savename '_Final_scaled_avg_ref'], 'png');
close all;

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
