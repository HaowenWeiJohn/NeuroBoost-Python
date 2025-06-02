clear; close all; clc;

% Define the pipeline name for tracking/logging purposes
PIPELINE_NAME = 'TESA_ICA';   

%% Variables
activity_prob=5; %[max] absolute thresold or activity probability limit(s) (in std. dev.) for automatic channel rejection
epoching_long=[-1 1]; %epoching in sec
demeaning_interval=[-1000 999]; %interval in ms for demeaning
downsample=1000; %Hz to reach for the downsampling
tr_rej_stnddev_first=5; %number of stnd dev to consider when computing joint probability and remove trial (trial rejection)
baseline_long=[-1000 -2]; %full period used for baseline correction
trigger_label= 'TMS'; 
interp_interval=[-1 6]; %interval where to cut and interpolate TMS pulse
low_pass_filt=90; %low pass filter limit in Hz
hp_filt=1; %high pass filter cut-off fr in Hz
notch_filt=[58 62]; %notch filter limit in Hz % 60 Hz in US
final_ref=[]; %reref for the last image. [] for avg reref


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

% -------------------------------------------------------------------------
%    Define the dataset‑save name (e.g., stimulation condition)
current_datasets_savename = 'Pos10_80';   % ← adjust if you change file names

% -------------------------------------------------------------------------

%    Construct the full path to the BrainVision header file (.vhdr)
vhdr_file = fullfile(data_root, [current_datasets_savename '.vhdr']);
assert(isfile(vhdr_file), ...
       'Session file "%s" does not exist.', vhdr_file);

current_output_folder = char( fullfile( ...
        ROOT_DIR, EXPERIMENT_NAME, PARTICIPANT_ID, ...
        'output', PIPELINE_NAME, current_datasets_savename) );
if ~exist(current_output_folder, 'dir'); mkdir(current_output_folder); end

%% 3. Load EEG data into EEGLAB -------------------------------------------
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
EEG = pop_loadbv(data_root, [current_datasets_savename '.vhdr']);
EEG = eeg_checkset(EEG);
EEG.setname = sprintf('%s_%s_%s', ...
        EXPERIMENT_NAME, PARTICIPANT_ID, current_datasets_savename);

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


%% Create a copy in the dataset in the intermediate folder
pop_saveset(...
    EEG,...
    'filename', [current_datasets_savename '_SetCreation.set'],...
    'filepath', current_output_folder);

%% find TMS pulse
EEG = pop_tesa_findpulse(EEG,'Cz', ...
                         'refract',4, ...
                         'rate',1e4, ...
                         'tmsLabel', trigger_label, ...
                         'plots','on');
EEG = eeg_checkset(EEG);


EEG.allchan = EEG.chanlocs;
EEG = pop_rejchan(EEG, 'elec', 1:size(EEG.data,1) ,'threshold',activity_prob,'norm','on','measure','kurt');

%change set name and save this step on disk
EEG.setname = [current_datasets_savename '_ChanRem'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);


%% Epoch data
% epoch data around the trigger S127 (TMS pulse) , -1000ms
% +1000ms
EEG = pop_epoch( EEG, { trigger_label }, epoching_long , 'epochinfo', 'yes');

%change set name and save this step on disk
EEG.setname = [current_datasets_savename '_epochs'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% Demeaning Correction
%apply demeaning, subtraction of the mean voltage of the whole
%epoch to each point in the epoch. This steps increase the
%reliability of the ICA
EEG = pop_rmbase( EEG, demeaning_interval);

%change set name and save this step on disk
EEG.setname = [current_datasets_savename '_Demeaning'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);


%% Remove TMS pulse artifact
% remove data around the TMS pulse
EEG = pop_tesa_removedata( EEG, interp_interval);

%change set name and save this step on disk
EEG.setname = [current_datasets_savename '_TMSpulseREM'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% Interpolate missing data around TMS pulse
% interpolation of missing data aroun TMS pulse using a cubic
% interpolation and 1 ms before an after the removed window
EEG = pop_tesa_interpdata( EEG, 'cubic', [1,1] );

%change set name and save this step on disk
EEG.setname = [current_datasets_savename '_TMSpInt'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% Downsample data
% now that the TMS pulse artifact is removed it is possible to
% demeaning without inserting aliasing artifacts
EEG = pop_resample( EEG, downsample);

EEG.setname = [current_datasets_savename '_DownSamp'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% Remove bad trials
EEG = pop_jointprob(EEG,1, 1:size(EEG.data,1),tr_rej_stnddev_first,tr_rej_stnddev_first,0,0); %automatic selection of bad epochs
disp(['number of rejected epochs: ' find(EEG.reject.rejjp)]) %print numebr of rejected epochs
bad_epochs=find(EEG.reject.rejjp); %create a new variables containing the indeces of bad epochs and store it in the next step
save([current_output_folder '\bad_epochs_' current_datasets_savename], 'bad_epochs');
EEG = pop_rejepoch( EEG, EEG.reject.rejjp, 0); %reject marked epochs

%% Replace interpolated data around TMS pulse with constant amplitude data (-1 to 6 ms)
% it is because interpolated data should not be present when
% performing an ICA (because it is like you are adding
% information to the dataset)
EEG = pop_tesa_removedata( EEG, interp_interval );

EEG.setname = [current_datasets_savename '_TMSpulse0'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% Remove TMS-evoked muscle activity (using FastICA and semi-auto component selection)
EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );

EEG.setname = [current_datasets_savename '_ICA-TMSmuscle'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);


%parameters for automatic component selection

EEG = pop_tesa_compselect( EEG, ...
    'compCheck','off',...
    'comps', 15, ...
    'figSize','small',...
    'plotTimeX',[-100 399],...
    'plotFreqX',[1 100],...
    'tmsMuscle','on',...
    'tmsMuscleThresh',8,...
    'tmsMuscleWin',[11 30],...
    'tmsMuscleFeedback','off',...
    'blink','off',...
    'blinkThresh',2.5,...
    'blinkElecs',{'Fp1','Fp2'},...
    'blinkFeedback','off',...
    'move','off',...
    'moveThresh',2,...
    'moveElecs',{'F7','F8'},...
    'moveFeedback','off',...
    'muscle','off',...
    'muscleThresh',0.6,...
    'muscleFreqWin',[30 100],...
    'muscleFeedback','off',...
    'elecNoise','off',...
    'elecNoiseThresh',4,...
    'elecNoiseFeedback','off' );


EEG.setname = [current_datasets_savename '_CompSel1'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);


%% load dataset afeter first component selection

EEG = pop_loadset(...
    'filename', [current_datasets_savename  '_CompSel1.set'],...
    'filepath',  current_output_folder);

%%  data removal (-1 to 6 ms)
% this step is suggested from the TESA script example, even if
% i don't understand why it is needed. I suppose is to be sure
% there are high frequency artifact at the edge of the removed
% window
EEG = pop_tesa_removedata( EEG, interp_interval );

EEG.setname = [current_datasets_savename '_TMSp0ext'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);


%% Interpolate missing data around TMS pulse
%interpolate back remove data that were substituted with
%costant values. This is because for filtering it is better to
%do not have sharp edges in the epoch.
EEG = pop_tesa_interpdata( EEG, 'cubic', [5,5] );

EEG.setname = [current_datasets_savename '_TMSp0Int'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);


%% Bandpass (1-90 Hz) and bandstop (48-52 Hz) filter data
EEG = pop_tesa_filtbutter( EEG, hp_filt, low_pass_filt, 4, 'bandpass' );
EEG = pop_tesa_filtbutter( EEG, notch_filt(1), notch_filt(2), 4, 'bandstop' );

EEG.setname = [current_datasets_savename '_notch&90_1_BandpassFilter'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% Replace interpolated data around TMS pulse with constant amplitude data (-1 to 6 ms)
EEG = pop_tesa_removedata( EEG, interp_interval );

EEG.setname = [current_datasets_savename '_TMSp0'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% Remove all other artifact (using FastICA and semi-auto component selection)
EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );

EEG.setname = [current_datasets_savename '_ICA2'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

EEG = pop_tesa_compselect( EEG,...
    'compCheck','off',...
    'comps',[],...
    'figSize','medium',...
    'plotTimeX',[-100 399],...
    'plotFreqX',[1 100],...
    'tmsMuscle','on',...
    'tmsMuscleThresh',8,...
    'tmsMuscleWin',[11 30],...
    'tmsMuscleFeedback','off',...
    'blink','on',...
    'blinkThresh',2.5,...
    'blinkElecs',{'Fp1','Fp2'},...
    'blinkFeedback','off',...
    'move','on',...
    'moveThresh',2,...
    'moveElecs',{'F7','F8'},...
    'moveFeedback','off',...
    'muscle','on',...
    'muscleThresh',0.6,...
    'muscleFreqWin',[30 100],...
    'muscleFeedback','off',...
    'elecNoise','on',...
    'elecNoiseThresh',4,...
    'elecNoiseFeedback','off' );


EEG.setname = [current_datasets_savename '_CompSel2'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);


%% load dataset afeter second component selection
EEG = pop_loadset(...
    'filename', [current_datasets_savename  '_CompSel2.set'],...
    'filepath', current_output_folder);

%% Interpolate missing data around TMS pulse
EEG = pop_tesa_interpdata( EEG, 'cubic', [5,5] );

EEG.setname = [current_datasets_savename '_After2ndCompSel_TMSpInt'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% Interpolate missing channels
%channel removed because were bad are now interpolated
EEG = pop_interp(EEG, EEG.allchan, 'spherical');

EEG.setname = [current_datasets_savename '_ChanInt'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% Re-reference to average of all electrodes
%Data were recorded with online reference in Tp9.
EEG = pop_reref( EEG, final_ref);

EEG.setname = [current_datasets_savename '_AvgReref'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% Baseline Correction (-1000 ms to -2 ms)
EEG = pop_rmbase( EEG, baseline_long);

EEG.setname = [current_datasets_savename '_baselineCorr'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);


%% Save point
EEG.setname= [current_datasets_savename 'TESA_Processed'];
EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);

%% Plot the results
% final TEP  scaled

% – Baseline: use the same window you defined earlier (baseline_long)
% – Show only –100…400 ms, scale –20…30 µV, give the plot a clear title
visualize_eeg_evoked( ...
        EEG, ...                 % pre-processed dataset
        baseline_long, ...       % e.g. [-1000 -2] ms
        [-100 400], ...          % x-axis limits
        [-20 30], ...            % y-axis limits
        'Final scaled avg reref');

%% Save the figure
outfile = fullfile(current_output_folder, ...
                   [current_datasets_savename '_Final_scaled_avg_ref']);

saveas(gcf, [outfile '.fig']);   % MATLAB-figure file
saveas(gcf, [outfile '.png']);   % high-res PNG

close(gcf);          % tidy up
diary off;           % stop logging (harmless if diary wasn’t on)



