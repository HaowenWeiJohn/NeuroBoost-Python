%==========================================================================
%  STEP-2  ────────────────────────────────────────────────────────────────
%  Basic cleaning:  • detect TMS pulses
%                   • reject noisy channels (“band” / bad electrodes)
%                   • excise & interpolate the -2 … +5 ms pulse artifact
%                   • demean the data
%                   • down-sample to 1000 Hz
%  (c) 2025  –  TMS-EEG-Experiment
%==========================================================================

clear; close all; clc;

PIPELINE_NAME = 'TESA_ICA';          % ← same as step-1
STEP_NAME     = 'step2_basic_processing';

%% 0.  environment --------------------------------------------------------
ENV_PATH  = fullfile('..','.env');
if isfile(ENV_PATH), env = loadenv(ENV_PATH); end     %#ok<*NOPTS>

ROOT_DIR        = getenv('TMS_EEG_ROOT_DIR');
EXPERIMENT_NAME = getenv('EXPERIMENT_NAME');
PARTICIPANT_ID  = getenv('PARTICIPANT_ID');

assert(~isempty(ROOT_DIR) && isfolder(ROOT_DIR), ...
    'ROOT_DIR "%s" not found.',ROOT_DIR);

%% 1.  load Step-1 dataset -------------------------------------------------
in_dir  = fullfile(ROOT_DIR,EXPERIMENT_NAME,PARTICIPANT_ID, ...
                  'output',PIPELINE_NAME);

% the file name was set in step-1:  <EXPERIMENT>-<ID>-<session>-step1_create_set.set
fset = dir(fullfile(in_dir,'*-step1_create_set.set'));
assert(~isempty(fset),'Step-1 dataset not found.');

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
EEG = pop_loadset('filename',fset(1).name,'filepath',in_dir);
EEG = eeg_checkset(EEG);

fprintf('Loaded dataset  : %s\n',EEG.setname);
fprintf('Sampling rate   : %.0f Hz\n',EEG.srate);

%% 2.  detect the TMS pulses ----------------------------------------------
% Uses the largest deflection in a reference channel to mark each pulse.
% Adjust 'chanDetect' if another channel is more reliable.
EEG = pop_tesa_findpulse(EEG, ...
        'rate',EEG.srate, ...        % default search resolution
        'chanDetect','Cz', ...       % reference channel
        'tmsReference','on', ...     % create EEG.event marks
        'plot','off');

%% 3.  reject very noisy (“band”) electrodes ------------------------------
% Kurtosis-based automatic rejection (use robust z-score; threshold = 5 SD)
[EEG, badChans] = pop_rejchan(EEG, ...
        'threshold',5,'norm','on','measure','kurt');
if ~isempty(badChans)
    fprintf('Rejected channels: %s\n', strjoin({EEG.urchanlocs(badChans).labels},', '));
end
EEG = eeg_checkset(EEG);

%% 4.  excise & interpolate the pulse artifact (-2 … +5 ms) ---------------
EEG = pop_tesa_removedata(EEG, ...
        'timelim',[-2 5], ...        % window (ms) relative to each pulse
        'keeptrigs','on', ...        % preserve event markers
        'method','linear');          % cubic also OK; linear safest pre-filter

%% 5.  demean (remove mean over entire recording) -------------------------
EEG = pop_tesa_demean(EEG);

%% 6.  down-sample to 1000 Hz ---------------------------------------------
targetFS = 1000;                     % Hz
if EEG.srate ~= targetFS
    EEG = pop_resample(EEG,targetFS);
end

%% 7.  housekeeping & save -------------------------------------------------
out_dir = fullfile(ROOT_DIR,EXPERIMENT_NAME,PARTICIPANT_ID, ...
                   'output',PIPELINE_NAME);
if ~exist(out_dir,'dir'); mkdir(out_dir); end

save_fname = sprintf('%s-%s-%s-%s.set', ...
                EXPERIMENT_NAME,PARTICIPANT_ID,EEG.session,STEP_NAME);

EEG.setname  = save_fname;
EEG.filename = save_fname;
EEG.filepath = out_dir;
pop_saveset(EEG,'filename',save_fname,'filepath',out_dir);

fprintf('\nStep-2 finished.  Dataset saved to:\n  %s\n',fullfile(out_dir,save_fname));
fprintf('Next pipeline stage ready.\n');
