%==========================================================================
% STEP-2  ·  Basic processing with plots
%   1  Detect TMS pulses                (pop_tesa_findpulse)
%   2  Reject noisy “band” electrodes   (pop_rejchan)
%   3  Epoch  –1 … +1 s around pulses   (pop_epoch)
%   4  Remove –2 … +5 ms pulse window   (pop_tesa_removedata + interp)
%   5  Demean each epoch                (pop_rmbase)
%   6  Down-sample to 1 kHz             (pop_resample)
%   7  Save dataset  →  <…>-step2_basic_processing.set
%   *  Visualise evoked potentials after steps 5 & 6 and save figures
%==========================================================================

clear; close all; clc;

PIPELINE_NAME = 'TESA_ICA';
STEP_NAME     = 'step2_basic_processing';

activity_prob = 5;

% -------------------------------------------------------------------------
% >>> EDIT THIS IF YOU PROCESS A DIFFERENT SESSION <<<
session = 'Pos10_80';          % e.g. 'Pos05_60', 'Baseline', …
% -------------------------------------------------------------------------

%% 0 ─ Environment ---------------------------------------------------------
ENV_PATH = fullfile('..','.env');
if isfile(ENV_PATH), loadenv(ENV_PATH); end

ROOT_DIR        = getenv('TMS_EEG_ROOT_DIR');
EXPERIMENT_NAME = getenv('EXPERIMENT_NAME');
PARTICIPANT_ID  = getenv('PARTICIPANT_ID');
assert(isfolder(ROOT_DIR),'ROOT_DIR not found');

%% 1 ─ Load step-1 dataset -------------------------------------------------
in_dir   = fullfile(ROOT_DIR,EXPERIMENT_NAME,PARTICIPANT_ID,'output',PIPELINE_NAME);
step1set = sprintf('%s-%s-%s-step1_create_set.set', ...
                   EXPERIMENT_NAME,PARTICIPANT_ID,session);
assert(isfile(fullfile(in_dir,step1set)), ...
      'Step-1 file "%s" not found in %s',step1set,in_dir);

[ALLEEG,EEG,CURRENTSET,ALLCOM] = eeglab;
EEG = pop_loadset('filename',step1set,'filepath',in_dir);
EEG = eeg_checkset(EEG);

fprintf('Loaded  %s  (%.0f Hz)\n',EEG.setname,EEG.srate);

%% 2 ─ Detect TMS pulses ---------------------------------------------------
EEG = pop_tesa_findpulse(EEG,'Cz', ...
                         'refract',4, ...
                         'rate',1e4, ...
                         'tmsLabel','TMS', ...
                         'plots','on');
EEG = eeg_checkset(EEG);

%% 3 ─ Reject noisy channels (“band electrodes”) --------------------------
origLabels = {EEG.chanlocs.labels};
[EEG,badIdx] = pop_rejchan(EEG, 'elec', 1:size(EEG.data,1) ,'threshold',activity_prob,'norm','on','measure','kurt');
if ~isempty(badIdx)
    fprintf('Rejected channels: %s\n', strjoin(origLabels(badIdx),', '));
end
EEG = eeg_checkset(EEG);

%% 4 ─ Epoch  –1 … +1 s around “TMS” events -------------------------------
EEG = pop_epoch(EEG,{'TMS'},[-1 1],'epochinfo','yes');
EEG = eeg_checkset(EEG);

%% 5 ─ Remove –2 … +5 ms artifact & interpolate ---------------------------
EEG = pop_tesa_removedata(EEG,[-2 5]);
EEG = pop_tesa_interpdata(EEG,'cubic',[1 1]);
EEG = eeg_checkset(EEG);

% ----- Visualise after pulse interpolation --------------------------------
visualize_eeg_evoked(EEG, [-750 0], [-100 400], [-50 50], 'After TMS pulse interpolate');
figName = sprintf('%s-%s-%s-%s-after_tms_pulse_interpolate.png', ...
                  EXPERIMENT_NAME, PARTICIPANT_ID, session, STEP_NAME);

out_dir = fullfile(ROOT_DIR,EXPERIMENT_NAME,PARTICIPANT_ID,'output',PIPELINE_NAME);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
saveas(gcf, fullfile(out_dir, figName));
close(gcf);

%% 6 ─ Demean each epoch ---------------------------------------------------
EEG = pop_rmbase(EEG,[]);
EEG = eeg_checkset(EEG);

% ----- Visualise after demean --------------------------------------------
visualize_eeg_evoked(EEG, [-750 0], [-100 400], [-50 50], 'After demean');
figName = sprintf('%s-%s-%s-%s-after_demean.png', ...
                  EXPERIMENT_NAME, PARTICIPANT_ID, session, STEP_NAME);

saveas(gcf, fullfile(out_dir, figName));
close(gcf);

%% 7 ─ Down-sample to 1 kHz ------------------------------------------------
if EEG.srate ~= 1000
    EEG = pop_resample(EEG,1000);
end





%% 8 ─ Save dataset --------------------------------------------------------
save_fname = sprintf('%s-%s-%s-%s.set', ...
                     EXPERIMENT_NAME,PARTICIPANT_ID,session,STEP_NAME);

EEG.setname  = save_fname;
EEG.filename = save_fname;
EEG.filepath = out_dir;
pop_saveset(EEG,'filename',save_fname,'filepath',out_dir);

fprintf('\nStep-2 complete  ➜  %s\n',fullfile(out_dir,save_fname));
