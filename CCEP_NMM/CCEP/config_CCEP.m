%% config_CCEP

addpath(genpath('git_rep/BasicCode_ECoG_DvB'))
addpath(genpath('git_rep/eeglab/'))     
addpath('git_rep/fieldtrip/')
ft_defaults

%% patient characteristics
cfg.sub_labels = {['sub-' input('Patient number (RESPXXXX or REC2StimXX): ','s')]};
cfg.ses_label = input('Session number (ses-X): ','s');
cfg.task_label = 'task-SPESclin';
cfg.run_label = {['run-' input('Run [daydayhhminmin]: ','s')]};

%% pre-allocation
% preprocessing step
cfg.dir = 'no'; % if you want to take negative/positive stimulation into account
cfg.amp = 'no'; % if you want to take stimulation current into account

% epoch length
cfg.epoch_length = 4; % in seconds
cfg.epoch_prestim = 2; % in seconds, with 4 seconds total resulting in -2:2

% detect CCEPs
cfg.amplitude_thresh = 2.6;
cfg.n1_peak_range = 100;

%% set path
if contains(cfg.sub_labels,'REC2Stim')
    % REC2Stim study
    cfg.CCEPpath = '/Fridge/REC2Stimstudy/derivatives/CCEP/';
    cfg.dataPath = '/Fridge/REC2Stimstudy/';
elseif contains(cfg.sub_labels,'RESP')
    % RESPect database
    cfg.CCEPpath = '/Fridge/users/dorien/derivatives/CCEP/';
    cfg.dataPath = '/Fridge/chronic_ECoG/';
end