%% config_CCEP


%% pre-allocation
% preprocessing step
cfg.dir = 'yes'; 
% yes, if you want to take negative/positive stimulation into account, 
% no if you want to use all stimuli in C1-C2 and C2-C1 and make no difference between those two

cfg.dir_avg = 'yes';
% yes, if you take negative/positive stimulation into account, but you do want
% to average some of both for detection of cceps
% no, if you use all stimuli in C1-C2 and C2-C1 separately, but want to
% average both separately for detection of cceps
% if cfg.dir = 'no', cfg.dir_avg is not used in the script

% epoch length
cfg.epoch_length = 4; % in seconds
cfg.epoch_prestim = 2; % in seconds, with 4 seconds total resulting in -2:2

% number of events needed in one stimulation pair to be used in further
% analysis
cfg.minstim = 5;

% detect CCEPs
cfg.amplitude_thresh = 2.6;
cfg.n1_peak_range = 100;

