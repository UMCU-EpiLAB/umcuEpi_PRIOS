%% config_CCEP

%% patient characteristics
cfg.sub_labels = {['sub-' input('Patient number (PRIOSXX): ','s')]};
cfg.ses_label = input('Session number (ses-X): ','s');
cfg.task_label = 'task-SPES*';
%cfg.run_label = input('Runlabel(run-XXXXXX): ','s');

%% pre-allocation
% preprocessing step
cfg.dir = 'yes'; 
% yes, if you want to take negative/positive stimulation into account, 
% no if you want to use all stimuli in C1-C2 and C2-C1 and make no difference between those two
cfg.amp = 'no';  
% yes, if you want to take stimulation current into account, 
% no if you want to use all stimuli in C1-C2 4mA and C1-C2 8 mA and make no difference between those two


% epoch length
cfg.epoch_length = 4; % in seconds
cfg.epoch_prestim = 2; % in seconds, with 4 seconds total resulting in -2:2

% number of events needed in one stimulation pair to be used in further
% analysis
cfg.minstim = 5;

% detect CCEPs
cfg.amplitude_thresh = 2.6;
cfg.n1_peak_range = 100;

