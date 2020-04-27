%% config_CCEP

%% patient characteristics
cfg.sub_labels = {['sub-' input('Patient number (RESPXXXX): ','s')]};
cfg.ses_label = input('Session number (ses-X): ','s');
cfg.task_label = 'task-SPES*';
%cfg.run_label = input('Runlabel(run-XXXXXX): ','s');

%% pre-allocation
% preprocessing step
cfg.dir = 'yes'; % if you want to take negative/positive stimulation into account
cfg.amp = 'no';  % if you want to take stimulation current into account

% epoch length
cfg.epoch_length = 4; % in seconds
cfg.epoch_prestim = 2; % in seconds, with 4 seconds total resulting in -2:2

% detect CCEPs
cfg.amplitude_thresh = 2.6;
cfg.n1_peak_range = 100;

