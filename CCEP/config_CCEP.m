%% config_CCEP

%% patient characteristics
cfg.sub_labels = {['sub-' input('Patient number (RESPXXXX: ','s')]};
cfg.ses_label = input('Session number (ses-X): ','s');
cfg.task_label = 'task-SPES*';
%cfg.run_label = input('Runlabel(run-XXXXXX): ','s');
%cfg.run_label = {[input('Runlabel (run-XXXXXX): ','s')]};


% files = dir(fullfile(myDataPath.dataPath,cfg.sub_labels{1}, cfg.ses_label,'ieeg',...
%     [cfg.sub_labels{1} '_' cfg.ses_label '_' cfg.task_label '_*'  '_events.tsv']));
% names = {files.name};
% strings = cellfun(@(x) x(strfind(names{1},'run-'):strfind(names{1},'run-')+9), names, 'UniformOutput', false);
% stringsz = [repmat('%s, ',1,size(strings,2)-1),'%s'];
% 
% cfg.run_label = {input(sprintf(['Choose one of these runs: \n' stringsz '\n'],strings{:}),'s')}; % Chosen run is in cfg.run_label


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

