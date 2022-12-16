% function to read BIDS-data into dataBase
% author: Dorien van Blooijs
% date: June 2019

% INPUTS:
% cfg: 
%   is a struct with the following fields:
% - sub_labels
%   cell{1 x subjects} containing the names of the subjects
% - ses_label
%   string containing the name of the session
% - task_label
%   string containing the name of the task
% - run_label
%   cell{1 x runs} containing the names of the runs
% 
% myDataPath: 
%   is a struct with the following fields:
% - dataPath
%   string containing the folder path where the BIDS files are located

% OUTPUTS:
% dataBase: 
%   is a struct that is created in this function. To the database structure
%   the following matrices will be added:
% - sub_label
%   string containing the name of the subject
% - ses_label
%   string containing the name of the session
% - task_label
%   string containg the name of the task
% - run_label
%   string containing the name of the run(s)
% - dataName
%   name of the eeg-file from which data is extracted
% - ccep_header
%   struct containing, among other things, the field: Fs (sample frequency)
% - tb_events
%   table containing information regarding events in the eeg-file (see BIDS
%   structure)
% - tb_channels
%   table containing information regarding the recording channels (see BIDS
%   structure)
% - tb_electrodes
%   table containing information regarding the hardware electrodes (see
%   BIDS structure)
% - ch
%   cell[channels x 1] containing the channel names 
% - data
%   matrix[channels x samples] containg the eeg-signals 

function dataBase = load_ECoGdata(cfg,myDataPath)

dataPath = myDataPath.dataPath;
sub_labels = cfg.sub_labels;
ses_label = cfg.ses_label;
run_label = cfg.run_label;

dataBase = struct([]);
for iRun = 1:size(run_label,2)
    D = dir(fullfile(dataPath,[sub_labels{1}],ses_label,'ieeg',...
        [sub_labels{1} '_' ses_label '_' cfg.task_label ,'_',run_label{iRun}, '_ieeg.eeg']));

    if size(D,1) == 0
        error('%s does not exist',fullfile(dataPath,sub_labels{:},ses_label,'ieeg',...
            [sub_labels{:} '_' ses_label '_' cfg.task_label ,'_',run_label{iRun}, '_ieeg.eeg']))
    end

    % determine run_label
    dataName = fullfile(D(1).folder, D(1).name);

    % use real ses label instead of ses-*
    ses_temp = extractBetween(dataName,'ses-','/');
    ses_label = ['ses-',ses_temp{1}];

    % use real task label instead of task-SPES*
    task_temp = extractBetween(dataName,'task-','_');
    task_label = ['task-',task_temp{1}];

    % load data
    ccep_data = ft_read_data(dataName,'dataformat','brainvision_eeg');
    ccep_header = ft_read_header(dataName);

    % load events
    D = dir(fullfile(dataPath,[sub_labels{1}],ses_label,'ieeg',...
        [sub_labels{1} '_' ses_label '_' task_label ,'_',run_label{iRun},'_events.tsv']));

    eventsName = fullfile(D(1).folder, D(1).name);

    tb_events = readtable(eventsName,'FileType','text','Delimiter','\t');

    % load electrodes
    D = dir(fullfile(dataPath,[sub_labels{1}],ses_label,'ieeg',...
        [sub_labels{1} '_' ses_label ,'_electrodes.tsv']));

    elecsName = fullfile(D(1).folder, D(1).name);

    tb_electrodes = readtable(elecsName,'FileType','text','Delimiter','\t');
    idx_elec_incl = ~strcmp(tb_electrodes.group,'other');
    tb_electrodes = tb_electrodes(idx_elec_incl,:);                         % remove all the electrodes with comment other in column group

    % load channels
    D = dir(fullfile(dataPath,[sub_labels{1}], ses_label,'ieeg',...
        [sub_labels{1} '_' ses_label '_' task_label '_' run_label{iRun} '_channels.tsv']));

    channelsName = fullfile(D(1).folder, D(1).name);

    tb_channels = readtable(channelsName,'FileType','text','Delimiter','\t');
    idx_ch_incl = strcmp(tb_channels.type,'ECOG')|strcmp(tb_channels.type,'SEEG'); % remove all the electrodes without comment ECOG or SEEG
% FIXTHIS: ook depth moet ik eruit halen!

    tb_channels = tb_channels(idx_ch_incl,:);
    ch_incl = tb_channels.name;

    data = ccep_data(idx_ch_incl,:);

    dataBase(iRun).sub_label = sub_labels{1};
    dataBase(iRun).ses_label = ses_label;
    dataBase(iRun).task_label = task_label;
    dataBase(iRun).run_label = run_label{iRun};
    dataBase(iRun).dataName = dataName;
    dataBase(iRun).ccep_header = ccep_header;
    dataBase(iRun).tb_events = tb_events;
    dataBase(iRun).tb_channels = tb_channels;
    dataBase(iRun).tb_electrodes = tb_electrodes;
    dataBase(iRun).ch = ch_incl;
    dataBase(iRun).data = data;

end
end