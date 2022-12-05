%% detection and visual scoring of CCEPs
% author: Dorien van Blooijs, Sifra Blok

% 1. select patient
% 2. load ECoG data, 
% 3. split into stimulation trials, 
% 4. merge files when multiple runs
% 5. re-reference the data,
% 6. only include stimulus pairs that are stimulated in both SPES clin and
%   SPES-propofol
% 7. detect N1-peaks after stimulation, 
% 8. save these for further analysis and display ->
% sub-PRIOSXX_ses_X_task-SPESXXXX_CCEP.mat

close all 
clear; 
clc;

%% set paths

% add current path from folder which contains this script
rootPath = matlab.desktop.editor.getActiveFilename;
RepoPath = fileparts(rootPath);
matlabFolder = strfind(RepoPath,'matlab');
addpath(genpath(RepoPath(1:matlabFolder+6)));

myDataPath = PRIOS_setLocalDataPath(1);

% housekeeping 
clear rootPath RepoPath matlabFolder

%% Select patient

cfg.sub_labels = {['sub-' input('Patient number (PRIOSXX): ','s')]};
cfg.ses_label = 'ses-*';
cfg.task_label = 'task-SPES*';

%% Load ECOG data

% Find files with SPESclin and SPESprop 
files = dir(fullfile(myDataPath.dataPath,cfg.sub_labels{1}, cfg.ses_label,'ieeg',...
    [cfg.sub_labels{1} '_' cfg.ses_label '_task-SPES*_run-*_events.tsv']));
fileNames = {files.name};

% Find run_labels
% pre-allocation
run_labels = cell(size(fileNames));

for iRun = 1:size(fileNames,2)
    run_labels{iRun} = fileNames{iRun}(strfind(fileNames{iRun},'run-'):strfind(fileNames{iRun},'run-')+9);
end

% Load all data (both SPESclin and SPESprop)
cfg.run_label = run_labels;
dataBase = load_ECoGdata(cfg, myDataPath);

fprintf('Data of subject %s is loaded. \n',cfg.sub_labels{1})

% housekeeping
clear iRun files run_labels fileNames
    
%% prepocess into epochs time-locked to stimulus
% pre-process the raw data into epochs that are time-locked to the stimulus
% artefact. Two databases are created: one for the propofol SPES
% (dataBase_prop) and one for the clinical SPES (dataBase_clin). 
% These databases are compared in later steps. 

% configurations of epochs: epoch length
cfg.epoch_length = 4; % in seconds
cfg.epoch_prestim = 2; % in seconds, with 4 seconds total resulting in -2:2

clear dataBase_clin dataBase_prop
countClin = 1; countProp = 1;

% save all stimuli of clinical SPES or propofol SPES
for iRun = 1:size(dataBase,2)

    if strcmpi(dataBase(iRun).task_label,'task-SPESclin')

       cfg.minstim = 5; % number of events needed in one stimulation pair to be used in further analysis
       
       dataBase_clin(countClin) = preprocess_ECoG_spes(dataBase(iRun),cfg);  

        fprintf('...%s, %s, %s, has been epoched and averaged... \n', ...
            dataBase(iRun).sub_label, dataBase(iRun).run_label, ...
            dataBase(iRun).task_label)
        
        countClin = countClin +1;

    elseif strcmpi(dataBase(iRun).task_label,'task-SPESprop')

        cfg.minstim = 1; % number of events needed in one stimulation pair to be used in further analysis
        
        dataBase_prop(countProp) = preprocess_ECoG_spes(dataBase(iRun),cfg);   

        fprintf('...%s, %s, %s, has been epoched and averaged... \n', ...
            dataBase(iRun).sub_label, dataBase(iRun).run_label, ...
            dataBase(iRun).task_label)

        countProp = countProp + 1;
    end
end

disp('Completed preprocessing the data')

% housekeeping
clear countProp countClin iRun

%% When a SPES is run in multiple runs, merge them.

if size(dataBase_clin,2)> 1    

%     dataBase_clinOrig = dataBase_clin; % enables check of the original databases
    dataBase_clin = merge_runs(dataBase_clin); 
    
elseif size(dataBase_prop,2) > 1 

%     dataBase_propOrig = dataBase_prop; % enables check of the original databases
    dataBase_prop = merge_runs(dataBase_prop);
     
end

%% Re-reference data
% The reference calculatede below is the median of the signals that have a
% post- and pre-stim variance that is smaller than the pre-stim variance of
% the median signal of all signals (CAR). At least 5% of the signals have
% to be averaged in the reference signal. 

% Lowest variance is used to avoid introducing CCEP's in signals that
% originally do not have a CCEP.

cfg.reref = 'y';

tic
dataBase_clin = rereference_with_lowest_var(dataBase_clin);
dataBase_prop = rereference_with_lowest_var(dataBase_prop);
toc

disp('Data is rereferenced and all trials of each stimulus pair are averaged.')

%% Check whether SPESclin and SPESprop contain the same stimulation pairs
% Stimpairs and electrodes which are different in the clinical and propofol
% SPES are removed.

[dataBase_clin, dataBase_prop] = similar_stimpairs(dataBase_clin, dataBase_prop);

%% Use the automatic N1 detector to detect ccep 

% parameters used in detector
cfg.amplitude_thresh = 2.6;
cfg.n1_peak_range = 100;

% detect CCEPs
dataBase_clin = detect_n1peak_ECoG_ccep(dataBase_clin,cfg);
dataBase_prop = detect_n1peak_ECoG_ccep(dataBase_prop,cfg);

disp('Detection of CCEPs is completed')

%% Save CC epoch sorted for later usage
% Saving the CCEP file takes approx 5 minutes

filename_clin = [dataBase_clin.sub_label,'_',dataBase_clin.ses_label,'_',dataBase_clin.task_label,'_CCEP.mat'];
filename_prop = [dataBase_prop.sub_label,'_',dataBase_prop.ses_label,'_',dataBase_prop.task_label,'_CCEP.mat'];

filefolder = fullfile(myDataPath.dataPath, 'derivatives', 'CCEPs');
if ~exist(filefolder,'dir')
    mkdir(filefolder)
end

% remove the unnecessary fields in the dataBase structs
dataBase_prop = rmfield(dataBase_prop,{'tb_events','tb_channels', ...
    'tb_electrodes','data','Burstsup','Seizure','cc_epoch_sorted', ...
    'cc_epoch_sorted_avg','ref'});

dataBase_clin = rmfield(dataBase_clin,{'tb_events','tb_channels', ...
    'tb_electrodes','data','Burstsup','Seizure','cc_epoch_sorted', ...
    'cc_epoch_sorted_avg','ref'});

% Save CCEPs for later analysis
save(fullfile(filefolder,filename_clin),'-struct','dataBase_clin','-v7.3');
save(fullfile(filefolder,filename_prop),'-struct','dataBase_prop','-v7.3');


%% END OF SCRIPT