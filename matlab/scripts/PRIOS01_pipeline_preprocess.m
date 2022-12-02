%% detection and visual scoring of CCEPs
% author: Dorien van Blooijs, Sifra Blok

% 1. load ECoG data, 
% 2. split into stimulation trials, 
% 3. re-reference the data,
% 4. detect N1-peaks after stimulation, 
% 5. visual check of detected N1s, 
% 6. save these for further analysis and display

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
    
%% CCEP for Single Pulse Electrical Stimulations
% pre-process the raw data into epochs that are time-locked to the stimulus
% artefact. Two databases are created: one for the propofol SPES
% (dataBase_prop) and one for the clinical SPES (dataBase_clin). 
% These databases are compared in later steps. 

% configurations of epochs
cfg.dir = 'yes'; 
% yes, if you want to take negative/positive stimulation into account, 
% no if you want to use all stimuli in C1-C2 and C2-C1 and make no difference between those two
cfg.amp = 'no';  
% yes, if you want to take stimulation current into account, 
% no if you want to use all stimuli in C1-C2 4mA and C1-C2 8 mA and make no difference between those two
% epoch length
cfg.epoch_length = 4; % in seconds
cfg.epoch_prestim = 2; % in seconds, with 4 seconds total resulting in -2:2

clear dataBase_clin dataBase_prop
countClin = 1; countProp = 1;

% FIXTHIS: preprocess_ECoG_spes moet nog flink uitgedund worden en
% herschreven! 
% cc_epoch_sorted = [80x5x98x8192]

% cc_epoch_sorted_select = [80x49x10x8192] --> deze wordt gerereferenced
% denk ik --> cc_epoch_sorted_select_reref, en dan wordt hier de average
% van genomen --> cc_epoch_sorted_select_reref_avg, en deze worden in
% select_n1_latency gebruikt om te plotten.

% cc_epoch_sorted_avg = [80x49x8192]

% --> de vraag is: waar wordt verder mee gewerkt?

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
%FIXTHIS: ik heb in merge_runs al een description toegevoegd, maar deze
%moet ik nog updaten als ik preprocess_ECoG_spes heb opgeschoond. 

if size(dataBase_clin,2)> 1    

    dataBase_clinOrig = dataBase_clin; % enables check of the original databases
    dataBase_clin = merge_runs(dataBase_clin); 
    
elseif size(dataBase_prop,2) > 1 

    dataBase_propOrig = dataBase_prop;
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

% FIXTHIS: ik moet deze functie nog opschonen!
tic
dataBase_clin = rereference_with_lowest_var(dataBase_clin, cfg);
dataBase_prop = rereference_with_lowest_var(dataBase_prop, cfg);
toc

disp('Data is rereferenced and all trials of each stimulus pair are averaged.')

%% Check whether SPESclin and SPESprop contain the same stimulation pairs
% Stimpairs and electrodes which are different in the clinical and propofol
% SPES are removed.

% FIXTHIS: ik moet deze functie nog opschonen. 
[dataBase_clin, dataBase_prop] = similar_stimpairs(dataBase_clin, dataBase_prop);

%% Use the automatic N1 detector to detect ccep 

% FIXTHIS: ik moet deze functie nog opschonen.

% detect CCEPs
cfg.amplitude_thresh = 2.6;
cfg.n1_peak_range = 100;

dataBase_clin = detect_n1peak_ECoG_ccep(dataBase_clin,cfg);
dataBase_prop = detect_n1peak_ECoG_ccep(dataBase_prop,cfg);

disp('Detection of CCEPs is completed')

%% Save CC epoch sorted for later usage
% Saving the CCEP file takes approx 5 minutes

filename_clin = [dataBase_clin.sub_label,'_',dataBase_clin.ses_label,'_',dataBase_clin.task_label,'_CCEP.mat'];
filename_prop = [dataBase_prop.sub_label,'_',dataBase_clin.ses_label,'_',dataBase_prop.task_label,'_CCEP.mat'];

filefolder = fullfile(myDataPath.dataPath, 'derivatives', 'CCEPs');
if ~exist(filefolder,'dir')
    mkdir(filefolder)
end

% Save CCEPs for later analysis
save(fullfile(filefolder,filename_clin),'-struct','dataBase_clin','-v7.3');
save(fullfile(filefolder,filename_prop),'-struct','dataBase_prop','-v7.3');

%% Visually check detected cceps

% FIXTHIS: ik moet deze functie nog opschonen. 

% Check the average signal in which a CCEP was detected

% load checked N1s if visual rating has started earlier
if exist(fullfile(myDataPath.CCEPpath, dataBase_clin.sub_label, ...
        dataBase_clin.ses_label, dataBase_clin.task_label,...
        [dataBase_clin.sub_label, '_', dataBase_clin.ses_label,'_',...
        dataBase_clin.task_label,'_N1sChecked.mat']),'file')

    load(fullfile(myDataPath.CCEPpath, dataBase_clin.sub_label, ...
        dataBase_clin.ses_label, dataBase_clin.task_label,...
        [dataBase_clin.sub_label, '_', dataBase_clin.ses_label,'_',...
        dataBase_clin.task_label,'_N1sChecked.mat']));
   dataBase_clin.ccep = ccep; 
end
   
% continue with the stimulation pair after the last saved stimulation pair
if sum(strcmp(fieldnames(dataBase_clin.ccep), 'checkUntilStimp')) == 1
    endstimp = ccep.checkUntilStimp;
else
    endstimp = 0;
end

visualRating_ccep(dataBase_clin, cfg, endstimp, myDataPath);

disp('visual rating of CCEPs from SPES clinical is completed')      

%% Visual check of propofol-SPES

if exist(fullfile(myDataPath.CCEPpath, dataBase_prop.sub_label, ...
        dataBase_prop.ses_label, dataBase_prop.task_name,...
        [dataBase_prop.sub_label, '_', dataBase_prop.ses_label,'_',...
        dataBase_prop.task_name,'_N1sChecked.mat']),'file')
    
   load(fullfile(myDataPath.CCEPpath, dataBase_prop.sub_label, ...
        dataBase_prop.ses_label, dataBase_prop.task_name,...
        [dataBase_prop.sub_label, '_', dataBase_prop.ses_label,'_',...
        dataBase_prop.task_name,'_N1sChecked.mat']));  
   dataBase_prop.ccep = ccep;
end

% continue with the stimulation pair after the last saved stimulation pair
if sum(strcmp(fieldnames(dataBase_prop.ccep), 'checkUntilStimp')) == 1
    endstimp_prop = dataBase_prop.ccep.checkUntilStimp;
else
    endstimp_prop = 0;
end

visualRating_ccep(dataBase_prop, cfg, endstimp_prop, myDataPath);

disp('visual rating of CCEPs from SPES propofol is completed')      

%% END OF SCRIPT