%% visual check of detected N1s

% 1. set paths
% 2. select patient
% 3. load sub-PRIOSXX_ses-X_task-SPESXXXX_CCEP.mat
% 4. visually check detected CCEPs
% 5. save as sub-PRIOSXX_ses-X_task-SPESXXXX_N1sChecked_XX.mat

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
cfg.ses_label = 'ses-1';
cfg.task_label = ['task-' input('Task [SPESclin/SPESprop]: ','s')];

%% load data

fileFolder = fullfile(myDataPath.dataPath,'derivatives','CCEPs');
fileName = [cfg.sub_labels{1}, '_', cfg.ses_label,'_', cfg.task_label, '_CCEP.mat'];

if exist(fullfile(fileFolder, ...
        fileName),'file')
    tic
    dataBase = load(fullfile(fileFolder, fileName));
    toc

else

    error('%s was not found. Please run first PRIOS01_pipeline_preprocess.m',...
        fullfile(fileFolder,fileName))

end

%% Visually check detected cceps
% Check the average signal in which a CCEP was detected

cfg.observer = input('Observer [SB/DvB]: ','s');
fileNameN1s = [dataBase.sub_label, '_', dataBase.ses_label,'_',...
        dataBase.task_label,'_N1sChecked_', cfg.observer, '.mat'];

% load checked N1s if visual rating has started earlier
if exist(fullfile(myDataPath.CCEPpath, 'checkedN1s',...
        fileNameN1s),'file')

    dataBase.ccep = load(fullfile(myDataPath.CCEPpath, 'checkedN1s', ...
        fileNameN1s));

end
   
% continue with the stimulation pair after the last saved stimulation pair
if sum(strcmp(fieldnames(dataBase.ccep), 'checkUntilStimp')) == 1
    endstimp = dataBase.ccep.checkUntilStimp;
else
    endstimp = 0;
end

visualRating_ccep(dataBase, cfg, endstimp, myDataPath);

fprintf('Visual rating of %s is completed \n',fileNameN1s)      

%% end script