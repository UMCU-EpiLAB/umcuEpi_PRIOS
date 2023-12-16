% compare number of CCEPs to each stimulus pair in propofol and
% clinical SPES

clear
close all
clc

%% general variables

cmap = parula(8);
pFDR = 0.05; % for FDR correction

%% set paths

% add current path from folder which contains this script
rootPath = matlab.desktop.editor.getActiveFilename;
RepoPath = fileparts(rootPath);
matlabFolder = strfind(RepoPath,'matlab');
addpath(genpath(RepoPath(1:matlabFolder+6)));

myDataPath = PRIOS_setLocalDataPath(1);

% housekeeping 
clear rootPath RepoPath matlabFolder

%% select subjects

all_sublabels = {'sub-PRIOS01','sub-PRIOS02','sub-PRIOS03','sub-PRIOS04',...
    'sub-PRIOS05','sub-PRIOS09'};

%% load visually checked N1s

clear dataBase

for nSub = 1:size(all_sublabels,2)

    % load visually checked N1s of SPES-clinical
    if ~exist(fullfile(myDataPath.CCEPpath,'checkedN1s',...
            [all_sublabels{nSub},'_ses-1_task-SPESclin_N1sChecked_comb.mat']),'file')

        error('You should first run PRIOS01_pipeline_preprocess.m by two observers, PRIOS02_visualCheckN1s.m and PRIOS03_pipeline_agreement.m before you can run this script.')

    else

        tmp = load(fullfile(myDataPath.CCEPpath,'checkedN1s',...
            [all_sublabels{nSub},'_ses-1_task-SPESclin_N1sChecked_comb.mat']));
       
        dataBase(nSub).sub_label = all_sublabels{nSub};
        dataBase(nSub).cc_stimchans = tmp.cc_stimchans;
        dataBase(nSub).cc_stimsets = tmp.cc_stimsets;
        dataBase(nSub).spesClin  = tmp;

    end

    % load visually checked N1s of SPES-propofol
    if ~exist(fullfile(myDataPath.CCEPpath,'checkedN1s',...
            [all_sublabels{nSub},'_ses-1_task-SPESprop_N1sChecked_comb.mat']),'file')

        error('You should first run PRIOS01_pipeline_preprocess.m by two observers, PRIOS02_visualCheckN1s.m and PRIOS03_pipeline_agreement.m before you can run this script.')

    else

        tmp = load(fullfile(myDataPath.CCEPpath,'checkedN1s',...
            [all_sublabels{nSub},'_ses-1_task-SPESprop_N1sChecked_comb.mat']));

        dataBase(nSub).sub_label = all_sublabels{nSub};
        dataBase(nSub).cc_stimchans = tmp.cc_stimchans;
        dataBase(nSub).cc_stimsets = tmp.cc_stimsets;
        dataBase(nSub).spesProp = tmp;
    end
end

% housekeeping
clear nSub all_sublabels tmp

%% display number of channels and number of stimulus pairs
clc

for nSub = 1:size(dataBase,2)
    
    fprintf('Subject: %s, number of channels (clin/prop): %2.0f/%2.0f, number of stimulus pairs: %2.0f/%2.0f \n',...
        dataBase(nSub).sub_label, size(dataBase(nSub).spesClin.ch,1), size(dataBase(nSub).spesProp.ch,1),...
        size(dataBase(nSub).spesClin.cc_stimsets,1),size(dataBase(nSub).spesProp.cc_stimsets,1))
end