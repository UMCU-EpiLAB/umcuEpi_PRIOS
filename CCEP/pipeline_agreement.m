clear; clc

%% Choose patient
config_CCEP

%% set paths
%myDataPath = setLocalDataPath(cfg);
myDataPath.CCEPpath = '/Fridge/users/sifra/derivatives/CCEP/' ;
myDataPath.dataPath = '/Fridge/chronic_ECoG/';

% set paths
addpath(genpath('/home/sifra/git_repositories/eeglab/')) ;    
addpath('/home/sifra/git_repositories/fieldtrip');
ft_defaults

localDataPath.CCEPpath = '/Fridge/users/sifra/derivatives/CCEP/'; % /Fridge/users/sifra/derivatives/CCEP
localDataPath.dataPath = '/Fridge/chronic_ECoG/';
%% select run

% choose between available runs
files = dir(fullfile(myDataPath.dataPath,cfg.sub_labels{1}, cfg.ses_label,'ieeg',...
    [cfg.sub_labels{1} '_' cfg.ses_label '_' cfg.task_label '_*'  '_events.tsv']));
names = {files.name};
strings = cellfun(@(x) x(strfind(names{1},'run-'):strfind(names{1},'run-')+9), names, 'UniformOutput', false);
stringsz = [repmat('%s, ',1,size(strings,2)-1),'%s'];

cfg.run_label = {input(sprintf(['Choose one of these runs: \n' stringsz '\n'],strings{:}),'s')}; % Chosen run is in cfg.run_label

clear files names strings stringsz

%% load data

dataBase = load_ECoGdata(cfg,myDataPath);

%% CCEP for 2 and 10 stimulations
stimulations = [2,5];
for K = 1:length(stimulations)
    [stim_dataBase(K)] = preprocess_ECoG_spes_test(dataBase,cfg,stimulations(K));

    % detect ccep
    [stim_database(K)] = detect_n1peak_ECoG_ccep(stim_dataBase(K), cfg);
    stim_database(K).ccep.amplitude_thresh = cfg.amplitude_thresh;
    stim_database(K).ccep.n1_peak_range = cfg.n1_peak_range;
    stim_database(K).ccep.cc_stimsets = stim_dataBase(K).cc_stimsets;
    stim_database(K).ccep.ch = stim_dataBase(K).ch;
    stim_database(K).ccep.stimpnames = stim_dataBase(K).stimpnames;
    stim_database(K).ccep.stimchans = stim_dataBase(K).cc_stimchans;  
end
 disp('Detection of ERs is completed')


%% save ccep
targetFolder = [myDataPath.CCEPpath, stim_database(1).sub_label,'/',stim_database(1).ses_label,'/', stim_database(1).run_label,'/'];

% Create the folder if it doesn't exist already.
if ~exist(targetFolder, 'dir')
    mkdir(targetFolder);
end

start_filename = strfind(stim_database(1).dataName,'/');
stop_filename = strfind(stim_database(1).dataName,'_ieeg');

    
if stim_database(1).stimnum == 2
    fileName=[stim_database(1).dataName(start_filename(end)+1:stop_filename-1),'_CCEP_2stims.mat'];

    ccep = stim_database.ccep;
    ccep.dataName = stim_database(1).dataName;

    save([targetFolder,fileName], 'ccep');
end

if stim_database(2).stimnum == 5
    fileName=[stim_database(1).dataName(start_filename(end)+1:stop_filename-1),'_CCEP_10stims.mat'];

    ccep = stim_database.ccep;
    ccep.dataName = stim_database(1).dataName;

    save([targetFolder,fileName], 'ccep');
end


fprintf('CCEPs is saved in %s%s \n',targetFolder,fileName)

%% detmine the agreement between 2 and 10 stims per run
clearvars -except localDataPath cfg
matrix_agreement = determine_agreement(localDataPath,cfg)
