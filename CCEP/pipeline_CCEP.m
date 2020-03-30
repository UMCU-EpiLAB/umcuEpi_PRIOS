%% pipeline CCEP
% author: Dorien van Blooijs
% date: September 2019

%% pre-allocation
clear 
config_CCEP

%% load data

dataBase = load_ECoGdata(cfg);

%% preprocess ccep

dataBase = preprocess_ECoG_spes(dataBase,cfg);

%% detect ccep

[dataBase] = detect_n1peak_ECoG_ccep(dataBase, cfg);
dataBase.ccep.amplitude_thresh = cfg.amplitude_thresh;
dataBase.ccep.n1_peak_range = cfg.n1_peak_range;
dataBase.ccep.cc_stimsets = dataBase.cc_stimsets;
dataBase.ccep.ch = dataBase.ch;

disp('Detection of ERs is completed')

%% visually check detected cceps

dataBase = visualRating_ccep(dataBase,cfg);

%% save ccep

targetFolder = [cfg.CCEPpath, dataBase(1).sub_label,'/',dataBase(1).ses_label,'/', dataBase(1).run_label,'/'];

% Create the folder if it doesn't exist already.
if ~exist(targetFolder, 'dir')
    mkdir(targetFolder);
end

start_filename = strfind(dataBase(1).dataName,'/');
stop_filename = strfind(dataBase(1).dataName,'_ieeg');

fileName=[dataBase(1).dataName(start_filename(end)+1:stop_filename-1),'_CCEP.mat'];

ccep = dataBase.ccep;
ccep.dataName = dataBase(1).dataName;

save([targetFolder,fileName], 'ccep');

fprintf('CCEP is saved in %s%s \n',targetFolder,fileName)
