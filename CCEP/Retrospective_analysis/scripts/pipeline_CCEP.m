%% pipeline CCEP
% author: Dorien van Blooijs
% date: September 2019

%% pre-allocation
clear 

config_CCEP

%% set paths
cfg.mode = 'retro';
myDataPath = setLocalDataPath(cfg);

%% select run

% choose between available runs
files = dir(fullfile(myDataPath.dataPath,cfg.sub_labels{1}, cfg.ses_label,'ieeg',...
    [cfg.sub_labels{1} '_' cfg.ses_label '_' cfg.task_label '_*'  '_events.tsv']));
names = {files.name};
strings = cell(size(names));
for i=1:size(names,2)
    strings{i} = names{i}(strfind(names{i},'run-'):strfind(names{i},'run-')+9);
end
stringsz = [repmat('%s, ',1,size(strings,2)-1),'%s'];

cfg.run_label = {input(sprintf(['Choose one of these runs: \n' stringsz '\n'],strings{:}),'s')}; % Chosen run is in cfg.run_label

clear files names strings stringsz

%% load data

dataBase = load_ECoGdata(cfg,myDataPath);

%% preprocess ccep
% eventnum(1)
dataBase = preprocess_ECoG_spes(dataBase,cfg);

%% detect ccep

[dataBase] = detect_n1peak_ECoG_ccep(dataBase, cfg);
dataBase.ccep.amplitude_thresh = cfg.amplitude_thresh;
dataBase.ccep.n1_peak_range = cfg.n1_peak_range;
dataBase.ccep.cc_stimsets = dataBase.cc_stimsets;
dataBase.ccep.ch = dataBase.ch;
dataBase.ccep.stimpnames = dataBase.stimpnames;
dataBase.ccep.stimchans = dataBase.cc_stimchans;

disp('Detection of ERs is completed')

%% visualize CCEPs per electrode
cfg.save_fig = str2double(input('Do you want to save the figures? [yes = 1, no = 0]: ','s'));
plot_ccep_av(dataBase,cfg);

%% visually check detected ccepsyy

dataBase = visualRating_ccep(dataBase);
% correct: y
% incorrect: n or enter

%% save ccep

targetFolder = [myDataPath.CCEPpath, dataBase(1).sub_label,'/',dataBase(1).ses_label,'/', dataBase(1).run_label,'/'];

% Create the folder if it doesn't exist already.
if ~exist(targetFolder, 'dir')
    mkdir(targetFolder);
end

start_filename = strfind(dataBase(1).dataName,'/');
stop_filename = strfind(dataBase(1).dataName,'_ieeg');

fileName=[dataBase(1).dataName(start_filename(end)+1:stop_filename-1),'_CCEP_2stims.mat'];

ccep = dataBase.ccep;
ccep.dataName = dataBase(1).dataName;

save([targetFolder,fileName], 'ccep');

fprintf('CCEP is saved in %s%s \n',targetFolder,fileName)