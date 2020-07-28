clear; 

%% Choose patient
config_CCEP

%% set paths
% Adapt for RESP or PRIOS patients!
myDataPath = setLocalDataPath(cfg);

%% select run
% choose between available runs
files = dir(fullfile(myDataPath.dataPath,cfg.sub_labels{1}, cfg.ses_label,'ieeg',...
    [cfg.sub_labels{1} '_' cfg.ses_label '_' cfg.task_label '_*'  '_events.tsv']));
names = {files.name};
strings = cell(size(names));
for n = 1:size(names,2)
    strings{n} = names{n}(strfind(names{n},'run-'):strfind(names{n},'run-')+9);
end
stringsz = [repmat('%s, ',1,size(strings,2)-1),'%s'];

cfg.run_label = {input(sprintf(['Choose one of these runs: \n' stringsz '\n'],strings{:}),'s')}; % Chosen run is in cfg.run_label

if ~contains(cfg.run_label,'run-')
   error('"run-" is missing in run_label') 
end

clear files names strings stringsz

%% load data

dataBase = load_ECoGdata(cfg,myDataPath);

%% CCEP for 2 and 10 stimulations

% avg_stim : write down the number of stimuli you want to average

% save only first stimulus in both directions
avg_stim = 1;
dataBase2stim = preprocess_ECoG_spes(dataBase,cfg,avg_stim);

% save all stimuli (5) in each direction
avg_stim = 5;
dataBaseallstim = preprocess_ECoG_spes(dataBase,cfg,avg_stim);

% check whether similar stimuli are present in the same stimulus pair
chan = 13; stim=1;
figure, 
subplot(2,1,1),
plot(squeeze(dataBase2stim.cc_epoch_sorted_select_avg(chan,stim,:,:))','Color',[0.8 0.8 0.8],'LineWidth',1)
hold on
plot(squeeze(dataBase2stim.cc_epoch_sorted_avg(chan,stim,:)),'k','LineWidth',2)
hold off
title('two stimuli')

subplot(2,1,2),
plot(squeeze(dataBaseallstim.cc_epoch_sorted_select_avg(chan,stim,:,:))','Color',[0.8 0.8 0.8],'LineWidth',1)
hold on
plot(squeeze(dataBaseallstim.cc_epoch_sorted_avg(chan,stim,:)),'k','LineWidth',2)
hold off
title('all stimuli')

%% Visually check all averaged signals for ERs

% EDIT: hier moet je dit weg halen, of je moet het optioneel maken met iets van: 
% s = input('do you want to visually check ERs? [y/n] ','s')
% if strcmp(s,'y')
% dan wel visueel checken
% end

% dataBaseallstim = plot_all_signals(dataBaseallstim);
% dataBase2stim = plot_all_signals(dataBase2stim);
% 
% % Save 
% targetFolder = [myDataPath.CCEPpath, dataBase(1).sub_label,'/',dataBase(1).ses_label,'/', dataBase(1).run_label,'/'];
% 
% vis_checked_all = [dataBaseallstim.sub_label, '_',dataBase(1).run_label, '_avg_10stims_vis.mat'];
% ccep = dataBaseallstim.ccep;
% save([targetFolder,vis_checked_all], 'ccep')
% 
% vis_checked_2 = [dataBase2stim.sub_label, '_', dataBase(1).run_label, '_avg_2stims_vis.mat'];
% ccep = dataBase2stim.ccep;
% save([targetFolder,vis_checked_2], 'ccep')
% fprintf('Checked files are saved in %s \n',targetFolder);

%% detect ccep in only first stimulus in both directions and all stimuli
dataBase2stim = detect_n1peak_ECoG_ccep(dataBase2stim,cfg);
dataBaseallstim = detect_n1peak_ECoG_ccep(dataBaseallstim,cfg);
dataBase2stim.NmbrofStims = '2_stims';
dataBaseallstim.NmbrofStims = '10_stims';

disp('Detection of ERs is completed')

%% save ccep
savefiles = input('Do you want to save the ccep-structures? [y/n] ','s');

targetFolder = [fullfile(myDataPath.CCEPpath, dataBase(1).sub_label,dataBase(1).ses_label,dataBase(1).run_label),'/'];

% Create the folder if it doesn't exist already.
if ~exist(targetFolder, 'dir')
    mkdir(targetFolder);
end

[~,filename,~] = fileparts(dataBase(1).dataName);

% save 2 stims
fileName=[extractBefore(filename,'_ieeg'),'_CCEP_2stims.mat'];
ccep2 = dataBase2stim.ccep;
ccep2.stimchans = dataBase2stim.cc_stimchans_all;
ccep2.stimpnames = dataBase2stim.stimpnames_all;
ccep2.stimsets = dataBase2stim.cc_stimsets_all;
ccep2.dataName = dataBase2stim.dataName;
ccep2.ch = dataBase2stim.ch;

if strcmp(savefiles,'y')
    save([targetFolder,fileName], 'ccep2');
    save([myDataPath.CCEP_allpat,fileName], 'ccep2');
end

% save all stims
fileName5=[extractBefore(filename,'_ieeg'),'_CCEP_10stims.mat'];
ccep10 = dataBaseallstim.ccep;
ccep10.stimchans = dataBaseallstim.cc_stimchans_all;
ccep10.stimpnames = dataBaseallstim.stimpnames_all;
ccep10.stimsets = dataBaseallstim.cc_stimsets_all;
ccep10.dataName = dataBaseallstim.dataName;
ccep10.ch = dataBaseallstim.ch;

if strcmp(savefiles,'y')
    save([targetFolder,fileName5], 'ccep10');
    save([myDataPath.CCEP_allpat,fileName5], 'ccep10');

    fprintf('CCEPs 2stims and 10stims is saved in %s \n and %s',targetFolder , myDataPath.CCEP_allpat);
end
