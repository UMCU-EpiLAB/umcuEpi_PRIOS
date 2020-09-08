clear; 

%% Choose patient
config_CCEP

%% set paths
% Adapt for RESP or PRIOS patients!
cfg.mode = 'retro';
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
tt = dataBase2stim.tt;

% save all stimuli (5) in each direction
avg_stim = 5;
dataBaseallstim = preprocess_ECoG_spes(dataBase,cfg,avg_stim);

% check whether similar stimuli are present in the same stimulus pair
chan = 13; stim=1;
figure, 
subplot(2,1,1),
plot(tt,squeeze(dataBase2stim.cc_epoch_sorted_select_avg(chan,stim,:,:))','Color',[0.8 0.8 0.8],'LineWidth',1)
hold on
plot(tt,squeeze(dataBase2stim.cc_epoch_sorted_avg(chan,stim,:)),'k','LineWidth',2)
hold off
title('two stimuli')
xlabel('time (s)')
xlim([-.2 1.0])

            
subplot(2,1,2),
plot(tt,squeeze(dataBaseallstim.cc_epoch_sorted_select_avg(chan,stim,1:5,:))','Color','r','LineWidth',1)
hold on
plot(tt,squeeze(dataBaseallstim.cc_epoch_sorted_select_avg(chan,stim,6:10,:))','Color','b','LineWidth',1)
hold on
plot(tt,squeeze(dataBaseallstim.cc_epoch_sorted_avg(chan,stim,:)),'k','LineWidth',2)
hold off
title('all stimuli')
xlabel('time (s)')
xlim([-.2 1.0])


figure()
plot(tt,squeeze(dataBaseallstim.cc_epoch_sorted_select_avg(chan,stim,1:5,:))','Color','r','LineWidth',1)
hold on
plot(tt,squeeze(dataBaseallstim.cc_epoch_sorted_select_avg(chan,stim,6:10,:))','Color','b','LineWidth',1)
hold on
plot(tt,squeeze(dataBaseallstim.cc_epoch_sorted_avg(chan,stim,:)),'k','LineWidth',2)
hold off
title('all stimuli')
xlabel('time (s)')
xlim([-.1 0.1])


%% Use the automatic N1 detector to detect ccep 
dataBase2stim = detect_n1peak_ECoG_ccep(dataBase2stim,cfg);
dataBaseallstim = detect_n1peak_ECoG_ccep(dataBaseallstim,cfg);
dataBase2stim.NmbrofStims = '2_stims';
dataBaseallstim.NmbrofStims = '10_stims';

disp('Detection of ERs is completed')

%% Visually check detected cceps
VisCheck = input('Do you want to visually check the detected CCEPs? [y/n] ','s');

if strcmp(VisCheck,'y')
    dataBaseallstim = visualRating_ccep(dataBaseallstim);
    dataBase2stim = visualRating_ccep(dataBase2stim);

    % Save the values for the agreement per run (2 and 10 stims)
    targetFolder = [myDataPath.CCEPpath, dataBase(1).sub_label,'/',dataBase(1).ses_label,'/', dataBase(1).run_label,'/'];

    checked_all = [dataBaseallstim.sub_label, '_',dataBase(1).run_label, '_CCEP_', dataBaseallstim.NmbrofStims ,'_checked.mat'];
    save([targetFolder,checked_all], 'ccep')

    checked_2 = [dataBase2stim.sub_label, '_', dataBase(1).run_label, '_CCEP_', dataBase2stim.NmbrofStims ,'_checked.mat'];
    save([targetFolder,checked_2], 'ccep')
    fprintf('Checked files are saved in %s \n',targetFolder);

    % Perform the determine agreement again.
    [agreement_check] = determine_agreement_checked(myDataPath,cfg);

    fprintf('Overall agreement = %1.2f, positive agreement = %1.2f, negative agreement = %1.2f \n',...
    agreement_check.OA, agreement_check.PA, agreement_check.NA)
end

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
ccep2.stimchans_all = dataBase2stim.cc_stimchans_all;
ccep2.stimchans_avg = dataBase2stim.cc_stimchans_avg;
ccep2.stimpnames_all = dataBase2stim.stimpnames_all;
ccep2.stimpnames_avg = dataBase2stim.stimpnames_avg;
ccep2.stimsets_all = dataBase2stim.cc_stimsets_all;
ccep2.stimsets_avg = dataBase2stim.cc_stimsets_avg;
ccep2.dataName = dataBase2stim.dataName;
ccep2.ch = dataBase2stim.ch;
ccep2.tt = dataBase2stim.tt;
%ccep2.epoch_sorted_avg = dataBase2stim.cc_epoch_sorted_avg;
%ccep2.epoch_sorted_select_avg = dataBase2stim.cc_epoch_sorted_select_avg;

if strcmp(savefiles,'y')
    save([targetFolder,fileName], 'ccep2');
    save([myDataPath.CCEP_allpat,fileName], 'ccep2');
end

% save all stims
fileName5=[extractBefore(filename,'_ieeg'),'_CCEP_10stims.mat'];
ccep10 = dataBaseallstim.ccep;
ccep10.stimchans_all = dataBaseallstim.cc_stimchans_all;
ccep10.stimchans_avg = dataBaseallstim.cc_stimchans_avg;
ccep10.stimpnames_all = dataBaseallstim.stimpnames_all;
ccep10.stimpnames_avg = dataBaseallstim.stimpnames_avg;
ccep10.stimsets_all = dataBaseallstim.cc_stimsets_all;
ccep10.stimsets_avg = dataBaseallstim.cc_stimsets_avg;
ccep10.dataName = dataBaseallstim.dataName;
ccep10.ch = dataBaseallstim.ch;
ccep10.tt = dataBaseallstim.tt;
%ccep10.epoch_sorted_avg = dataBaseallstim.cc_epoch_sorted_avg;          % Deze epoch_sorted heb ik nodig voor plot_all_ccep_and_av maar hierdoor duurt het opslaan mega lang
%ccep10.epoch_sorted_select_avg = dataBaseallstim.cc_epoch_sorted_select_avg;

if strcmp(savefiles,'y')
    save([targetFolder,fileName5], 'ccep10');
    save([myDataPath.CCEP_allpat,fileName5], 'ccep10');

    fprintf('CCEPs 2stims and 10stims is saved in %s \n and %s',targetFolder , myDataPath.CCEP_allpat);
end


