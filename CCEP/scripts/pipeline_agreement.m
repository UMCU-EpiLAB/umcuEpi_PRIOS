clear; 
% test of het pushen naar github werkt.
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

TotalPosN1 = (size(dataBaseallstim.cc_stimsets_avg,1)*length(dataBaseallstim.ch)) - (size(dataBaseallstim.cc_stimsets_avg,1)*2);

%% Visually check all averaged signals for ERs

dataBaseallstim = plot_all_signals(dataBaseallstim);
dataBase2stim = plot_all_signals(dataBase2stim);

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
    
    fprintf('CCEPs 2stims and 10stims is saved in %s \n',targetFolder);
end

%% determine the agreement between 2 and 10 stims per run
% The determine_agreement function is not only determining the agreement
% when 2 sessions are compared. It could be possible to compare more, but
% then the values for W, Z and XandY should be changed. 
close all

% load file of continue with variables if they are already in workspace
if exist('ccep2','var') && exist('ccep10','var') % if you want to continue with the previous part
    runs(1).name = fileName5;
    runs(1).ccep = ccep10;
    
    runs(2).name = fileName;
    runs(2).ccep = ccep2;
    
else
    files = dir(fullfile(myDataPath.CCEPpath, cfg.sub_labels{1}, 'ses-*' ,cfg.run_label{1},...
        [cfg.sub_labels{1} '_ses-*_' cfg.task_label '_*'  '_CCEP_*.mat']));
    
    runs = struct;
    if isempty(files)
        fprintf('WARNING: No runs are found');
    else
        for i = 1:size(files,1)
            names = fullfile(files(i).folder, files(i).name);
            ccep = load(names);
            runs(i).name = names;
            runs(i).ccep = ccep.ccep;
        end
    end
end

agreement = determine_agreement(runs);

fprintf('Overall agreement = %1.2f, positive agreement = %1.2f, negative agreement = %1.2f \n',...
    agreement.agreement_run.OA, agreement.agreement_run.PA, agreement.agreement_run.NA)

% fprintf('Overall agreement = %1.2f, positive agreement = %1.2f, negative agreement = %1.2f \n',...
%     agreement_stim.OA, agreement_stim.PA, agreement_stim.NA)

%% Determine the location of the ones (ER vs. No-ER)

[FindOnes, LocOnes, stimchans] = find_ones(dataBaseallstim,agreement.agreement_run);
 
%% visually check detected cceps

dataBaseallstim = visualRating_ccep(dataBaseallstim, stimchans);
dataBase2stim = visualRating_ccep(dataBase2stim, stimchans);
% correct: y
% incorrect: n or enter

% Save the values for the agreement per run (2 and 10 stims)
targetFolder = [myDataPath.CCEPpath, dataBase(1).sub_label,'/',dataBase(1).ses_label,'/', dataBase(1).run_label,'/'];

checked_all = [dataBaseallstim.sub_label, '_',dataBase(1).run_label, '_CCEP_', dataBaseallstim.NmbrofStims ,'_checked.mat'];
ccep = dataBaseallstim.ccep;
save([targetFolder,checked_all], 'ccep')

checked_2 = [dataBase2stim.sub_label, '_', dataBase(1).run_label, '_CCEP_', dataBase2stim.NmbrofStims ,'_checked.mat'];
ccep = dataBase2stim.ccep;
save([targetFolder,checked_2], 'ccep')
fprintf('Checked files are saved in %s \n',targetFolder);

% Perform the determine agreement again.
[agreement_check] = determine_agreement_checked(myDataPath,cfg);

fprintf('Overall agreement = %1.2f, positive agreement = %1.2f, negative agreement = %1.2f \n',...
    agreement_check.OA, agreement_check.PA, agreement_check.NA)


%% Plot all 10 stimuli and the average for the 10 stims and the 2 stims
% This combines the two scripts above
tic;
dataBaseallstim.save_fig = str2double(input('Do you want to save the figures? [yes = 1, no = 0]: ','s'));
plot_all_ccep_and_av(dataBaseallstim,dataBase2stim, myDataPath, LocOnes, stimchans, agreement);
toc;

%% Calculate agreement parameters
close all;
dataBase2stim = rewrite_Amat(dataBase2stim,agreement.Amat2);
dataBaseallstim = rewrite_Amat(dataBaseallstim,agreement.Amat10);

%% Determine the indegree, outdegree, Betweenness centrality, the number of ERs per stimpair and the number of ERs per electrode
% The variables are saved in an excel in the run folder of the subject number
close all;
<<<<<<< HEAD
agreement_parameter = agreement_parameters(Amat10,Amat2, dataBaseallstim, dataBase2stim,myDataPath);
=======
agreement_parameter = agreement_parameters(agreement, dataBaseallstim, dataBase2stim,stimchans);
>>>>>>> upstream/master





%% Plot the average signal of the 2 stims or 10 stims
% 
% dataBaseallstim.save_fig = str2double(input('Do you want to save the figures? [yes = 1, no = 0]: ','s'));
% plot_ccep_av_stimp(dataBaseallstim,dataBase2stim, myDataPath, stimchans, LocOnes, TotOnesStim, dif_mat)
% 
% fprintf('All CCEPS average are saved');
% 
% %% Plot the all 10 stimuli per stimulation pair
% dataBaseallstim.save_fig = str2double(input('Do you want to save the figures? [yes = 1, no = 0]: ','s'));
% plot_all_ccep(dataBaseallstim, myDataPath, LocOnes, stimchans, dif_mat);
