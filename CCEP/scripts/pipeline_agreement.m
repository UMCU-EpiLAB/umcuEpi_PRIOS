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
strings = cellfun(@(x) x(strfind(names{1},'run-'):strfind(names{1},'run-')+9), names, 'UniformOutput', false);
stringsz = [repmat('%s, ',1,size(strings,2)-1),'%s'];

cfg.run_label = {input(sprintf(['Choose one of these runs: \n' stringsz '\n'],strings{:}),'s')}; % Chosen run is in cfg.run_label

if ~contains(cfg.run_label,'run-')
   error('"run-" is missing in run_label') 
end

clear files names strings stringsz

%% load data

dataBase = load_ECoGdata(cfg,myDataPath);

%% CCEP for 2 and 10 stimulations

% avg_stim = 1: we only want to use the first stimulation of all stimuli 
% in one stimulation pair. So instead of averaging all stimuli in one 
% stimulation pair, we only 'average' the first stimulus. If you want to 
% use all stimuli, use avg_stim = [];

% save only first stimulus in both directions
avg_stim = 1;
dataBase2stim = preprocess_ECoG_spes(dataBase,cfg,avg_stim);

% save all stimuli
avg_stim = [];
dataBaseallstim = preprocess_ECoG_spes(dataBase,cfg,avg_stim);

%% Visually check all averaged signals for ERs

dataBaseallstim = plot_all_signals(dataBaseallstim);
dataBase2stim = plot_all_signals(dataBase2stim);

% Save 
targetFolder = [myDataPath.CCEPpath, dataBase(1).sub_label,'/',dataBase(1).ses_label,'/', dataBase(1).run_label,'/'];

vis_checked_all = [dataBaseallstim.sub_label, '_',dataBase(1).run_label, '_avg_10stims_vis.mat'];
ccep = dataBaseallstim.ccep;
save([targetFolder,vis_checked_all], 'ccep')

vis_checked_2 = [dataBase2stim.sub_label, '_', dataBase(1).run_label, '_avg_2stims_vis.mat'];
ccep = dataBase2stim.ccep;
save([targetFolder,vis_checked_2], 'ccep')
fprintf('Checked files are saved in %s \n',targetFolder);

%% detect ccep in only first stimulus in both directions and all stimuli
dataBase2stim = detect_n1peak_ECoG_ccep(dataBase2stim,cfg);
dataBaseallstim = detect_n1peak_ECoG_ccep(dataBaseallstim,cfg);
dataBase2stim.NmbrofStims = '2_stims';
dataBaseallstim.NmbrofStims = '10_stims';

disp('Detection of ERs is completed')

%% save ccep
targetFolder = [fullfile(myDataPath.CCEPpath, dataBase(1).sub_label,dataBase(1).ses_label,dataBase(1).run_label),'/'];

% Create the folder if it doesn't exist already.
if ~exist(targetFolder, 'dir')
    mkdir(targetFolder);
end

start_filename = strfind(dataBase(1).dataName,'/');
stop_filename = strfind(dataBase(1).dataName,'_ieeg');

% save 2 stims
fileName=[dataBase2stim.dataName(start_filename(end)+1:stop_filename-1),'_CCEP_2stims.mat'];
ccep = dataBase2stim.ccep;
ccep.stimchans = dataBase2stim.cc_stimchans;
ccep.stimpnames = dataBase2stim.stimpnames;
ccep.stimsets = dataBase2stim.cc_stimsets;
ccep.dataName = dataBase2stim.dataName;
ccep.ch = dataBase2stim.ch;
save([targetFolder,fileName], 'ccep');

% save all stims
fileName5=[dataBaseallstim.dataName(start_filename(end)+1:stop_filename-1),'_CCEP_10stims.mat'];
ccep = dataBaseallstim.ccep;
ccep.stimchans = dataBaseallstim.cc_stimchans;
ccep.stimpnames = dataBaseallstim.stimpnames;
ccep.stimsets = dataBaseallstim.cc_stimsets;
ccep.dataName = dataBaseallstim.dataName;
ccep.ch = dataBaseallstim.ch;
save([targetFolder,fileName5], 'ccep');

fprintf('CCEPs 2stims and 10stims is saved in %s \n',targetFolder);


%% determine the agreement between 2 and 10 stims per run
% The determine_agreement function is not only determining the agreement
% when 2 sessions are compared. It could be possible to compare more, but
% then the values for W, Z and XandY should be changed. 
close all

[agreement_run, agreement_stim, compare_mat, dif_mat, TotERs10, TotERs2, TotOnesStim, Amat10, Amat2] = determine_agreement(myDataPath,cfg);

fprintf('Overall agreement = %1.2f, positive agreement = %1.2f, negative agreement = %1.2f \n',...
    agreement_run.OA, agreement_run.PA, agreement_run.NA)

% fprintf('Overall agreement = %1.2f, positive agreement = %1.2f, negative agreement = %1.2f \n',...
%     agreement_stim.OA, agreement_stim.PA, agreement_stim.NA)

%% Determine the location of the ones (ER vs. No-ER)

[FindOnes, LocOnes, stimchans] = find_ones(dataBaseallstim,agreement_run);
 
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

%% Plot all 10 stimuli and the average for the 10 stims and the 2 stims
% This combines the two scripts above
tic;
dataBaseallstim.save_fig = str2double(input('Do you want to save the figures? [yes = 1, no = 0]: ','s'));
plot_all_ccep_and_av(dataBaseallstim,dataBase2stim, myDataPath, LocOnes, stimchans, dif_mat, TotOnesStim);
toc;

%% Calculate agreement parameters
close all;
dataBase2stim = rewrite_Amat(dataBase2stim,Amat2);
dataBaseallstim = rewrite_Amat(dataBaseallstim,Amat2);

%%
[indegree, outdegree, BC, rank_stimp, rank_elec]  = agreement_parameters(Amat10, dataBaseallstim, stimchans);

[indegree_2, outdegree_2, BC_2, rank_stimp_2, rank_elec_2] = agreement_parameters(Amat2, dataBase2stim, stimchans);

targetFolder = [myDataPath.CCEPpath, dataBase(1).sub_label,'/',dataBase(1).ses_label,'/', dataBase(1).run_label,'/'];
fileName = ['degrees_',dataBase(1).sub_label,'.xlsx'];
sheet = 'Indegree all';
writetable(rank_elec,[targetFolder, fileName],'sheet',sheet,'WriteRowNames',true)

sheet =  'Indegree 2';
writetable(rank_elec_2,[targetFolder, fileName],'sheet',sheet,'WriteRowNames',true)
 
sheet = 'Outdegree all';
writetable(rank_stimp,[targetFolder, fileName],'sheet',sheet,'WriteRowNames',true)

sheet =  'Outdegree 2';
writetable(rank_stimp_2,[targetFolder, fileName],'sheet',sheet,'WriteRowNames',true)



