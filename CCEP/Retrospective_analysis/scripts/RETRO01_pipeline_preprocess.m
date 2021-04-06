%% this script loads all runs of one subject, 
% divides the data into epochs around stimuli,
% averages all 10 stimuli or two stimuli (first positive (F1-F2) and first negative (F2-F1) stimulus)
% detects N1s after each averaged stimulus
% and saves it

clear; 

% set paths
cfg.mode = 'retro';

% patient characteristics
cfg.sub_labels = {['sub-' input('Patient number (RESPXXXX): ','s')]};
cfg.ses_label = input('Session number (ses-X): ','s');
cfg.task_label = 'task-SPESclin';

myDataPath = setLocalDataPath(cfg);         % When retrospective analysis is run, folders of prospective are removed from path.    

% Choose patient
config_CCEP

%% Load ECOG data
% Find runs
files = dir(fullfile(myDataPath.dataPath,cfg.sub_labels{1}, cfg.ses_label,'ieeg',...
    [cfg.sub_labels{1} '_' cfg.ses_label '_' cfg.task_label '_*'  '_events.tsv']));
names = {files.name};

% pre-allocation
strings = cell(size(names));
for n = 1:size(names,2)
    strings{n} = names{n}(strfind(names{n},'run-'):strfind(names{n},'run-')+9);
end

% Load data 
for R = 1:size(strings,2)
    cfg.run_label = strings(R);
    dataBase(R) = load_ECoGdata(cfg,myDataPath); 
end

fprintf('...Runs of Subject %s have run...\n',cfg.sub_labels{1})


%% Filter
% When this is used, dataBase.data will change into the fltered data
% DataBase.raw_data will not be changed and will keep the raw data

% 50 Hz filter and 120 low pass filter
dataBase = filter_bedArt(dataBase);
% 
fprintf('Subject %s is filtered. \n',cfg.sub_labels{1})
% 
%% CCEP for 2 and 10 stimulations
% avg_stim : write down the number of stimuli you want to average

% save only first stimulus in both directions
for i = 1:size(dataBase,2)
    avg_stim = 1;
    dataBase2stim(i,:) = preprocess_ECoG_spes(dataBase(i),cfg,avg_stim);
    tt = dataBase2stim.tt;
end

% save all stimuli (5) in each direction
for i = 1:size(dataBase,2)
    avg_stim = 5;
    dataBaseallstim(i,:) = preprocess_ECoG_spes(dataBase(i),cfg,avg_stim);
end
% When SPES was ran in multiple runs it has to be merge to combine all stimulations into one file.
% Since dataBase2stim and dataBaseallstim are based on the same stimulation, it does not matter whether you determine the size of
% dataBase2stims or dataBaseallstims
if size(dataBase2stim,1) >1         
    dataBase2stim = merge_runs(dataBase2stim);
    dataBaseallstim = merge_runs(dataBaseallstim);       
end

fprintf('...%s has been preprocessed... \n',dataBase(1).sub_label)

%% Do a quick check by visualizing the stimuli of 2 stims and 10 stims.
chan = 16; stim=51;
figure, 
sgtitle(sprintf('Patient %s, Stimulation pair %s, on %s',dataBase2stim.sub_label, dataBase2stim.stimpnames_avg{stim},dataBase2stim.ch{chan}))

subplot(2,1,1),

plot(tt,squeeze(dataBase2stim.cc_epoch_sorted_select_avg(chan,stim,:,:))','Color',[0.8 0.8 0.8],'LineWidth',1)
hold on
plot(tt,squeeze(dataBase2stim.cc_epoch_sorted_avg(chan,stim,:)),'k','LineWidth',2)
hold off
title('2 signals')
xlabel('time (s)')
xlim([-.2 1.0])
ylim([-300 300])
xlabel('Time (s)')
ylabel('Potential (\muV)')

           
subplot(2,1,2),
plot(tt,squeeze(dataBaseallstim.cc_epoch_sorted_select_avg(chan,stim,1:5,:))','Color',[0.8 0.8 0.8],'LineWidth',1)
hold on
plot(tt,squeeze(dataBaseallstim.cc_epoch_sorted_select_avg(chan,stim,6:10,:))','Color',[0.8 0.8 0.8],'LineWidth',1)
hold on
plot(tt,squeeze(dataBaseallstim.cc_epoch_sorted_avg(chan,stim,:)),'k','LineWidth',2)
hold off
title('All signals')
xlabel('time (s)')
xlim([-.2 1.0])
ylim([-300 300])

xlabel('Time (s)')
ylabel('Potential (\muV)')


%%
figure('Position',[256,618,677,357])
plot(tt,squeeze(dataBase2stim.cc_epoch_sorted_select_avg(chan,stim,1,:))','Color','b','LineWidth',1)
hold on
plot(tt,squeeze(dataBase2stim.cc_epoch_sorted_select_avg(chan,stim,2,:))','Color','r','LineWidth',1)
hold on
plot(tt,squeeze(dataBase2stim.cc_epoch_sorted_avg(chan,stim,:)),'k','LineWidth',3)
hold off
title(sprintf('First positive and negative response \n%s, evoked by %s',dataBase2stim.ch{chan},dataBase2stim.stimpnames_avg{stim}))
line([0 0],[-1200 350],'LineStyle','--','LineWidth',1.5,'color','k')

xlabel('time (s)')
xlim([-0.1 1])
ylim([-1200 350])
xlabel('Time (s)')
ylabel('Potential (\muV)')
lgd = legend('Positive direction','Negative direction');
lgd.FontWeight = 'bold';
lgd.Location = 'SouthEast';


%% Use the automatic N1 detector to detect ccep 
dataBase2stim = detect_n1peak_ECoG_ccep(dataBase2stim,cfg);
dataBaseallstim = detect_n1peak_ECoG_ccep(dataBaseallstim,cfg);
dataBase2stim.NmbrofStims = '2_stims';
dataBaseallstim.NmbrofStims = '10_stims';

disp('Detection of ERs is completed')                

%% visually check the automatically detected CCEPs
% 
% dataBase2stim = visualRating_ccep(dataBase2stim);
% dataBaseallstim = visualRating_ccep(dataBaseallstim);


%% save ccep
savefiles = input('Do you want to save the ccep-structures? [y/n] ','s');

targetFolder = [fullfile(myDataPath.CCEPpath, dataBase(1).sub_label,dataBase(1).ses_label,dataBase(1).run_label),'/'];

% Create the folder if it doesn't exist already.
if ~exist(targetFolder, 'dir')
    mkdir(targetFolder);
end

[~,filename,~] = fileparts(dataBase(1).dataName);

% save 2 stims
fileName=[extractBefore(filename,'_ieeg'),'_CCEP_2stims_filtered.mat'];
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
ccep2.SOZ = dataBase(1).tb_electrodes.soz;

if strcmp(savefiles,'y')
    save([targetFolder,fileName], 'ccep2');
    save([myDataPath.CCEP_allpat,fileName], 'ccep2');
end

% save all stims
fileName5=[extractBefore(filename,'_ieeg'),'_CCEP_10stims_filtered.mat'];
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
ccep10.SOZ = dataBase(1).tb_electrodes.soz;


if strcmp(savefiles,'y')
    save([targetFolder,fileName5], 'ccep10');
    save([myDataPath.CCEP_allpat,fileName5], 'ccep10');

    fprintf('CCEPs 2stims and 10stims is saved in %s \n and %s',targetFolder , myDataPath.CCEP_allpat);
end


 
 