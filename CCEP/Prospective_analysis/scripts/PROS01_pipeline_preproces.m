clear; 

% Choose patient
config_CCEP

% set paths
cfg.mode = 'pros';
myDataPath = setLocalDataPath(cfg);

%% Load ECOG data
% Find if there are multiple runs
files = dir(fullfile(myDataPath.dataPath,cfg.sub_labels{1}, cfg.ses_label,'ieeg',...
    [cfg.sub_labels{1} '_' cfg.ses_label '_' cfg.task_label '_*'  '_events.tsv']));
names = {files.name};
strings = cell(size(names));

% Find run_labels
for n = 1:size(names,2)
    strings{n} = names{n}(strfind(names{n},'run-'):strfind(names{n},'run-')+9);
end

% Load data (also possible for multiple runs)
for R = 1:size(strings,2)
    cfg.run_label = strings(R);
    dataBase(R) = load_ECoGdata(cfg,myDataPath);
end

    fprintf('Both runs of subject %s have run. \n',cfg.sub_labels{1})
    
    
%% Filter
% When this is used, dataBase.data will change into the fltered data
% DataBase.raw_data will not be changed and will keep the raw data
dataBase = filter_bedArt(dataBase);

fprintf('Both runs of subject %s are filtered. \n',cfg.sub_labels{1})
%% CCEP for SPES-clin stimulations
% save all stimuli of clinical SPES

for i = 1:size(dataBase,2)
    dataBase(i).task_name = dataBase(i).dataName(strfind(dataBase(i).dataName,'task-'):strfind(dataBase(i).dataName,'_run')-1); 

    if ismember(dataBase(i).task_name,'task-SPESclin')
%        avg_stim_clin = 5;
       cfg.minstim = 5;
%        cfg.max_stim = 5;
       dataBase_clin(i,:) = preprocess_ECoG_spes(dataBase(i),cfg);
      
    elseif ismember(dataBase(i).task_name,'task-SPESprop')
%         avg_stim = 1;                                                         %%% This was necessary for Retro, NOT for prospective
        cfg.minstim = 1;
%         cfg.max_stim = 1;
        dataBase_prop(i,:) = preprocess_ECoG_spes(dataBase(i),cfg);     

    end
end

fprintf('...%s has been epoched and averaged... \n',cfg.sub_labels{1})

% Remove empty rows in the dataBase structs, sometimes these are formed
% when there are multiple runs.
dataBase_clin = dataBase_clin(all(~cellfun(@isempty,struct2cell(dataBase_clin))));
dataBase_prop = dataBase_prop(all(~cellfun(@isempty,struct2cell(dataBase_prop))));


% When a SPES is ran in multiple runs, merge them.
if size(dataBase_clin,1)>1          % When SPES was ran in multiple runs
     dataBase_clin = merge_runs(dataBase_clin); 
elseif size(dataBase_prop,1) >1 
     dataBase_prop = merge_runs(dataBase_prop);
     
end

tt = dataBase_clin.tt;

% check whether similar stimuli are present in the same stimulus pair
chan = 6; stim=12;
figure, 
subplot(2,1,1),
plot(tt,squeeze(dataBase_prop.cc_epoch_sorted_select_avg(chan,stim,:,:))','Color',[0.8 0.8 0.8],'LineWidth',1)
hold on
plot(tt,squeeze(dataBase_prop.cc_epoch_sorted_avg(chan,stim,:)),'k','LineWidth',2)
hold off
title('SPES prop')
xlabel('time (s)')
xlim([-.2 1.0])
            
subplot(2,1,2),
plot(tt,squeeze(dataBase_clin.cc_epoch_sorted_select_avg(chan,stim,1:5,:))','Color','r','LineWidth',1)
hold on
plot(tt,squeeze(dataBase_clin.cc_epoch_sorted_select_avg(chan,stim,6:10,:))','Color','b','LineWidth',1)
hold on
plot(tt,squeeze(dataBase_clin.cc_epoch_sorted_avg(chan,stim,:)),'k','LineWidth',2)
hold off
title('SPES clin')
xlabel('time (s)')
xlim([-.2 1.0])


figure()
plot(tt,squeeze(dataBase_clin.cc_epoch_sorted_select_avg(chan,stim,1:5,:))','Color','r','LineWidth',1)
hold on
plot(tt,squeeze(dataBase_clin.cc_epoch_sorted_select_avg(chan,stim,6:10,:))','Color','b','LineWidth',1)
hold on
plot(tt,squeeze(dataBase_clin.cc_epoch_sorted_avg(chan,stim,:)),'k','LineWidth',2)
hold off
title('SPES clin')
xlabel('time (s)')
xlim([-.1 0.1])


%% Check whether SPESclin and SPESprop contain the same stimulation pairs
[dataBase_clin, dataBase_prop] = similar_stimpairs(dataBase_clin, dataBase_prop);

%% Use the automatic N1 detector to detect ccep 

dataBase_clin = detect_n1peak_ECoG_ccep(dataBase_clin,cfg);
dataBase_prop = detect_n1peak_ECoG_ccep(dataBase_prop,cfg);

disp('Detection of ERs is completed')

%% Visually check detected cceps
% Check the average signal in which an ER was detected
VisCheck = input('Do you want to visually check the detected CCEPs? [y/n] ','s');

if strcmp(VisCheck,'y')
    dataBase_clin = visualRating_ccep(dataBase_clin);
    dataBase_prop = visualRating_ccep(dataBase_prop);
end

disp('CCEPs are checked')      

%% save ccep
savefiles = input('Do you want to save the ccep-structures? [y/n] ','s');

for i = 1:size(dataBase,2)
    if dataBase(i).task_name == 'task-SPESclin'
        targetFolder_clin = [fullfile(myDataPath.CCEPpath, dataBase(i).sub_label,dataBase(i).ses_label,dataBase(i).task_name),'/'];
        [~,filename_clin,~] = fileparts(dataBase(i).dataName);
    elseif dataBase(i).task_name == 'task-SPESprop'
        targetFolder_prop = [fullfile(myDataPath.CCEPpath, dataBase(i).sub_label,dataBase(i).ses_label,dataBase(i).task_name),'/'];
        [~,filename_prop,~] = fileparts(dataBase(i).dataName);
    end
end

% Create the folder if it doesn't exist already.
if ~exist(targetFolder_prop, 'dir')
    mkdir(targetFolder_prop);
end

% [~,filename,~] = fileparts(dataBase(1).dataName);

% save propofol SPES
fileName_prop=[extractBefore(filename_prop,'_ieeg'),'_CCEP_prop.mat'];
ccep_prop = dataBase_prop.ccep;
ccep_prop.stimchans_all = dataBase_prop.cc_stimchans_all;
ccep_prop.stimchans_avg = dataBase_prop.cc_stimchans_avg;
ccep_prop.stimpnames_all = dataBase_prop.stimpnames_all;
ccep_prop.stimpnames_avg = dataBase_prop.stimpnames_avg;
ccep_prop.stimsets_all = dataBase_prop.cc_stimsets_all;
ccep_prop.stimsets_avg = dataBase_prop.cc_stimsets_avg;
ccep_prop.dataName = dataBase_prop.dataName;
ccep_prop.ch = dataBase_prop.ch;
ccep_prop.tt = dataBase_prop.tt;
%ccep2.epoch_sorted_avg = dataBase2stim.cc_epoch_sorted_avg;
%ccep2.epoch_sorted_select_avg = dataBase2stim.cc_epoch_sorted_select_avg;

if strcmp(savefiles,'y')
    save([targetFolder_prop,fileName_prop], 'ccep_prop');
    save([myDataPath.CCEP_allpat,fileName_prop], 'ccep_prop');
end


% Create the folder if it doesn't exist already.
if ~exist(targetFolder_clin, 'dir')
    mkdir(targetFolder_clin);
end

% save all stims
fileName_clin=[extractBefore(filename_clin,'_ieeg'),'_CCEP_clin.mat'];
ccep_clin = dataBase_clin.ccep;
ccep_clin.stimchans_all = dataBase_clin.cc_stimchans_all;
ccep_clin.stimchans_avg = dataBase_clin.cc_stimchans_avg;
ccep_clin.stimpnames_all = dataBase_clin.stimpnames_all;
ccep_clin.stimpnames_avg = dataBase_clin.stimpnames_avg;
ccep_clin.stimsets_all = dataBase_clin.cc_stimsets_all;
ccep_clin.stimsets_avg = dataBase_clin.cc_stimsets_avg;
ccep_clin.dataName = dataBase_clin.dataName;
ccep_clin.ch = dataBase_clin.ch;
ccep_clin.tt = dataBase_clin.tt;
%ccep10.epoch_sorted_avg = dataBaseallstim.cc_epoch_sorted_avg;          % Deze epoch_sorted heb ik nodig voor plot_all_ccep_and_av maar hierdoor duurt het opslaan mega lang
%ccep10.epoch_sorted_select_avg = dataBaseallstim.cc_epoch_sorted_select_avg;

if strcmp(savefiles,'y')
    save([targetFolder_clin,fileName_clin], 'ccep_clin');
    save([myDataPath.CCEP_allpat,fileName_clin], 'ccep_clin');

    fprintf('CCEPs are saved for SPESprop en SPESclin for subject %s \n' , dataBase(1).sub_label);
end



%% Plot the average signal of all electrodes per stimulation pair
dataBase_clin.save_fig = input('Do you want plot all average signals to the stimulus per stimpair? [y/n] ','s');
if strcmp(dataBase_clin.save_fig, 'y')
    plot_all_ccep(dataBase_clin, dataBase_prop, myDataPath)
end
% plot_all_ccep( dataBase_prop, myDataPath)

