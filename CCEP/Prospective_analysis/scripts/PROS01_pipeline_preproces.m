clear; 
clc;

% set paths
cfg.mode = 'pros';
myDataPath = setLocalDataPath(cfg);

% Choose patient
config_CCEP

%% Load ECOG data
% Find SPESclin and SPESprop 
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
    dataBase(R) = load_ECoGdata(cfg,myDataPath); %#ok<SAGROW>
end

fprintf('Data of subject %s is loaded. \n',cfg.sub_labels{1})
    
    
%% Filter
% When this is used, dataBase.data will change into the filtered data
% DataBase.raw_data will not be changed and will keep the raw data
dataBase = filter_bedArt(dataBase);

fprintf('Both runs of subject %s are filtered. \n',cfg.sub_labels{1})

%% CCEP for SPES-clin stimulations
% save all stimuli of clinical SPES

for i = 1:size(dataBase,2)
    dataBase(i).task_name = dataBase(i).dataName(strfind(dataBase(i).dataName,'task-'):strfind(dataBase(i).dataName,'_run')-1); 

    if ismember(dataBase(i).task_name,'task-SPESclin')
       cfg.minstim = 5;
       dataBase_clin(i,:) = preprocess_ECoG_spes(dataBase(i),cfg);       %#ok<SAGROW>
      
    elseif ismember(dataBase(i).task_name,'task-SPESprop')
        cfg.minstim = 1;
        dataBase_prop(i,:) = preprocess_ECoG_spes(dataBase(i),cfg);      %#ok<SAGROW>

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


%% Check whether SPESclin and SPESprop contain the same stimulation pairs
% Stimpairs and electrodes which are different in the clinical and propofol
% SPES are removed.
[dataBase_clin, dataBase_prop] = similar_stimpairs(dataBase_clin, dataBase_prop);


% Visualise the signal of specified stimpair and channel
chan =1; stim=2;

tt = dataBase_clin.tt;
figure('Position',[515,333,1034,707]); 
subplot(2,1,1),
plot(tt,squeeze(dataBase_prop.cc_epoch_sorted_select_avg(chan,stim,:,:))','Color',[0.8 0.8 0.8],'LineWidth',1)
hold on
plot(tt,squeeze(dataBase_prop.cc_epoch_sorted_avg(chan,stim,:)),'k','LineWidth',2)
hold off
title(sprintf('SPES prop, %s, %s, %s',dataBase(1).sub_label, dataBase_prop.stimpnames_avg{stim},  dataBase_prop.ch{chan}))
xlabel('time (s)'); xlim([-.2 0.5]); ylabel('Potential \muV');            
% Create patch to indicate the 9 ms interval
ax = gca; ylimits = ax.YTick;
patch([0 0.009 0.009 0],[min(ylimits) min(ylimits) max(ylimits) max(ylimits)],[0.6,0.2,0.2], 'EdgeAlpha',0)
alpha(0.1)                % set patches transparency
ylim([-600 600])
xlim([-0.01 0.3])

subplot(2,1,2),
plot(tt,squeeze(dataBase_clin.cc_epoch_sorted_select_avg(chan,stim,:,:))','Color',[0.8 0.8 0.8],'LineWidth',1)
hold on
plot(tt,squeeze(dataBase_clin.cc_epoch_sorted_avg(chan,stim,:)),'k','LineWidth',2)
hold off
title(sprintf('SPES clin, %s, %s, %s',dataBase(1).sub_label, dataBase_clin.stimpnames_avg{stim},  dataBase_clin.ch{chan}))
xlabel('time (s)'); xlim([-.2 0.5]); ylabel('Potential \muV');
% Create patch to indicate the 9 ms interval
ax = gca; ylimits = ax.YTick;
patch([0 0.009 0.009 0],[min(ylimits) min(ylimits) max(ylimits) max(ylimits)],[0.6,0.2,0.2], 'EdgeAlpha',0)
alpha(0.1)                % set patches transparency
ylim([-600 600])
xlim([-0.01 0.3])




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


%% Visually detect N2's
% Check the signals with a checked N1 if they have an N2.
% Can only be performed when ERs are already visually checked and the
% n1_peak_amplitude_check and n1_peak_sample_check files are saved.

VisCheck_n2 = input('Do you want to visually detect N2s? [y/n] ','s');


if strcmp(VisCheck_n2,'y')
    dataBase_clin = visualRating_N2(dataBase_clin);
    dataBase_prop = visualRating_N2(dataBase_prop);
end

disp('N2 peaks are checked')      

%% Visually check all stimuli of a SPES
% The ERs detected with the detector are shown with a blue dot.
% This is very time consuming! only perform when the ER-detector results
% are not as expected.
VisCheck_all = input('Do you want to visually check ALL signals? [y/n] ','s');

if strcmp(VisCheck_all,'y')
    dataBase_clin = visualRating_all_ccep(dataBase_clin);
    dataBase_prop = visualRating_all_ccep(dataBase_prop);
end

disp('All responses are checked') 

%% save ccep
savefiles = input('Do you want to save the ccep-structures? [y/n] ','s');

for i = 1:size(dataBase,2)
    if isequal(dataBase(i).task_name, 'task-SPESclin')
        targetFolder_clin = [fullfile(myDataPath.CCEPpath, dataBase(i).sub_label,dataBase(i).ses_label,dataBase(i).task_name),'/'];
        [~,filename_clin,~] = fileparts(dataBase(i).dataName);
    elseif isequal(dataBase(i).task_name, 'task-SPESprop')
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
% When visual check is performed and data is checked, then save with
% correct name (also ensures that checked file is not overwritten when
% file is completely run)
if isfield(dataBase_prop.ccep, 'n1_peak_amplitude_check')
    fileName_prop=[extractBefore(filename_prop,'_ieeg'),'_CCEP_prop_filt_check19022021.mat'];      
else 
    fileName_prop=[extractBefore(filename_prop,'_ieeg'),'_CCEP_prop_Check_21012021.mat'];    
end

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

if strcmp(savefiles,'y')
    save([targetFolder_prop,fileName_prop], 'ccep_prop');
    save([myDataPath.CCEP_allpat,fileName_prop], 'ccep_prop');
end


% Create the folder if it doesn't exist already.
if ~exist(targetFolder_clin, 'dir')
    mkdir(targetFolder_clin);
end

% save Clinical-SPES
% When N1s are visually checked save with check in the name
if isfield(dataBase_clin.ccep, 'n1_peak_amplitude_check')
    fileName_clin=[extractBefore(filename_clin,'_ieeg'),'_CCEP_clin_filt_check.mat'];           
    
else % When N1s are not visually checked.
    fileName_clin=[extractBefore(filename_clin,'_ieeg'),'_CCEP_clin_Check_21012021.mat'];              
end

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

if strcmp(savefiles,'y')
    save([targetFolder_clin,fileName_clin], 'ccep_clin');
    save([myDataPath.CCEP_allpat,fileName_clin], 'ccep_clin');

    fprintf('CCEPs are saved for SPESprop en SPESclin for subject %s \n' , dataBase(1).sub_label);
end


%% Plot the average signal of all electrodes per stimulation pair
% Create figures with a plot per stimulation pair with the averaged response per electrode 
% Easy compare between clinical-SPES and propofol-SPES

dataBase_clin.save_fig = input('Do you want plot all average signals to the stimulus per stimpair? [y/n] ','s');
if strcmp(dataBase_clin.save_fig, 'y')
    plot_all_ccep(dataBase_clin, dataBase_prop, myDataPath)
end



%% Determine the amplitude and latency of the P1 and the highest point before N1
% Necessary to determine the rise and fall times of the N1.
% Amplitude of the P1 is not correct. Latency is.

P1_latency(dataBase_clin, dataBase_prop,cfg, myDataPath);

disp('P1_latency is saved to be later used in PROS02_pipeline_agreement.') 

%% Unique occurence 
% Determine how often each stimulation pair is stimulated during
% propofol-SPES and clinical-SPES
% null hypothesis that x is normally distributed, results in 1 when the null hypothesis is rejected 
NorDisClin = lillietest(occ)     ;            
NorDisProp = lillietest(occ_prop);

XX = dataBase(1).tb_events.electrical_stimulation_site;
[uniqueXX, ~, J]=unique(XX) ;
occ = histcounts(J, 1:numel(uniqueXX));


mean_occ = median(occ);
N1_rise_fall(dataBase_clin, dataBase_prop, cfg,myDataPath);

XX_prop = dataBase(2).tb_events.electrical_stimulation_site;
[uniqueXX_prop, ~, J_prop]=unique(XX_prop) ;
occ_prop = histcounts(J_prop, 1:numel(uniqueXX_prop));

mean_occ_prop = mean(occ_prop);

