%% detection and visual scoring of ERs
% author: Dorien van Blooijs, Sifra Blok
% date: September 2019, April 2020
% load ECoG data, split into stimulation trials, re-reference the data,
% detect N1-peaks after stimulation, visual check of detected N1s, save
% these for further analysis and display

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
    [cfg.sub_labels{1} '_' cfg.ses_label '_*'  '_events.tsv']));
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

% housekeeping
clear R n files strings
    
%% CCEP for SPES-clin stimulations
clear dataBase_clin dataBase_prop

% save all stimuli of clinical SPES or propofol SPES
for i = 1:size(dataBase,2)
    dataBase(i).task_name = dataBase(i).dataName(strfind(dataBase(i).dataName,'task-'):strfind(dataBase(i).dataName,'_run')-1); 

    if ismember(dataBase(i).task_name,'task-SPESclin')
       cfg.minstim = 5;
       dataBase_clin(i,:) = preprocess_ECoG_spes(dataBase(i),cfg);       %#ok<SAGROW>
        fprintf('...%s, %s, %s, has been epoched and averaged... \n',dataBase(i).sub_label, dataBase(i).run_label, dataBase(i).task_name)

    elseif ismember(dataBase(i).task_name,'task-SPESprop')
        cfg.minstim = 1;
        dataBase_prop(i,:) = preprocess_ECoG_spes(dataBase(i),cfg);      %#ok<SAGROW>
        fprintf('...%s, %s, %s, has been epoched and averaged... \n',dataBase(i).sub_label, dataBase(i).run_label, dataBase(i).task_name)

    end
end


% Remove empty rows in the dataBase structs, sometimes these are formed
% when there are multiple runs.
dataBase_clin = dataBase_clin(all(~cellfun(@isempty,struct2cell(dataBase_clin))));
dataBase_prop = dataBase_prop(all(~cellfun(@isempty,struct2cell(dataBase_prop))));

% When a SPES is ran in multiple runs, merge them.
if size(dataBase_clin,1)>1          % When SPES was ran in multiple runs
    dataBase_clinOrig = dataBase_clin; % enables check of the original databases
     dataBase_clin = merge_runs(dataBase_clin); 

elseif size(dataBase_prop,1) >1 
    dataBase_propOrig = dataBase_prop;
    dataBase_prop = merge_runs(dataBase_prop);
     
end

% housekeeping
clear i names

%% Re-reference data
% The reference calculatede below is the median of the signals that have a
% post- and pre-stim variance that is smaller than the pre-stim variance of
% the median signal of all signals (CAR). At least 5% of the signals have
% to be averaged in the reference signal. 

% Lowest variance is used to avoid introducing CCEP's in signals that
% originally do not have a CCEP.

cfg.reref = input('Do you want to rereference the data with an average of the 10 signals with the lowest variance? [y/n]: ','s'); % ('yes' = re-reference, 'no' = no re-reference)

if strcmp(cfg.reref,'y')
    tic
    dataBase_clin = rereference_with_lowest_var(dataBase_clin, cfg);
    
    dataBase_prop = rereference_with_lowest_var(dataBase_prop, cfg);
    
    toc

    disp('Data is rereferenced and all trials of each stimulus pair are averaged.')

else
   
    disp('Data is not re-referenced and all trials of each stimulus pair are averaged.')
end

%% Check whether SPESclin and SPESprop contain the same stimulation pairs
% Stimpairs and electrodes which are different in the clinical and propofol
% SPES are removed.
[dataBase_clin, dataBase_prop] = similar_stimpairs(dataBase_clin, dataBase_prop);

%% Use the automatic N1 detector to detect ccep 

dataBase_clin = detect_n1peak_ECoG_ccep(dataBase_clin,cfg);
dataBase_prop = detect_n1peak_ECoG_ccep(dataBase_prop,cfg);

disp('Detection of ERs is completed')

%% Visually check detected cceps

% Check the average signal in which an ER was detected
VisCheck = input('Do you want to visually check the detected clinical-SPES CCEPs? [y/n] ','s');

% load checked N1s if visual rating has started earlier
if exist(fullfile(myDataPath.CCEPpath, dataBase_clin.sub_label, ...
        dataBase_clin.ses_label, dataBase_clin.task_label,...
        [dataBase_clin.sub_label, '_', dataBase_clin.ses_label,'_',...
        dataBase_clin.task_label,'_N1sChecked.mat']),'file')
    
   dataBase_clin.ccep = load(fullfile(myDataPath.CCEPpath, dataBase_clin.sub_label, ...
        dataBase_clin.ses_label, dataBase_clin.task_label,...
        [dataBase_clin.sub_label, '_', dataBase_clin.ses_label,'_',...
        dataBase_clin.task_label,'_N1sChecked.mat']));   
end


if strcmpi(VisCheck,'y')
   
    % continue with the stimulation pair after the last saved stimulation pair
    if sum(strcmp(fieldnames(dataBase_clin.ccep), 'checkUntilStimp')) == 1
        endstimp_clin = dataBase_clin.ccep.checkUntilStimp;
    else
        endstimp_clin = 0;
    end
    
    dataBase_clin = visualRating_ccep(dataBase_clin, cfg, endstimp_clin, myDataPath);
    
end

disp('CCEPs are checked')      


%% Visual check of propofol-SPES
VisCheck = input('Do you want to visually check the detected propofol-SPES CCEPs? [y/n] ','s');

if exist(fullfile(myDataPath.CCEPpath, dataBase_prop.sub_label, ...
        dataBase_prop.ses_label, dataBase_prop.task_name,...
        [dataBase_prop.sub_label, '_', dataBase_prop.ses_label,'_',...
        dataBase_prop.task_name,'_N1sChecked.mat']),'file')
    
   dataBase_prop.ccep = load(fullfile(myDataPath.CCEPpath, dataBase_prop.sub_label, ...
        dataBase_prop.ses_label, dataBase_prop.task_name,...
        [dataBase_prop.sub_label, '_', dataBase_prop.ses_label,'_',...
        dataBase_prop.task_name,'_N1sChecked.mat']));   
end

if strcmpi(VisCheck,'y')
    
     % continue with the stimulation pair after the last saved stimulation pair
    if sum(strcmp(fieldnames(dataBase_prop.ccep), 'checkUntilStimp')) == 1
        endstimp_prop = dataBase_prop.ccep.checkUntilStimp;
    else
        endstimp_prop = 0;
    end
    
    dataBase_prop = visualRating_ccep(dataBase_prop, cfg, endstimp_prop, myDataPath);
    
end
  
disp('CCEPs are checked')      

%% Display doubtfull observations
% When still the visual check off the N1's was just performed, ccep_prop is still a variable. 
% When the doubtfull observations are looked at at a later time, then you
% can load them.
% if exist('ccep_clin','var')
%     doubt_clin = ccep_clin.obs_tab;
% else
%     saved_ccep_clin = load(fullfile(myDataPath.CCEPpath, dataBase_clin.sub_label, ...
%         dataBase_clin.ses_label, dataBase_clin.task_name,...
%         [dataBase_clin.sub_label, '_', dataBase_clin.ses_label,'_',...
%         dataBase_clin.task_name,'_',dataBase_clin.run_label,'_CCEP_clin_reref_check.mat']));  
%      
%     ccep_clin = saved_ccep_clin.ccep_clin;
%     doubt_clin = ccep_clin.obs_tab;
% 
% end
%    
% ccep_clin = recheck_doubtfull_cceps(doubt_clin, dataBase_clin, ccep_clin, cfg);
% 
% 
% if exist('ccep_prop','var')
%     doubt_prop = ccep_prop.obs_tab;
% else
%     saved_ccep_prop = load(fullfile(myDataPath.CCEPpath, dataBase_prop.sub_label, ...
%         dataBase_prop.ses_label, dataBase_prop.task_name,...
%         [dataBase_prop.sub_label, '_', dataBase_prop.ses_label,'_',...
%         dataBase_prop.task_name,'_',dataBase_prop.run_label,'_CCEP_prop_reref_check.mat']));   
%     
%     ccep_prop = saved_ccep_prop.ccep_prop;
%     doubt_prop = ccep_prop.obs_tab;
%  
% end
% 
% ccep_prop = recheck_doubtfull_cceps(doubt_prop, dataBase_prop, ccep_prop, cfg);
% 
% 


%% Visually detect N2's
% Check the signals with a checked N1 if they have an N2.
% Can only be performed when ERs are already visually checked and the
% n1_peak_amplitude_check and n1_peak_sample_check files are saved.

%%% IK HEB DIT NOG NIET GEUPDATE MET DE REREF DATA, AANGEZIEN DE N2 INFO MISSCHIEN NOG WAT MINDER RELEVANT IS NU VOOR HET CONGRES. 

% VisCheck_n2 = input('Do you want to visually detect N2s? [y/n] ','s');
% 
% 
% if strcmp(VisCheck_n2,'y')
%     dataBase_clin = visualRating_N2(dataBase_clin);
%     dataBase_prop = visualRating_N2(dataBase_prop);
% end
% 
% disp('N2 peaks are checked')      

%% save ccep
savefiles = input('Do you want to save the ccep-structures? [y/n] ','s');

for i = 1:size(dataBase,2)
    if isequal(dataBase(i).task_label, 'task-SPESclin')
        targetFolder = [fullfile(myDataPath.CCEPpath, dataBase(i).sub_label,dataBase(i).ses_label,dataBase(i).task_label),'/'];
        [~,filename,~] = fileparts(dataBase(i).dataName);

        dataBase_temp = dataBase_clin;
        modus = 'clin';

    elseif isequal(dataBase(i).task_label, 'task-SPESprop')
        targetFolder = [fullfile(myDataPath.CCEPpath, dataBase(i).sub_label,dataBase(i).ses_label,dataBase(i).task_label),'/'];
        [~,filename,~] = fileparts(dataBase(i).dataName);

        dataBase_temp = dataBase_prop;
        modus = 'prop';
    end

    % Create the folder if it doesn't exist already.
    if ~exist(targetFolder, 'dir')
        mkdir(targetFolder);
    end

    % save SPES
    % When N1s are visually checked save with check in the name
    if isfield(dataBase_temp.ccep, 'n1_peak_amplitude_check')
        fileName=[extractBefore(filename,'_run'),'_N1sChecked.mat'];

    else % When N1s are not visually checked.
        fileName=[extractBefore(filename,'_run'),'_N1s.mat'];
    end

    ccep = dataBase_temp.ccep;
    ccep.stimchans_all = dataBase_temp.cc_stimchans_all;
    ccep.stimchans_avg = dataBase_temp.cc_stimchans_avg;
    ccep.stimpnames_all = dataBase_temp.stimpnames_all;
    ccep.stimpnames_avg = dataBase_temp.stimpnames_avg;
    ccep.stimsets_all = dataBase_temp.cc_stimsets_all;
    ccep.stimsets_avg = dataBase_temp.cc_stimsets_avg;
    ccep.dataName = dataBase_temp.dataName;
    ccep.ch = dataBase_temp.ch;
    ccep.tt = dataBase_temp.tt;
    ccep.dir = cfg.dir;
    ccep.amp = cfg.amp;
    ccep.epoch_length = cfg.epoch_length;
    ccep.epoch_prestim = cfg.epoch_prestim;
    ccep.reref = cfg.reref;

    if strcmp(savefiles,'y')
        save([targetFolder,fileName], 'ccep');
%         save([myDataPath.CCEP_allpat,fileName], 'ccep');

    end

end

fprintf('CCEPs are saved for SPESprop en SPESclin for subject %s \n' , dataBase(1).sub_label);


%% Plot the average signal of all electrodes per stimulation pair
% Create figures with a plot per stimulation pair with the averaged response per electrode 
% Easy compare between clinical-SPES and propofol-SPES
% 
% dataBase_clin.save_fig = input('Do you want plot all average signals to the stimulus per stimpair? [y/n] ','s');
% if strcmp(dataBase_clin.save_fig, 'y')
%     plot_all_ccep(dataBase_clin, dataBase_prop, myDataPath)
% end



%% Determine the amplitude and latency of the P1 and the highest point before N1
% Necessary to determine the rise and fall times of the N1.
% Amplitude of the P1 is not correct. Latency is.

% P1_latency(dataBase_clin, dataBase_prop,cfg, myDataPath);
% 
% disp('P1_latency is saved to be later used in PROS02_pipeline_agreement.') 

%% Unique occurence 
% Determine how often each stimulation pair is stimulated during
% propofol-SPES and clinical-SPES
% % null hypothesis that x is normally distributed, results in 1 when the null hypothesis is rejected 
% NorDisClin = lillietest(occ)     ;            
% NorDisProp = lillietest(occ_prop);
% 
% XX = dataBase(1).tb_events.electrical_stimulation_site;
% [uniqueXX, ~, J]=unique(XX) ;
% occ = histcounts(J, 1:numel(uniqueXX));
% 
% 
% mean_occ = median(occ);
% N1_rise_fall(dataBase_clin, dataBase_prop, cfg,myDataPath);
% 
% XX_prop = dataBase(2).tb_events.electrical_stimulation_site;
% [uniqueXX_prop, ~, J_prop]=unique(XX_prop) ;
% occ_prop = histcounts(J_prop, 1:numel(uniqueXX_prop));
% 
% mean_occ_prop = mean(occ_prop);
% 
