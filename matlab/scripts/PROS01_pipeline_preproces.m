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


%% Remove empty rows in the dataBase structs, sometimes these are formed
% when there are multiple runs.
dataBase_clin = dataBase_clin(all(~cellfun(@isempty,struct2cell(dataBase_clin))));
dataBase_prop = dataBase_prop(all(~cellfun(@isempty,struct2cell(dataBase_prop))));

% When a SPES is ran in multiple runs, merge them.
if size(dataBase_clin,1)>1                                  % When SPES was ran in multiple runs
    dataBase_clinOrig = dataBase_clin;                      % enables check of the original databases
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



%% Save CC epoch sorted to be later used
% Saving the CCEP file takes approx 5 minutes

filename_clin = [dataBase_clin.sub_label,'_',dataBase_clin.task_label,'_CCEP.mat'];
filename_prop = [dataBase_prop.sub_label,'_',dataBase_prop.task_label,'_CCEP.mat'];

filefolder = fullfile(myDataPath.dataPath, 'derivatives', '/', 'CCEPs' ,'/');
if ~exist(filefolder,'dir')
    mkdir(filefolder)
end

% Save CCEPs for later analysis
save(fullfile(filefolder,filename_clin),'dataBase_clin','-v7.3');
save(fullfile(filefolder,filename_prop),'dataBase_prop','-v7.3');




%% Visually check detected cceps

% Check the average signal in which an ER was detected
VisCheck = input('Do you want to visually check the detected clinical-SPES CCEPs? [y/n] ','s');

% load checked N1s if visual rating has started earlier
if exist(fullfile(myDataPath.CCEPpath, dataBase_clin.sub_label, ...
        dataBase_clin.ses_label, dataBase_clin.task_label,...
        [dataBase_clin.sub_label, '_', dataBase_clin.ses_label,'_',...
        dataBase_clin.task_label,'_N1sChecked.mat']),'file')

    load(fullfile(myDataPath.CCEPpath, dataBase_clin.sub_label, ...
        dataBase_clin.ses_label, dataBase_clin.task_label,...
        [dataBase_clin.sub_label, '_', dataBase_clin.ses_label,'_',...
        dataBase_clin.task_label,'_N1sChecked.mat']));
   dataBase_clin.ccep = ccep; 
end


if strcmpi(VisCheck,'y')
   
    % continue with the stimulation pair after the last saved stimulation pair
    if sum(strcmp(fieldnames(dataBase_clin.ccep), 'checkUntilStimp')) == 1
        endstimp = ccep.checkUntilStimp;
    else
        endstimp = 0;
    end
    
    visualRating_ccep(dataBase_clin, cfg, endstimp, myDataPath);
 
else %% already save the automatically detected CCEP's
    fileName=[extractBefore(filename,'_run'),'_N1s.mat'];
    save([targetFolder,fileName], 'ccep');
        
end

disp('CCEPs are checked')      


%% Visual check of propofol-SPES
VisCheck = input('Do you want to visually check the detected propofol-SPES CCEPs? [y/n] ','s');

if exist(fullfile(myDataPath.CCEPpath, dataBase_prop.sub_label, ...
        dataBase_prop.ses_label, dataBase_prop.task_name,...
        [dataBase_prop.sub_label, '_', dataBase_prop.ses_label,'_',...
        dataBase_prop.task_name,'_N1sChecked.mat']),'file')
    
   load(fullfile(myDataPath.CCEPpath, dataBase_prop.sub_label, ...
        dataBase_prop.ses_label, dataBase_prop.task_name,...
        [dataBase_prop.sub_label, '_', dataBase_prop.ses_label,'_',...
        dataBase_prop.task_name,'_N1sChecked.mat']));  
   dataBase_prop.ccep = ccep;
end

if strcmpi(VisCheck,'y')
    
     % continue with the stimulation pair after the last saved stimulation pair
    if sum(strcmp(fieldnames(dataBase_prop.ccep), 'checkUntilStimp')) == 1
        endstimp_prop = dataBase_prop.ccep.checkUntilStimp;
    else
        endstimp_prop = 0;
    end
    
    visualRating_ccep(dataBase_prop, cfg, endstimp_prop, myDataPath);
    
end
  
disp('CCEPs are checked')      


