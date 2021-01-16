clear; 
clc;

% Choose patient
config_CCEP

% set paths
cfg.mode = 'pros';
myDataPath = setLocalDataPath(cfg);

%% Load ECOG data
% Only load the clinical SPES

files = dir(fullfile(myDataPath.dataPath,cfg.sub_labels{1}, cfg.ses_label,'ieeg',...
    [cfg.sub_labels{1} '_' cfg.ses_label '_task-SPESclin_*'  '_events.tsv']));
names = {files.name};
strings = cell(size(names));

% Find run_labels
for n = 1:size(names,2)
    strings{n} = names{n}(strfind(names{n},'run-'):strfind(names{n},'run-')+9);
end

% Load data 
for R = 1:size(strings,2)
    cfg.run_label = strings(R);
    dataBase(R) = load_ECoGdata(cfg,myDataPath);
end

fprintf('Data of subject %s is loaded. \n',cfg.sub_labels{1})
    
    %% CCEP for SPES-clin stimulations
% save all stimuli of clinical SPES
dataBase.task_name = dataBase.dataName(strfind(dataBase.dataName,'task-'):strfind(dataBase.dataName,'_run')-1); 
cfg.minstim = 5;
dataBase_clin = preprocess_ECoG_spes(dataBase,cfg);
   
fprintf('...%s has been epoched and averaged... \n',cfg.sub_labels{1})

%% Use the automatic N1 detector to detect ccep 

dataBase_clin = detect_n1peak_ECoG_ccep(dataBase_clin,cfg);

disp('Detection of ERs is completed')

%% Interobserver difference check
% script to determine the latency of the N1's which are incorrectly
% rejected by rater-2

% Only run this script for PRIOS04 or PRIOS05
dataBase = check_interobs(dataBase_clin, myDataPath);


%% Determine the Cohen's Kappa interobserver variability
% Determine this with checked files of two raters/observers
interobserverKappa(myDataPath);
