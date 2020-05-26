clear; 

%% Choose patient
config_CCEP

%% set paths
myDataPath = setLocalDataPath(cfg);

%% select run
% choose between available runs
files = dir(fullfile(myDataPath.dataPath,cfg.sub_labels{1}, cfg.ses_label,'ieeg',...
    [cfg.sub_labels{1} '_' cfg.ses_label '_' cfg.task_label '_*'  '_events.tsv']));
names = {files.name};
strings = cellfun(@(x) x(strfind(names{1},'run-'):strfind(names{1},'run-')+9), names, 'UniformOutput', false);
stringsz = [repmat('%s, ',1,size(strings,2)-1),'%s'];

cfg.run_label = {input(sprintf(['Choose one of these runs: \n' stringsz '\n'],strings{:}),'s')}; % Chosen run is in cfg.run_label

clear files names strings stringsz

%% load data

dataBase = load_ECoGdata(cfg,myDataPath);

%% CCEP for 2 and 10 stimulations
stimulations = [1,5];
for K = 1:length(stimulations)
    [stim_dataBase(K)] = preprocess_ECoG_spes(dataBase,cfg,stimulations(K));

    % detect ccep
    [stim_database(K)] = detect_n1peak_ECoG_ccep(stim_dataBase(K), cfg);
    stim_database(K).ccep.amplitude_thresh = cfg.amplitude_thresh;
    stim_database(K).ccep.n1_peak_range = cfg.n1_peak_range;
    stim_database(K).ccep.cc_stimsets = stim_dataBase(K).cc_stimsets;
    stim_database(K).ccep.ch = stim_dataBase(K).ch;
    stim_database(K).ccep.stimpnames = stim_dataBase(K).stimpnames;
    stim_database(K).ccep.stimchans = stim_dataBase(K).cc_stimchans; 
    stim_database(K).ccep.stimnum = stim_dataBase(K).stimnum;
end
 disp('Detection of ERs is completed')

%% save ccep
targetFolder = [myDataPath.CCEPpath, stim_database(1).sub_label,'/',stim_database(1).ses_label,'/', stim_database(1).run_label,'/'];

% Create the folder if it doesn't exist already.
if ~exist(targetFolder, 'dir')
    mkdir(targetFolder);
end

start_filename = strfind(stim_database(1).dataName,'/');
stop_filename = strfind(stim_database(1).dataName,'_ieeg');

    
if stim_database(1).stimnum == 1
    fileName2=[stim_database(1).dataName(start_filename(end)+1:stop_filename-1),'_CCEP_2stims.mat'];

    ccep = stim_database(1).ccep;
    ccep.dataName = stim_database(1).dataName;

    save([targetFolder,fileName2], 'ccep');
end

if stim_database(2).stimnum == 5
    fileName5=[stim_database(2).dataName(start_filename(end)+1:stop_filename-1),'_CCEP_10stims.mat'];

    ccep = stim_database(2).ccep;
    ccep.dataName = stim_database(1).dataName;

    save([targetFolder,fileName5], 'ccep');
end

fprintf('CCEPs is saved in %s%s \n',targetFolder);

%% detmine the agreement between 2 and 10 stims per run
% The determine_agreement function is not only determining the agreement
% when 2 sessions are compared. It coule be possible to compare more, but
% then the values for W, Z and XandY should be changed. 

[overall_agr, positive_agr, negative_agr,compare_mat] = determine_agreement(myDataPath,cfg,stim_database);
agreement.OA = overall_agr;
agreement.PA = positive_agr;
agreement.NA = negative_agr; 

%% Save the values for the agreement per run (2 and 10 stims)
targetFolder = [myDataPath.CCEPpath, stim_database(1).sub_label,'/',stim_database(1).ses_label,'/', stim_database(1).run_label,'/'];

% Create the folder if it doesn't exist already.
if ~exist(targetFolder, 'dir')
    mkdir(targetFolder);
end

Agreements = [stim_database(2).sub_label, '_', stim_database(2).run_label,'_agreement2_versus10.mat'];

save([targetFolder,Agreements], 'agreement');

fprintf('Agreemtents are saved in %s%s \n',targetFolder);

