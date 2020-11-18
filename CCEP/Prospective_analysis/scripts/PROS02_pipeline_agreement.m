% Script pipeline_preproces.m should be performed first to obtain the correct documents.

clear; 
% %% Choose patient
ccep_allPat.sub_labels = {'sub-PRIOS01','sub-PRIOS02','sub-PRIOS03','sub-PRIOS04','sub-PRIOS05','sub-PRIOS06'};
%ccep_allPat.name = {[input('Patient number type (PRIOSXX): ','s')]};

% set paths
ccep_allPat.mode = 'pros';
myDataPath = setLocalDataPath(ccep_allPat);

%% Load all ccep files in the folder CCEP_files_allPat
files = dir(fullfile(myDataPath.CCEP_allpat));

dataBase = struct;                      


for i = 1:size(ccep_allPat.sub_labels,2)
    dataBase(i).sub_label = ccep_allPat.sub_labels{i}; 
    
     % Find rows with the sub_label of interest 
    respLoc = find(contains({files(:).name},ccep_allPat.sub_labels{i}));        %find(contains({files(:).name},'merged'));
  
     % load all both the SPESclin and SPESprop of the patient
    for j=1:size(respLoc,2)                                                      % number of rows with the run_label of interest
       if contains(files(respLoc(j)).name,'clin_filt_check.') 
          load(fullfile(files(respLoc(j)).folder,files(respLoc(j)).name));
          dataBase(i).ccep_clin = ccep_clin;
          dataBase(i).filenameClin = files(respLoc(j)).name;
          
       elseif contains(files(respLoc(j)).name,'prop_filt_check.') 
          load(fullfile(files(respLoc(j)).folder,files(respLoc(j)).name));
          dataBase(i).ccep_prop = ccep_prop;   
          dataBase(i).filenameProp = files(respLoc(j)).name;
       end
    end
end
%file_clin = dir(fullfile(myDataPath.CCEP_allpat,['sub-',ccep_allPat.name{1} '*_CCEP_clin.mat']));
%file_prop = dir(fullfile(myDataPath.CCEP_allpat,['sub-',ccep_allPat.name{1} '*_CCEP_prop.mat']));

%load(fullfile(file_clin.folder, file_clin.name));
%load(fullfile(file_prop.folder, file_prop.name));

%% determine the agreement between 2 and 10 stims per run
% The determine_agreement function is not only determining the agreement
% when 2 sessions are compared. It could be possible to compare more, but
% then the values for W, Z and XandY should be changed. 

% load file of continue with variables if they are already in workspace
for subj = 1:size(dataBase,2)
    
    if exist('ccep_prop','var') && exist('ccep_clin','var') % if you want to continue with the previous part
        runs(1).name = dataBase(subj).filenameClin;
        runs(1).ccep = dataBase(subj).ccep_clin;
        runs(1).sub_label = dataBase(subj).sub_label;

        runs(2).name = dataBase(subj).filenameProp;
        runs(2).ccep = dataBase(subj).ccep_prop;
        runs(2).sub_label = dataBase(subj).sub_label;

    else
        fprintf('WARNING: No runs are found');
    end
    
    % Determine the agreement between the two matching runs
    agreement = determine_agreement(runs);          % Deze agreement nog toevoegen aan ccep! handig voor visualize Gridstructure
    
    dataBase(subj).agreement = agreement;
    
    fprintf('%s Overall agreement = %1.2f, positive agreement = %1.2f, negative agreement = %1.2f \n',...
        dataBase(subj).sub_label, agreement.agreement_stim.OA, agreement.agreement_stim.PA, agreement.agreement_stim.NA)
end


%% Rewrite the adjacency matrices to calculate agreement parameters
close all;

for subj = 1:size(dataBase,2)
    dataBase(subj).ccep_prop = rewrite_Amat(dataBase(subj).ccep_prop, dataBase(subj).agreement.AmatProp);
    dataBase(subj).ccep_clin = rewrite_Amat(dataBase(subj).ccep_clin, dataBase(subj).agreement.AmatClin);
end

%% Determine the indegree, outdegree, Betweenness centrality, the number of ERs per stimpair and the number of ERs per electrode
% The variables are saved in an excel in the run folder of the subject number
close all;

for subj = 1:size(dataBase,2)
    dataBase(subj).agreement_parameter = agreement_parameters(dataBase(subj).agreement, ...
        dataBase(subj).ccep_prop, dataBase(subj).ccep_clin, myDataPath);
    
    dataBase(subj).statistics = statistical_agreement(myDataPath, dataBase(subj).agreement_parameter, dataBase(subj).ccep_clin);
end

%% Visualise the agreement in a scatter plot
    scatter_networkPar(dataBase,myDataPath)

%% load electrodes positions (xlsx/electrodes.tsv)
% database (ccep10) is only used for channels and stimpairs and these are
% equal for 2 and 10, so does not matter which database is used.

for subj = 1:size(dataBase,2)
    visualise_gridstructure(myDataPath, dataBase(subj).ccep_clin, dataBase(subj).ccep_prop, dataBase(subj).agreement_parameter);
end


%% Make boxplots of the latency and amplitude of the N1 peaks. 
    boxplot_N1_peak(dataBase, myDataPath)


close all

