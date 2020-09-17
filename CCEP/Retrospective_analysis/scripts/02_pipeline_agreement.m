% Script pipeline_preproces.m should be performed first to obtain the correct documents.

clear; 
% %% Select all patients
cfg.sub_labels = {'sub-RESP0701','sub-RESP0702','sub-RESP0703','sub-RESP0706','sub-RESP0724','sub-RESP0728'}; %{[input('Patient number type (RESPXXXX or PRIOSXX): ','s')]};

% set paths
cfg.mode = 'retro';
myDataPath = setLocalDataPath(cfg);

%% Load all ccep files in the folder CCEP_files_allPat
% Be aware that for some patients, SPES is saved in multiple runs

files = dir(fullfile(myDataPath.CCEP_allpat));

% Find all run-labels (sometimes SPES is ran in multiple runs)
for k = 1:size(files,1)
    if contains(files(k).name, 'sub')  ;                                        % If row is not empty
        run_label{k,:} = extractBetween(files(k).name, 'clin_', '_CCEP');
    end
end

unique_runlabels = unique([run_label{:,:}])';                                   % Be aware that this is also ordered

% Create database with the CCEP information of all patients of all runs and
% both protocols (2 and 10 stims)
dataBase = struct;                      
for i=1:size(unique_runlabels,1)                                                % DataBase must have the size of the number of runs 
     dataBase(i).run_label = unique_runlabels{i};
        
    % Find rows with the run_label of interest to determine the sub_label
    runLoc = find(contains({files(:).name},unique_runlabels{i}));              
    dataBase(i).sub_label = extractBefore(files(runLoc(1)).name, '_ses')  ;      % sub_label is the same for every row the run_label is the same
   
    
    % load all both the 10 stimuli and 2 stimuli of the patient
    for j=1:size(runLoc,2)                                                      % number of rows with the run_label of interest
       if contains(files(runLoc(j)).name,'10stims') 
          load(fullfile(files(runLoc(j)).folder,files(runLoc(j)).name));
          dataBase(i).ccep10 = ccep10;
          dataBase(i).filename10 = files(runLoc(j)).name;
          
       elseif contains(files(runLoc(j)).name,'2stims') 
          load(fullfile(files(runLoc(j)).folder,files(runLoc(j)).name));
          dataBase(i).ccep2 = ccep2;   
          dataBase(i).filename2 = files(runLoc(j)).name;
       end
    end
end

% Sort the rows based on the subjects instead of the run_label order
[~,index] = sortrows({dataBase.sub_label}.'); 
dataBase = dataBase(index);

% small cleanup
clear runLoc k j files ccep10 ccep2 run_label

%% determine the agreement between 2 and 10 stims per run
% The determine_agreement function is not only determining the agreement
% when 2 sessions are compared. It could be possible to compare more, but
% then the values for W, Z and XandY should be changed. 
close 
clc

% CHECKEN OF DIT NOG WERKT NU ER PER PATIENT MEERDERE RUNS ZIJN!!!!
for subj = 1:size(dataBase,2)
    
    % Find the 2 runs matching.
    % DIT MOET DUS EIGENLIJK, VIND ALLE RUNS PER PATIENT, VOEG DIE SAMEN EN
    % BEPAAK DAN DE AGREEMENT
    runs(1).ccep = dataBase(subj).ccep10;
    runs(1).name = dataBase(subj).filename10;
    runs(1).sub_label = dataBase(subj).sub_label;
    runs(2).ccep = dataBase(subj).ccep2;
    runs(2).name = dataBase(subj).filename2;
    runs(1).sub_label = dataBase(subj).sub_label;
    
    agreement = determine_agreement(runs);          % Deze agreement nog toevoegen aan ccep! handig voor visualize Gridstructure
    
    dataBase(subj).agreement = agreement;
    
    fprintf('%s Overall agreement = %1.2f, positive agreement = %1.2f, negative agreement = %1.2f \n',...
        dataBase(subj).sub_label, agreement.agreement_run.OA, agreement.agreement_run.PA, agreement.agreement_run.NA)
end

% short clean up
clear subj agreement runs

%% Determine the location of the ones (ER vs. No-ER)
% ccep10 is only necessary for the channels and stimpairs and those are
% equal for 2 and 10 stimuli so does not matter which database is used.

for subj = 1:size(dataBase,2)
    ccep10 = dataBase(subj).ccep10;
    agreement = dataBase(subj).agreement;
    
    LocOnes = find_ones(ccep10, agreement.agreement_run);
    dataBase(subj).LocOnes = LocOnes;
end

% clean up
clear ccep10 agreement LocOnes subj

%% Plot all 10 stimuli and the average for the 10 stims and the 2 stims
% This does not work without epoch_sorted information
% ZIE NOTITIE IN PIPELINE_PREPROCES (SAVE CCEPS)
plot_fig = input('Do you want plot the 10 stimuli and the average signals? [y/n] ','s');

if strcmp(plot_fig,'y')
    ccep10.save_fig = str2double(input('Do you want to save the figures? [yes = 1, no = 0]: ','s'));
    plot_all_ccep_and_av(ccep10, ccep2, myDataPath, LocOnes, agreement);
end

%% Calculate agreement parameters

close all;

for subj = 1:size(dataBase,2)
    
    dataBase(subj).ccep2 = rewrite_Amat(dataBase(subj).ccep2, dataBase(subj).agreement.Amat2);
    dataBase(subj).ccep10 = rewrite_Amat(dataBase(subj).ccep10, dataBase(subj).agreement.Amat10);
end

%% Determine the indegree, outdegree, Betweenness centrality, the number of ERs per stimpair and the number of ERs per electrode
% The variables are saved in an excel in the run folder of the subject number
close all;
clc

for subj = 1:size(dataBase,2)
    dataBase(subj).agreement_parameter = agreement_parameters(dataBase(subj).agreement, ...
        dataBase(subj).ccep2, dataBase(subj).ccep10, myDataPath);
    
    dataBase(subj).statistics = statistical_agreement(myDataPath, dataBase(subj).agreement_parameter, dataBase(subj).ccep10);
end

%% Load electrodes positions (xlsx/electrodes.tsv)
plot_fig = 'n';                 % 'n' when not all ER responses per stim have to be plot, 'y' when you do want to plot all
close all;

for subj = 1:size(dataBase,2)
    visualise_gridstructure(myDataPath, dataBase(subj).ccep10, dataBase(subj).ccep2, dataBase(subj).agreement_parameter,plot_fig);
end
