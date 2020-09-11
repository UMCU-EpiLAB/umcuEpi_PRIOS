% Script pipeline_preproces.m should be performed first to obtain the correct documents.

clear; 
% %% Choose patient
cfg.sub_labels = {'sub-RESP0701','sub-RESP0702','sub-RESP0703','sub-RESP0706','sub-RESP0724','sub-RESP0728'}; %{[input('Patient number type (RESPXXXX or PRIOSXX): ','s')]};

% set paths
cfg.mode = 'retro';
myDataPath = setLocalDataPath(cfg);

%% Load all ccep files in the folder CCEP_files_allPat

files = dir(fullfile(myDataPath.CCEP_allpat));

dataBase = struct;
for i=1:size(cfg.sub_labels,2)
    dataBase(i).sub_label = cfg.sub_labels{i};
    
    respLoc = find(contains({files(:).name},cfg.sub_labels{i}));
    
    for j=1:size(respLoc,2)
       if contains(files(respLoc(j)).name,'10stims') 
          load(fullfile(files(respLoc(j)).folder,files(respLoc(j)).name));
          dataBase(i).ccep10 = ccep10;
          dataBase(i).filename10 = files(respLoc(j)).name;
          
       elseif contains(files(respLoc(j)).name,'2stims') 
          load(fullfile(files(respLoc(j)).folder,files(respLoc(j)).name));
          dataBase(i).ccep2 = ccep2;   
          dataBase(i).filename2 = files(respLoc(j)).name;
       end
    end
end

% small cleanup
clear respLoc i j files ccep10 ccep2

%% determine the agreement between 2 and 10 stims per run
% The determine_agreement function is not only determining the agreement
% when 2 sessions are compared. It could be possible to compare more, but
% then the values for W, Z and XandY should be changed. 
close 
clc

for subj = 1:size(dataBase,2)
    
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
    
    statistics = statistical_agreement(myDataPath, dataBase(subj).agreement_parameter, dataBase(subj).ccep10);
end

%% load electrodes positions (xlsx/electrodes.tsv)
% database (ccep10) is only used for channels and stimpairs and these are
% equal for 2 and 10, so does not matter which database is used.
visualise_gridstructure(myDataPath, ccep10, ccep2, agreement_parameter);

