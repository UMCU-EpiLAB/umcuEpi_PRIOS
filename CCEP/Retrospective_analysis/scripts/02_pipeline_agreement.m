% Script pipeline_preproces.m should be performed first to obtain the correct documents.

clear; 
% %% Choose patient
cfg.sub_labels = {'sub-RESP0701','sub-RESP0702','sub-RESP0703','sub-RESP0706','sub-RESP0724','sub-RESP0728'}; %{[input('Patient number type (RESPXXXX or PRIOSXX): ','s')]};

% set paths
cfg.mode = 'retro';
myDataPath = setLocalDataPath(cfg);

%% Load all ccep files in the folder CCEP_files_allPat

files = dir(fullfile(myDataPath.CCEP_allpat));


for k = 1:size(files,1)
    if contains(files(k).name, 'sub')                   % If row is not empty
        run_label{k,:} = extractBetween(files(k).name, 'clin_', '_CCEP')
        sub_label(k,:) = extractBefore(files(k).name, '_ses')
    end
end

unique_runlabels = unique([run_label{:,:}])';
[~,loc_sub_labels] = unique([sub_label(:,:)],'rows');

dataBase = struct;
for i=1:size(unique_runlabels,1)                          % Eigenlijk wil ik dit voor het aantal unique run_labels
    
    row_run_label = find(contains({files(:).name}, run_label{i,:}))
    
    for k = 1:size(row_run_label,2)
        run2sub(k,:) = extractBefore(files(row_run_label(k)).name, '_ses')
    end
    
    %dataBase(i).sub_label = sub_label(i,:);
    
    respLoc = find(contains({files(:).name}, sub_labels{i}));
 
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






% files_all = dir(fullfile(myDataPath.CCEP_allpat));
% files = size(find(contains({files_all(:).name},cfg.sub_labels)),2);         % Find number of rows containing subject info
% diff = size(files_all,1)-files;
% 
% % Find sublabels (necessary because some SPES sessions are ran in two runs)
% for k = diff+1:size(files_all,1)
%     sub_labels{k-diff,:} = extractBefore(files_all(k).name, '_ses');
% end
% 
% % Fill the dataBase with the filename and run-label of the subjects(label)
% dataBase = struct;
% for i=1:size(cfg.sub_labels,2)
%     dataBase(i).sub_label = cfg.sub_labels(i);
%     
%     respLoc = find(contains({files_all(:).name},sub_labels{i}));                % Find which rows contain infor about that patient
% 
%     for k = 1:size(respLoc,2)
%         run_label(k,:) = extractBetween(files_all(i).name, 'clin_', '_CCEP');
%     end
%     
%     % Find run label when SPESclin contains multiple runs
%     if size(run_label,1)>2                          % SPES is ran in more than 1 run
%         for j=1:size(respLoc,2)
%           if contains(files_all(respLoc(j)).name,'10stims') 
%             load(fullfile(files_all(respLoc(j)).folder,files_all(respLoc(j)).name));
%             dataBase(i).ccep10 = ccep10;
%             dataBase(i).filename10 = files_all(respLoc(j)).name;
%             dataBase(i).run_label = run_label(j); 
%           
%          elseif contains(files_all(respLoc(j)).name,'2stims') 
%              load(fullfile(files_all(respLoc(j)).folder,files_all(respLoc(j)).name));
%              dataBase(i).ccep2 = ccep2;   
%              dataBase(i).filename2 = files_all(respLoc(j)).name;
%              dataBase(i).run_label = run_label(j); 
%           end
%         end  
%         
%     else
%         for j=1:size(respLoc,2)
%           if contains(files_all(respLoc(j)).name,'10stims') 
%             load(fullfile(files_all(respLoc(j)).folder,files_all(respLoc(j)).name));
%             dataBase(i).ccep10 = ccep10;
%             dataBase(i).filename10 = files_all(respLoc(j)).name;
%             dataBase(i).run_label = run_label(1);                               % run_label(1) == run_label(2)
%           
%          elseif contains(files_all(respLoc(j)).name,'2stims') 
%              load(fullfile(files_all(respLoc(j)).folder,files_all(respLoc(j)).name));
%              dataBase(i).ccep2 = ccep2;   
%              dataBase(i).filename2 = files_all(respLoc(j)).name;
%              dataBase(i).run_label = run_label(1); 
%           end
%         end 
%     end
%     
%     
%     
%     % ALS SPES IS OPGEDEELD IN MEERDERE RUNS GAAAT HET FOUT!!! DAN WORDT
%     % HET OVERSCHRIJVEN EN IS HET DUS NIET MEER COMPLEET. GELDT NU VORO 703
%     % EN 724  
%   
% end

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
    
    dataBase(subj).statistics = statistical_agreement(myDataPath, dataBase(subj).agreement_parameter, dataBase(subj).ccep10);
end

%% Load electrodes positions (xlsx/electrodes.tsv)
plot_fig = 'n';                 % 'n' when not all ER responses per stim have to be plot, 'y' when you do want to plot all
close all;

for subj = 1:size(dataBase,2)
    visualise_gridstructure(myDataPath, dataBase(subj).ccep10, dataBase(subj).ccep2, dataBase(subj).agreement_parameter,plot_fig);
end
