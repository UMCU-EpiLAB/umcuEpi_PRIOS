% Script pipeline_preproces.m should be performed first to obtain the correct documents.

clear; 
% %% Choose patient
ccep_allPat.name = {[input('Patient number type (PRIOSXX): ','s')]};

% set paths
myDataPath = setLocalDataPath(ccep_allPat);

%% Load all ccep files in the folder CCEP_files_allPat
file10 = dir(fullfile(myDataPath.CCEP_allpat,['sub-',ccep_allPat.name{1} '*_CCEP_10*.mat']));
file2 = dir(fullfile(myDataPath.CCEP_allpat,['sub-',ccep_allPat.name{1} '*_CCEP_2*.mat']));

load(fullfile(file10.folder, file10.name));
load(fullfile(file2.folder, file2.name));


%% determine the agreement between 2 and 10 stims per run
% The determine_agreement function is not only determining the agreement
% when 2 sessions are compared. It could be possible to compare more, but
% then the values for W, Z and XandY should be changed. 

% load file of continue with variables if they are already in workspace
if exist('ccep2','var') && exist('ccep10','var') % if you want to continue with the previous part
    runs(1).name = file10.name;
    runs(1).ccep = ccep10;
    
    runs(2).name = file2.name;
    runs(2).ccep = ccep2;
    
else
    fprintf('WARNING: No runs are found');
end

agreement = determine_agreement(runs);          % Deze agreement nog toevoegen aan ccep! handig voor visualize Gridstructure

fprintf('Overall agreement = %1.2f, positive agreement = %1.2f, negative agreement = %1.2f \n',...
    agreement.agreement_run.OA, agreement.agreement_run.PA, agreement.agreement_run.NA)

% fprintf('Overall agreement = %1.2f, positive agreement = %1.2f, negative agreement = %1.2f \n',...
%     agreement_stim.OA, agreement_stim.PA, agreement_stim.NA)

%% Determine the location of the ones (ER vs. No-ER)
% ccep10 is only necessary for the channels and stimpairs and those are
% equal for 2 and 10 stimuli so does not matter which database is used.
LocOnes = find_ones(ccep10, agreement.agreement_run);
 

%% Plot all 10 stimuli and the average for the 10 stims and the 2 stims
% This does not work without epoch_sorted information, though, matlab
% cannot save since it is to big.
plot_fig = input('Do you want plot the 10 stimuli and the average signals? [y/n] ','s');

if strcmp(plot_fig,'y')
    ccep10.save_fig = str2double(input('Do you want to save the figures? [yes = 1, no = 0]: ','s'));
    plot_all_ccep_and_av(ccep10, ccep2, myDataPath, LocOnes, agreement);
end

%% Calculate agreement parameters
close all;
ccep2 = rewrite_Amat(ccep2, agreement.Amat2);
ccep10 = rewrite_Amat(ccep10, agreement.Amat10);

%% Determine the indegree, outdegree, Betweenness centrality, the number of ERs per stimpair and the number of ERs per electrode
% The variables are saved in an excel in the run folder of the subject number
close all;

agreement_parameter = agreement_parameters(agreement, ccep2, ccep10, myDataPath);


statistics = statistical_agreement(agreement_parameter,ccep10);
%% load electrodes positions (xlsx/electrodes.tsv)
% database (ccep10) is only used for channels and stimpairs and these are
% equal for 2 and 10, so does not matter which database is used.
visualise_gridstructure(myDataPath, ccep10, ccep2, agreement_parameter);



