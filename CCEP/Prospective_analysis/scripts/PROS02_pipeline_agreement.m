% Script PROS01_pipeline_preproces.m should be performed first to obtain the correct documents.

clear; 
% %% Choose patient
ccep_allPat.sub_labels = {'sub-PRIOS01','sub-PRIOS02','sub-PRIOS03','sub-PRIOS04','sub-PRIOS05','sub-PRIOS06'};
%ccep_allPat.name = {[input('Patient number type (PRIOSXX): ','s')]};

% set paths
ccep_allPat.mode = 'pros';
myDataPath = setLocalDataPath(ccep_allPat);

%% Compare data of two observers and determine the Cohen's Kappa interobserver variability
% Determine this with checked files of two raters/observers

% Exclude responses that were scored differently between observers

dataBase = interobserver_analysis(myDataPath);

% An excel is saved in myDataPath.CCEP_interObVar with the different responces between R1 and R2

%% Scriptje om eventueel de responses te laten zien die anders gescoord zijn door beide observers
% hiervoor moet dus wel ccep opgeslagen worden.


%% Determine the distance between electrodes to determine correlation between N1-latency and electrode distance







%% determine the agreement between runs
% The determine_agreement function is not only determining the agreement
% when 2 sessions are compared. It could be possible to compare more, but
% then the values for truetrue, truefalse and falsefalse should be changed. 

% load file of continue with variables if they are already in workspace
for subj = 1:size(dataBase,2)
    
    runs(1).ccep = dataBase(subj).ccep_clin;
    runs(2).ccep = dataBase(subj).ccep_prop;
    
    % Determine the agreement between the two matching runs
    dataBase(subj).agreement = determine_agreement(runs, myDataPath);         

    fprintf('%s Overall agreement = %1.2f, positive agreement = %1.2f, negative agreement = %1.2f \n',...
        dataBase(subj).ccep_clin.sub_label, dataBase(subj).agreement.OA, dataBase(subj).agreement.PA, dataBase(subj).agreement.NA)
end


%% Rewrite the adjacency matrices to calculate agreement parameters per electrode instead of per stimulation pair.
% Rewrite the adjacency matrices with stimluation pairs in the columns and the electrodes in the rows
% to an adjacency matrix with electrodes in the columns and rows. 
% This way the network characteristics can be determined per electrode
% instead of per stimulation pair.

close all;

for subj = 1:size(dataBase,2)
    dataBase(subj).ccep_prop = rewrite_Amat(dataBase(subj).ccep_prop, dataBase(subj).agreement.AmatProp);
    dataBase(subj).ccep_clin = rewrite_Amat(dataBase(subj).ccep_clin, dataBase(subj).agreement.AmatClin);
end

%% Determine the indegree, outdegree, Betweenness centrality, the number of ERs per stimpair and the number of ERs per electrode
% The variables are saved in an excel in the run folder of the subject number

for subj = 1:size(dataBase,2)
    dataBase(subj).agreement_parameter = agreement_parameters(dataBase(subj).agreement, ...
        dataBase(subj).ccep_prop, dataBase(subj).ccep_clin);
    
    [dataBase(subj).statistics, dataBase(subj).rank] = statistical_agreement(myDataPath, dataBase(subj).agreement_parameter, dataBase(subj).ccep_clin);
end


%% Visualise the agreement in a scatter plots

scatter_networkPar(dataBase,myDataPath)

%% Determine multiplication factor of the network parameters    
% Data is not normally distributed therefore the median is calculated
measure = {'ERs per stimp','Indegree','Outdegree','BC'};

Mult_factor = zeros(size(dataBase,2), size(measure,2));
        
for subj = 1:size(dataBase,2)

    for n = 1:size(measure,2)
        
        if strcmp(measure{n},'ERs per stimp')
             M_Clin = prctile(dataBase(subj).agreement_parameter.ERs_stimpClin,[25 50 75]);    % prctile treats NaNs as missing values and removes them
             M_Prop = prctile(dataBase(subj).agreement_parameter.ERs_stimpProp,[25 50 75]);       
        
        elseif strcmp(measure{n},'Indegree')
             M_Clin = prctile(dataBase(subj).agreement_parameter.indegreeN_Clin,[25 50 75]);
             M_Prop = prctile(dataBase(subj).agreement_parameter.indegreeN_Prop,[25 50 75]);
        
        elseif strcmp(measure{n},'Outdegree')
             M_Clin = prctile(dataBase(subj).agreement_parameter.outdegreeN_Clin,[25 50 75]);
             M_Prop = prctile(dataBase(subj).agreement_parameter.outdegreeN_Prop,[25 50 75]);
        
        elseif strcmp(measure{n},'BC')
             M_Clin = prctile(dataBase(subj).agreement_parameter.BCN_Clin,[25 50 75]);
             M_Prop = prctile(dataBase(subj).agreement_parameter.BCN_Prop,[25 50 75]);
        end
    
        
        Mult_factor(subj,n) = M_Prop(2)/M_Clin(2);       
    end
end

T = table(Mult_factor(:,1),Mult_factor(:,2),Mult_factor(:,3),Mult_factor(:,4), 'VariableNames',measure,'RowNames',{'PRIOS01','PRIOS02','PRIOS03','PRIOS04','PRIOS05','PRIOS06'});
                
for n=1:size(measure,2)

    Mult = sum(Mult_factor(:,n)) / 6;
    fprintf('Multiplication factor of the %s of the SPES-clin and SPES-prop = %1.1f \n', measure{n}, Mult);

end
    
%% Visualise the network characteristics in the gridstructure of the electrodes positions (xlsx/electrodes.tsv)
% Color scale can be specific for the measurement or the same color
% scale can be used for prop and clin. 
% When using the same: the absolute values can easier be compared
% When using specific colorscale, the ranking is easier comparable.

for subj = 1:size(dataBase,2)
    visualise_gridstructure(myDataPath, dataBase(subj).ccep_clin, dataBase(subj).agreement_parameter);
end

%% Make bar graph of number of ERs per SPES session per patient
% Function used to group/sort all scripts only used for visualisation of
% results for the report

% Also contains heat_map_grid.m & ERs_stimp_violin.m

vis_report(dataBase, myDataPath)

%% Make boxplots of the latency of the N1 peaks.
% Folder Violinplot-Matlab has to be added to the path. 
[table_latency, av_lat_elec] = boxplot_N1_peak(dataBase, myDataPath);



%% Rise and Fall times N1 peak
% vis_P1(myDataPath,dataBase);



%% Plot electrodes position on brain with N1-latency information
% Determine N1-latency per brain part

plot_electrodes_on_MRI(myDataPath, table_latency, dataBase, av_lat_elec)

