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
    respLoc = find(contains({files(:).name},ccep_allPat.sub_labels{i}));        
  
     % load all both the SPESclin and SPESprop of the patient
    for j=1:size(respLoc,2)                                                      
       if contains(files(respLoc(j)).name,'clin_filt_check.')               % Change to load the files of interest 
          load(fullfile(files(respLoc(j)).folder,files(respLoc(j)).name));
          dataBase(i).ccep_clin = ccep_clin;
          dataBase(i).filenameClin = files(respLoc(j)).name;
            
       elseif contains(files(respLoc(j)).name,'prop_filt_check.')           % Change to load the files of interest
          load(fullfile(files(respLoc(j)).folder,files(respLoc(j)).name));
          dataBase(i).ccep_prop = ccep_prop;   
          dataBase(i).filenameProp = files(respLoc(j)).name;
       end
    end
end


%% determine the agreement between 2 and 10 stims per run
% The determine_agreement function is not only determining the agreement
% when 2 sessions are compared. It could be possible to compare more, but
% then the values for truetrue, truefalse and falsefalse should be changed. 

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
close all;

for subj = 1:size(dataBase,2)
    dataBase(subj).agreement_parameter = agreement_parameters(dataBase(subj).agreement, ...
        dataBase(subj).ccep_prop, dataBase(subj).ccep_clin);
    
    [dataBase(subj).statistics, dataBase(subj).rank] = statistical_agreement(myDataPath, dataBase(subj).agreement_parameter, dataBase(subj).ccep_clin);
end

%% Create Violinplot of the ranking of the number of ERs per stimulation pair. 

    ERs_perStimp_violin(dataBase,myDataPath)   


%% Visualise the agreement in a scatter plot of the absolute values

    scatter_networkPar(dataBase,myDataPath)

%% Determine multiplication factor of the network parameters    
% Data is not normally distributed therefore the median is calculated
measure = {'ERs per stimp','Indegree','Outdegree','BC'};

Mult_factor = zeros(size(dataBase,2), size(measure,2));
        
for subj = 1:size(dataBase,2)

    for n=1:size(measure,2)
        
        if strcmp(measure{n},'ERs per stimp')
             M_Clin = median(dataBase(subj).agreement_parameter.ERs_stimpClin,'omitnan');
             M_Prop = median(dataBase(subj).agreement_parameter.ERs_stimpProp,'omitnan');       
        elseif strcmp(measure{n},'Indegree')
             M_Clin = median(dataBase(subj).agreement_parameter.indegreeN_Clin,'omitnan');
             M_Prop = median(dataBase(subj).agreement_parameter.indegreeN_Prop,'omitnan');
        elseif strcmp(measure{n},'Outdegree')
             M_Clin = median(dataBase(subj).agreement_parameter.outdegreeN_Clin,'omitnan');
             M_Prop = median(dataBase(subj).agreement_parameter.outdegreeN_Prop,'omitnan');
        elseif strcmp(measure{n},'BC')
             M_Clin = median(dataBase(subj).agreement_parameter.BCN_Clin,'omitnan');
             M_Prop = median(dataBase(subj).agreement_parameter.BCN_Prop,'omitnan');
        end
    
        
        Mult_factor(subj,n) = M_Clin/M_Prop;       
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

%% Visualise the network characteristics in a heatmap to later plot on the MRI
% Create a heatmap of the network characteristics with the outlay of the
% electrodes from the matlabSjabloon in Excel.

for subj = 1:size(dataBase,2)
    
    heat_map_grid(myDataPath, dataBase(subj).ccep_clin, dataBase(subj).agreement_parameter)

end
%% Make bar graph of number of ERs per SPES session per patient
% Function used to group/sort all scripts only used for visualisation of
% results for the report

vis_report(dataBase, myDataPath)


%% Make boxplots of the latency and amplitude of the N1 peaks.
% Folder Violinplot-Matlab has to be added to the path. 
    boxplot_N1_peak(dataBase, myDataPath)

%% Determine the Cohen's Kappa interobserver variability
% Determine this with checked files of two raters/observers
interobserverKappa(myDataPath);