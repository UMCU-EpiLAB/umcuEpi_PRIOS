% Script PROS01_pipeline_preproces.m should be performed first to obtain the correct documents.

clear; 
% %% Choose patient
ccep_allPat.sub_labels = {'sub-PRIOS01','sub-PRIOS02','sub-PRIOS03','sub-PRIOS04','sub-PRIOS05','sub-PRIOS06','sub-PRIOS09'};
%ccep_allPat.name = {[input('Patient number type (PRIOSXX): ','s')]};

% set paths
ccep_allPat.mode = 'pros';
myDataPath = setLocalDataPath(ccep_allPat);

%% Compare data of two observers and determine the Cohen's Kappa interobserver variability
% Determine this with checked files of two raters/observers

% Exclude responses that were scored differently between observers

dataBase = interobserver_analysis(myDataPath);

% An excel is saved in myDataPath.CCEP_interObVar with the different responces between R1 and R2


%% Find ERs in propofol and not in clinical
for pat = 1:size(dataBase,2)
    
    elec_in_prop = [];
    elec_in_clin = [];
    % Add the detections of clinical and propofol spes.
    % answere == 2 --> both ER, ==1 either one had a ER
    ones = ~isnan(dataBase(pat).ccep_clin.n1_peak_sample) + ~isnan(dataBase(pat).ccep_prop.n1_peak_sample);
    loc_ones = find(ones == 1);
    
    % Determine the value of all those 1's.
    % When non-NAN, than the propofol-SPES had a value and clinical-SPES
    % not
    ones_prop = ~isnan(dataBase(pat).ccep_prop.n1_peak_sample(loc_ones));
    ones_clin = ~isnan(dataBase(pat).ccep_clin.n1_peak_sample(loc_ones));

    elec_prop = loc_ones(ones_prop);
    elec_clin = loc_ones(ones_clin);
    % number of channels for this subject
    nr_ch = size(dataBase(pat).ccep_clin.ch,1);

    for i = 1:size(elec_prop,1)               
        % Determine stimulation pair
        stimp = ceil(elec_prop(i) /nr_ch);
        % Determine electrode
        elec = elec_prop(i) - nr_ch * (stimp-1);
        
%         name_elec_in_prop(i,1) = dataBase(pat).ccep_clin.stimpnames_avg(stimp) ;
%         name_elec_in_prop(i,2) = dataBase(pat).ccep_clin.ch(elec) ;
        elec_in_prop(i,1) = stimp ;
        elec_in_prop(i,2) = elec ;

    end

    dataBase(pat).elec_in_prop = elec_in_prop;

    for i = 1:size(elec_clin,1)
        % Determine stimulation pair
        stimp = ceil(elec_clin(i) /nr_ch);
        % Determine electrode
        elec = elec_clin(i) - nr_ch * (stimp-1);
        
%         name_elec_in_prop(i,1) = dataBase(pat).ccep_clin.stimpnames_avg(stimp) ;
%         name_elec_in_prop(i,2) = dataBase(pat).ccep_clin.ch(elec) ;
        elec_in_clin(i,1) = stimp ;
        elec_in_clin(i,2) = elec ;
        
    end

    dataBase(pat).elec_in_clin = elec_in_clin;


end


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

for n = 1:size(measure,2)        

    for subj = 1:size(dataBase,2)
        
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

T = table(Mult_factor(:,1),Mult_factor(:,2),Mult_factor(:,3),Mult_factor(:,4), 'VariableNames',measure,'RowNames',ccep_allPat.sub_labels);
                
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

close all;

%% Make boxplots of the latency of the N1 peaks.
% Folder Violinplot-Matlab has to be added to the path. 
[table_latency, av_lat_elec, dataBase] = boxplot_N1_peak(dataBase, myDataPath);

prctile([av_lat_elec(:,1);av_lat_elec(:,3);av_lat_elec(:,5);av_lat_elec(:,7);av_lat_elec(:,9);av_lat_elec(:,11);av_lat_elec(:,13)],[25 50 75])
prctile([av_lat_elec(:,2);av_lat_elec(:,4);av_lat_elec(:,6);av_lat_elec(:,8);av_lat_elec(:,10);av_lat_elec(:,12);av_lat_elec(:,14)],[25 50 75])


%% Plot electrodes position on brain with N1-latency information
% Determine N1-latency per brain part
close all

plot_electrodes_on_MRI(myDataPath, table_latency, dataBase, av_lat_elec)


%% Determine the distance between electrodes to determine correlation between N1-latency and electrode distance
close all
distance_elec_stimp(dataBase, myDataPath, av_lat_elec)
