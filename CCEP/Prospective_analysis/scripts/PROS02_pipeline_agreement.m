% Script PROS01_pipeline_preproces.m should be performed first to obtain the correct documents.

clear; 

%% Choose patient
ccep_allPat.sub_labels = {'sub-PRIOS01','sub-PRIOS02','sub-PRIOS03','sub-PRIOS04','sub-PRIOS05','sub-PRIOS06','sub-PRIOS09'};

% set paths
ccep_allPat.mode = 'pros';
myDataPath = setLocalDataPath(ccep_allPat);

%% Compare data of two observers and determine the Cohen's Kappa interobserver variability
% Determine this with checked files of two raters/observers

% Exclude responses that were scored differently between observers

dataBase = interobserver_analysis(myDataPath);

% An excel is saved in myDataPath.CCEP_interObVar with the different responces between R1 and R2


%% Make figure of average response of all responses with N1-peak
% fig_average_resp(dataBase, myDataPath)


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
        elec_in_prop(i,1) = stimp ;
        elec_in_prop(i,2) = elec ;

    end

    dataBase(pat).elec_in_prop = elec_in_prop;

    for i = 1:size(elec_clin,1)
        % Determine stimulation pair
        stimp = ceil(elec_clin(i) /nr_ch);
        % Determine electrode
        elec = elec_clin(i) - nr_ch * (stimp-1);
        elec_in_clin(i,1) = stimp ;
        elec_in_clin(i,2) = elec ;
        
    end

    dataBase(pat).elec_in_clin = elec_in_clin;

end


%% Create Adjacency matric of both sessions
% 
for subj = 1:size(dataBase,2)
    
    runs(1).ccep = dataBase(subj).ccep_clin;
    runs(2).ccep = dataBase(subj).ccep_prop;
    
    % Determine the agreement between the two matching runs
    dataBase(subj).Amat = determine_aMat(runs);         
end


%% Rewrite the adjacency matrices to calculate agreement parameters per electrode instead of per stimulation pair.
% Rewrite the adjacency matrices with stimluation pairs in the columns and the electrodes in the rows
% to an adjacency matrix with electrodes in the columns and rows. 
% This way the network characteristics can be determined per electrode
% instead of per stimulation pair.

close all;

for subj = 1:size(dataBase,2)
    dataBase(subj).ccep_prop = rewrite_Amat(dataBase(subj).ccep_prop, dataBase(subj).Amat.AmatProp);
    dataBase(subj).ccep_clin = rewrite_Amat(dataBase(subj).ccep_clin, dataBase(subj).Amat.AmatClin);
end

%% Determine the indegree, outdegree, Betweenness centrality, the number of ERs per stimpair and the number of ERs per electrode
% The variables are saved in an excel in the run folder of the subject number
clc

for subj = 1:size(dataBase,2)
    dataBase(subj).agreement_parameter = agreement_parameters(dataBase(subj).Amat, ...
        dataBase(subj).ccep_prop, dataBase(subj).ccep_clin);
    
    [dataBase(subj).statistics, dataBase(subj).rank] = statistical_agreement(dataBase(subj).agreement_parameter, dataBase(subj).ccep_clin);
end

%% Visualise the agreement in a scatter plots

scatter_networkPar(dataBase,myDataPath)

%% Make bar graph for all patients

barGraphStims(dataBase,myDataPath)


%% Determine multiplication factor of the network parameters    
% Data are not normally distributed therefore the median is calculated
clc

multiplication_fac(dataBase, ccep_allPat)

    
%% Make bar graph of number of ERs per SPES session per patient
% Function used to group/sort all scripts only used for visualisation of
% results for the report
% Also contains ERs_stimp_violin.m

vis_report(dataBase, myDataPath)

close all;

%% Make boxplots of the latency of the N1 peaks.
% Folder Violinplot-Matlab has to be added to the path. 
[dataBase] = boxplot_N1_peak(dataBase, myDataPath);

%% Plot electrodes position on brain with N1-latency information
% Determine N1-latency per brain part
close all

plot_electrodes_on_MRI(myDataPath, dataBase)


%% Determine the distance between electrodes to determine correlation between N1-latency and electrode distance
close all
distance_elec_stimp(dataBase, myDataPath, av_lat_elec)



