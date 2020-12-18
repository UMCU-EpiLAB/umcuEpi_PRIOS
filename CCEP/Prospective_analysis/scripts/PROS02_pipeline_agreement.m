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

%% Determine multiplication factor of the network parameters    
% Data is not normally distributed therefore the median is calculated

for subj = 1:size(dataBase,2)
   measure = {'ERs per stimp','Indegree','Outdegree','BC'};

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
    
%% load electrodes positions (xlsx/electrodes.tsv)
% database (ccep10) is only used for channels and stimpairs and these are
% equal for 2 and 10, so does not matter which database is used.

for subj = 1:size(dataBase,2)
    visualise_gridstructure(myDataPath, dataBase(subj).ccep_clin, dataBase(subj).ccep_prop, dataBase(subj).agreement_parameter);
end

%% Make bar graph of number of ERs per SPES session per patient
figure('Position',[407,689,939,373])
ax1 = axes('Position',[0.074,0.11,0.9,0.82]);

for subj = 1:size(dataBase,2)
   
    % prealloction of the column number   
    clin_colm = 2*subj-1;                      
    prop_colm = 2*subj; 

    ERs_tot(subj,clin_colm) = sum(sum(~isnan(dataBase(subj).ccep_clin.n1_peak_amplitude_check)));
    ERs_tot(subj,prop_colm) = sum(sum(~isnan(dataBase(subj).ccep_prop.n1_peak_amplitude_check)));
        
    priosLab = extractAfter(dataBase(subj).sub_label,'sub-');
        
end

x = {'PRIOS01','PRIOS02','PRIOS03','PRIOS04','PRIOS05','PRIOS06'};
 
X = categorical(x);
b =   bar(ax1,X,ERs_tot,1);         % bar(ax1,X,ERs_tot,1);

for i = 1:2:size(ERs_tot,2)
        b([i]).FaceColor(:) =  [0.5843 0.8157 0.9882];          %[0 0 1]
        b([i+1]).FaceColor(:) = [0.9882 0.6157 0.5843];                        %[1 0 0]
end
 

% Place the Number of ERs next to the column
for i = 1:2:12
    xtips1 = b(i).XEndPoints;
    ytips1 = b(i).YEndPoints;
    labels1 = string(b(i).YData);
    Zero = find(labels1 == '0');
    labels1(Zero) = NaN;
    text(xtips1,ytips1,labels1,'HorizontalAlignment','right',...
    'VerticalAlignment','bottom','FontSize',11,'fontweight','Bold') 
end

for i = 2:2:12
    xtips1 = b(i).XEndPoints;
    ytips1 = b(i).YEndPoints;
    labels1 = string(b(i).YData);
    Zero = find(labels1 == '0');
    labels1(Zero) = NaN;
    text(xtips1,ytips1,labels1,'HorizontalAlignment','left',...
    'VerticalAlignment','bottom','FontSize',11,'fontweight','Bold') 
end
legend('SPES-clin','SPES-Prop')
ylabel('Number of ERs');
title('Total number of ERs evoked per SPES session')

% Save figure
outlabel='ERs_per_stimp.jpg';
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'jpg')


%% Make boxplots of the latency and amplitude of the N1 peaks.
% Folder Violinplot-Matlab has to be added to the path. 
    boxplot_N1_peak(dataBase, myDataPath)


%% Determine the Cohen's Kappa interobserver variability
% Determine this with checked files of two raters/observers
interobserverKappa(myDataPath);
