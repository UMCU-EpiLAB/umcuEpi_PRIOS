%% comparing the network measures indegree, outdegree and betweenness
% centrality between the two SPES protocols: SPES clinical and SPES
% propofol.

% in this script, Figures4A and 4B are made.

% 1. set paths
% 2. select subjects
% 3. load visually checked N1_peak_latencies
% 4. rewrite adjacency matrix to [respElec x stimElec]
% 5. calculate the indegree, outdegree and betweenness centrality (bc)
% 6. statistical analysis with spearman correlation and FDR correction for
% multiple testing
% 7. make figure 4A bar grahps
% 8. make figure 4B scatter plots

clear
close all
clc

%% general variables

cmap = parula(8);
pFDR = 0.05; % for FDR correction

%% set paths

% add current path from folder which contains this script
rootPath = matlab.desktop.editor.getActiveFilename;
RepoPath = fileparts(rootPath);
matlabFolder = strfind(RepoPath,'matlab');
addpath(genpath(RepoPath(1:matlabFolder+6)));

myDataPath = PRIOS_setLocalDataPath(1);

% housekeeping
clear rootPath RepoPath matlabFolder

%% select subjects

all_sublabels = {'sub-PRIOS01','sub-PRIOS02','sub-PRIOS03','sub-PRIOS04',...
    'sub-PRIOS05','sub-PRIOS09'};

%% load visually checked N1s

clear dataBase

for nSub = 1:size(all_sublabels,2)

    if ~exist(fullfile(myDataPath.CCEPpath,'checkedN1s',...
            [all_sublabels{nSub},'_ses-1_task-SPESclin_N1sChecked_comb.mat']),'file')

        error('You should first run PRIOS01_pipeline_preprocess.m by two observers, PRIOS02_visualCheckN1s.m and PRIOS03_pipeline_agreement.m before you can run this script.')

    else

        tmp = load(fullfile(myDataPath.CCEPpath,'checkedN1s',...
            [all_sublabels{nSub},'_ses-1_task-SPESclin_N1sChecked_comb.mat']));

        dataBase(nSub).sub_label = all_sublabels{nSub};
        dataBase(nSub).dataName = tmp.dataName;
        dataBase(nSub).ch = tmp.ch;
        dataBase(nSub).cc_stimchans = tmp.cc_stimchans;
        dataBase(nSub).cc_stimsets = tmp.cc_stimsets;
        dataBase(nSub).spesClin  = tmp;

    end

    if ~exist(fullfile(myDataPath.CCEPpath,'checkedN1s',...
            [all_sublabels{nSub},'_ses-1_task-SPESprop_N1sChecked_comb.mat']),'file')

        error('You should first run PRIOS01_pipeline_preprocess.m by two observers, PRIOS02_visualCheckN1s.m and PRIOS03_pipeline_agreement.m before you can run this script.')

    else

        tmp = load(fullfile(myDataPath.CCEPpath,'checkedN1s',...
            [all_sublabels{nSub},'_ses-1_task-SPESprop_N1sChecked_comb.mat']));

        dataBase(nSub).sub_label = all_sublabels{nSub};
        dataBase(nSub).dataName = tmp.dataName;
        dataBase(nSub).ch = tmp.ch;
        dataBase(nSub).cc_stimchans = tmp.cc_stimchans;
        dataBase(nSub).cc_stimsets = tmp.cc_stimsets;
        dataBase(nSub).spesProp = tmp;
    end
end

% housekeeping
clear nSub all_sublabels tmp

%% Rewrite the adjacency matrices to calculate agreement parameters per electrode instead of per stimulation pair.
% Rewrite the adjacency matrices with stimluation pairs in the columns and
% the electrodes in the rows to an adjacency matrix with electrodes in the
% columns and rows [respElec x stimElec]. Now we can determine the network
% characteristics per electrode instead of per stimulation pair.

for nSub = 1:size(dataBase,2)
    dataBase(nSub).spesClin.elecAmat = rewrite_Amat(dataBase(nSub).spesClin);
    dataBase(nSub).spesProp.elecAmat = rewrite_Amat(dataBase(nSub).spesProp);
end

% housekeeping
clear nSub

%% Determine the indegree, outdegree, Betweenness centrality

% collect the network measures in a separate structure:
networkMeasures = struct;

for nSub = 1:size(dataBase,2)

    dataBase(nSub).spesClin = calcNetworkMeasures(dataBase(nSub).spesClin);
    dataBase(nSub).spesProp = calcNetworkMeasures(dataBase(nSub).spesProp);

    networkMeasures(nSub).sub_label = dataBase(nSub).sub_label;
    networkMeasures(nSub).indegreeClin = dataBase(nSub).spesClin.indegreeNorm;
    networkMeasures(nSub).indegreeProp = dataBase(nSub).spesProp.indegreeNorm;
    networkMeasures(nSub).outdegreeClin = dataBase(nSub).spesClin.outdegreeNorm;
    networkMeasures(nSub).outdegreeProp = dataBase(nSub).spesProp.outdegreeNorm;
    networkMeasures(nSub).bcClin = dataBase(nSub).spesClin.bcNorm;
    networkMeasures(nSub).bcProp = dataBase(nSub).spesProp.bcNorm;

end

% housekeeping
clear nSub

%% statistical analysis: comparing the network measures (indegree, 
% outdegree, bc) between SPES-clinical and SPES-propofol with spearman
% pairwise correlation and FDR correction for multiple testing

[rho(1),p(1)] = corr(vertcat(networkMeasures(:).indegreeClin), ...
    vertcat(networkMeasures(:).indegreeProp), ...
    'type','Spearman','rows','pairwise');

[rho(2),p(2)] = corr(vertcat(networkMeasures(:).outdegreeClin), ...
    vertcat(networkMeasures(:).outdegreeProp), ...
    'type','Spearman','rows','pairwise');

[rho(3),p(3)] = corr(vertcat(networkMeasures(:).bcClin), ...
    vertcat(networkMeasures(:).bcProp), ...
    'type','Spearman','rows','pairwise');

% FDR correction
[pSort,pInd] = sort(p(:));

m = length(p);
thisVal = NaN(size(pSort));
for kk = 1:length(pSort)
    thisVal(kk) = (kk/m)*pFDR;
end

pSig = p;
pSig(pInd) = pSort < thisVal;

% housekeeping
clear tbl nSub chi2 ccepsClin ccepsProp kk pSort pInd thisVal m

%% Figure 4A: Make bar graph for all patients
% display horizontal bar graphs in which the indegree, outdegree, and bc
% are sorted based on the measures during SPES-clinical. Each subject has
% its own color (defined by cmap).

close all
measure = {'indegree','outdegree','bc'};

% for each measure
for nMeasure = 1:3
    dataClin(1,:) = {networkMeasures(:).([measure{nMeasure},'Clin'])};
    dataProp(1,:) = {networkMeasures(:).([measure{nMeasure},'Prop'])};

    for nSub = 1:size(networkMeasures,2)
        dataClin{2,nSub} = nSub * ones(size(dataClin{1,nSub}));
        dataProp{2,nSub} = nSub * ones(size(dataProp{1,nSub}));
    end

    % combine measures from all separate subjects
    dataClinAll = [-1*vertcat(dataClin{1,:}), vertcat(dataClin{2,:})];
    dataPropAll = [vertcat(dataProp{1,:}), vertcat(dataProp{2,:})];

    % delete NaNs and Infs
    idxClin = isnan(dataClinAll(:,1)) | isinf(dataClinAll(:,1));
    idxProp = isnan(dataPropAll(:,1)) | isinf(dataPropAll(:,1));
    idx = idxClin | idxProp;
    dataClinAll(idx,:) = [];
    dataPropAll(idx,:) = [];

    % sort data based on SPES-clinical
    [~,I] = sort(dataClinAll(:,1),'descend');
    sortDataClinAll = dataClinAll(I,:);
    sortDataPropAll = dataPropAll(I,:);

    h = figure(nMeasure);
    hold on

    for m = 1:size(dataPropAll,1)

        barh(m,sortDataPropAll(m,1),'FaceColor',cmap(sortDataPropAll(m,2),:), ...
            'EdgeColor',cmap(sortDataPropAll(m,2),:),'FaceAlpha',0.2,'EdgeAlpha',0.2)
        barh(m,sortDataClinAll(m,1),'FaceColor',cmap(sortDataClinAll(m,2),:), ...
            'EdgeColor',cmap(sortDataClinAll(m,2),:))

    end
    hold off

    xmin = round(1.1*min(sortDataClinAll(:,1)),1);
    xlim([xmin -1*xmin])

    h.Units = 'normalized';
    h.Position = [0.35 0.4 0.4 0.5];

    % save figure
    figureName = ['figure4A_bar', measure{nMeasure}];
    print(h,'-vector','-depsc',fullfile(myDataPath.Figures,figureName))

    fprintf('Figure is saved as .eps in \n %s \n', ...
        fullfile(myDataPath.Figures,figureName))

    % housekeeping
    clear h dataClin dataClinAll dataProp dataPropAll sortDataPropAll sortDataClinAll nSub I idx figureName idxClin idxProp m xmin

end

% housekeeping
clear nMeasure

%% Figure 4B: network measures SPES-clinical versus SPES-propofol in a scatter plot

close all
szmkr = 15;

% for each measure (indegree, outdegree, bc)
for nMeasure = 1:3

    dataClin = {(networkMeasures(:).([measure{nMeasure},'Clin']))};
    dataClinAll = vertcat(dataClin{:});
    % remove NaNs and Infs
    idxInf = isinf(dataClinAll);
    idxNaN = isnan(dataClinAll);
    dataClinAll = dataClinAll(~idxInf & ~idxNaN);

    dataProp = {(networkMeasures(:).([measure{nMeasure},'Prop']))};
    dataPropAll = vertcat(dataProp{:});
    % remove NaNs and Infs
    idxInf = isinf(dataPropAll);
    idxNaN = isnan(dataPropAll);
    dataPropAll = dataPropAll(~idxInf & ~idxNaN);

    h = figure(nMeasure);
    for nSub = 1:size(dataClin,2)
        hold on
        scatter(dataClin{nSub}, dataProp{nSub},szmkr,cmap(nSub,:),'filled')

    end

    % display a fitted line when the correlation is significant (p <0.05)
    % after FDR correction
    if p(nMeasure) <0.05 && pSig(nMeasure) == 1
        [P,S] = polyfit(dataClinAll,dataPropAll,1);
        [y_fit, ~] = polyval(P,0:0.01:1,S);

        plot(0:0.01:1,y_fit,'k','LineWidth',2)
    end

    hold off
    axis square

    legend({networkMeasures(:).sub_label},'location','bestoutside')

    xlabel(sprintf('Normalized %s SPES-clinical',measure{nMeasure}))
    ylabel(sprintf('Normalized %s SPES-propofol',measure{nMeasure}))

    title(sprintf('rho = %1.3f, p = %1.3f',rho(nMeasure),p(nMeasure)))

    if ~strcmp(measure{nMeasure},'bc')
        xmax = round(max([dataClinAll; dataPropAll]),1);
    else
        xmax = 0.04;
    end

    xlim([0 xmax])
    ylim([0 xmax])

    h.Units = 'normalized';
    h.Position = [0.35 0.4 0.4 0.5];

    % save figure
    figureName = ['figure4B_scatter', measure{nMeasure}];
    print(h,'-vector','-depsc',fullfile(myDataPath.Figures,figureName))

    fprintf('Figure is saved as .eps in \n %s \n', ...
        fullfile(myDataPath.Figures,figureName))

    % housekeeping
    clear h idxInf idxNaN dataClin dataClinAll dataProp dataPropAll P S xmax y_fit figureName

end

% housekeeping
clear nSub nMeasure szmkr

%% end of script