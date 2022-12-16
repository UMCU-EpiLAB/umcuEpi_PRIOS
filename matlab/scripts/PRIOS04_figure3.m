%% compare number of CCEPs in propofol and clinical SPES
% In this script, figure 3 of the manuscript is made. 

% 1. set paths
% 2. select subjects
% 3. load visually checked N1s
% 4. calculate the number of response channels with evoked CCEPs for each
% stimulus pair during SPES-clinical and SPES-propofol and analyse
% statistics with a Wilcoxon-signed-rank test and perform FDR correction. 
% 5. print the statistics in command window
% 6. make figure showing the number of evoked responses per stimulus pair
% for each subject

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

    % load visually checked N1s of SPES-clinical
    if ~exist(fullfile(myDataPath.CCEPpath,'checkedN1s',...
            [all_sublabels{nSub},'_ses-1_task-SPESclin_N1sChecked_comb.mat']),'file')

        error('You should first run PRIOS01_pipeline_preprocess.m by two observers, PRIOS02_visualCheckN1s.m and PRIOS03_pipeline_agreement.m before you can run this script.')

    else

        tmp = load(fullfile(myDataPath.CCEPpath,'checkedN1s',...
            [all_sublabels{nSub},'_ses-1_task-SPESclin_N1sChecked_comb.mat']));
       
        dataBase(nSub).sub_label = all_sublabels{nSub};
        dataBase(nSub).cc_stimchans = tmp.cc_stimchans;
        dataBase(nSub).cc_stimsets = tmp.cc_stimsets;
        dataBase(nSub).spesClin  = tmp;

    end

    % load visually checked N1s of SPES-propofol
    if ~exist(fullfile(myDataPath.CCEPpath,'checkedN1s',...
            [all_sublabels{nSub},'_ses-1_task-SPESprop_N1sChecked_comb.mat']),'file')

        error('You should first run PRIOS01_pipeline_preprocess.m by two observers, PRIOS02_visualCheckN1s.m and PRIOS03_pipeline_agreement.m before you can run this script.')

    else

        tmp = load(fullfile(myDataPath.CCEPpath,'checkedN1s',...
            [all_sublabels{nSub},'_ses-1_task-SPESprop_N1sChecked_comb.mat']));

        dataBase(nSub).sub_label = all_sublabels{nSub};
        dataBase(nSub).cc_stimchans = tmp.cc_stimchans;
        dataBase(nSub).cc_stimsets = tmp.cc_stimsets;
        dataBase(nSub).spesProp = tmp;
    end
end

% housekeeping
clear nSub all_sublabels tmp

%% calculate the number of evoked responses per stimulus pair
% for SPES-clinical and SPES-propofol, analyse statistics with
% Wilcoxon-signed-rank test and apply FDR correction

for nSub = 1:size(dataBase,2)
    ccepsClin = ~isnan(dataBase(nSub).spesClin.n1_peak_sample_check);
    ccepsProp = ~isnan(dataBase(nSub).spesProp.n1_peak_sample_check);

    % number of evoked responses per stimulus pair
    dataBase(nSub).numRespClin = sum(ccepsClin,1);
    dataBase(nSub).numRespProp = sum(ccepsProp,1);

    p = signrank(dataBase(nSub).numRespClin, dataBase(nSub).numRespProp);

    dataBase(nSub).p_numResp = p;
    
end

% FDR correction
pVals = [dataBase(:).p_numResp];
[pSort,pInd] = sort(pVals(:));

m = length(pVals);
thisVal = NaN(size(pSort));
for kk = 1:length(pSort)
    thisVal(kk) = (kk/m)*pFDR;
end

pSig = pVals;
pSig(pInd) = pSort < thisVal;

% housekeeping
clear tbl p nSub chi2 ccepsClin ccepsProp kk pSort pInd pVals thisVal m

%% print statistics

clc
for nSub = 1:size(dataBase,2)
    fprintf('%s: SPES-clinical responses/stimulus pair (median, range): %2.1f (%2.1f-%2.1f) \n',...
        dataBase(nSub).sub_label, median(dataBase(nSub).numRespClin), ...
        min(dataBase(nSub).numRespClin), max(dataBase(nSub).numRespClin))
fprintf('%s: SPES-propofol responses/stimulus pair (median, range): %2.1f (%2.1f-%2.1f), p = %2.3f \n\n',...
        dataBase(nSub).sub_label, median(dataBase(nSub).numRespProp), ...
        min(dataBase(nSub).numRespProp), max(dataBase(nSub).numRespProp),...
        dataBase(nSub).p_numResp)

end

fprintf('Combined: SPES-clinical responses/stimulus pair (median, range): %2.1f (%2.1f-%2.1f) \n',...
    median([dataBase(:).numRespClin]), ...
        min([dataBase(:).numRespClin]), max([dataBase(:).numRespClin]))

fprintf('Combined: SPES-propofol responses/stimulus pair (median, range): %2.1f (%2.1f-%2.1f) \n\n',...
    median([dataBase(:).numRespProp]), ...
        min([dataBase(:).numRespProp]), max([dataBase(:).numRespProp]))

%% figure 3:
% display the number of evoked CCEPs per stimulus pair during SPES-clinical
% and SPES-propofol for each subjects and connect both with a line. 

close all
MkrSze = 10;

h = figure(1);
hold on
for nSub = 1:size(dataBase,2)

    scatter((nSub-0.2)*ones(size(dataBase(nSub).numRespClin)), ...
        dataBase(nSub).numRespClin, MkrSze, ...
        cmap(nSub,:),'filled')
    scatter((nSub+0.2)*ones(size(dataBase(nSub).numRespProp)), ...
        dataBase(nSub).numRespProp, MkrSze, ...
        cmap(nSub,:),'filled')
    for n = 1:size(dataBase(nSub).numRespClin,2)
        plot([nSub-0.2 nSub+0.2], ...
            [dataBase(nSub).numRespClin(n), dataBase(nSub).numRespProp(n)], ...
            'color',cmap(nSub,:),'LineWidth',0.5)
    end

    plot([nSub-0.2, nSub+0.2],...
        [median(dataBase(nSub).numRespClin), median(dataBase(nSub).numRespProp)],'k','LineWidth',2)

end

xlabel('Subjects')
ylabel('Number of CCEPs per stimulus pair')

h.Units = 'normalized';
h.Position = [0.35 0.5 0.65 0.4];

ax = gca;
ax.XTick = 1:6;
ax.XTickLabel = {dataBase(:).sub_label};
ylim([0 1.5*max([[dataBase(:).numRespClin], [dataBase(:).numRespProp]])])

% save figure
figureName = 'figure3_numCCEPs';
print(h,'-vector','-depsc',fullfile(myDataPath.Figures,figureName))

fprintf('Figure is saved as .eps in \n %s \n', ...
    fullfile(myDataPath.Figures,figureName))

% housekeeping
clear ax figureName h MkrSze n nSub

%% end of script