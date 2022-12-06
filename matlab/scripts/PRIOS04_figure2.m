% compare number of CCEPs to each stimulus pair in propofol and
% clinical SPES

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

%% number of evoked responses per stimulus pair

for nSub = 1:size(dataBase,2)
    % find all non-nans and add this binary matrix in ccepsClin/ccepsProp
    ccepsClin = ~isnan(dataBase(nSub).spesClin.n1_peak_sample_check);
    ccepsProp = ~isnan(dataBase(nSub).spesProp.n1_peak_sample_check);
    dataBase(nSub).ccepsClin = ccepsClin(:);
    dataBase(nSub).ccepsProp = ccepsProp(:);

    % how many cceps occur in both SPES-clin and SPES-prop, or in either
    % one of the protocols --> chi2-test
    [tbl,chi2,p] = crosstab(dataBase(nSub).ccepsClin, ...
        dataBase(nSub).ccepsProp);

    dataBase(nSub).tblClinProp = tbl;
    dataBase(nSub).chi2ClinProp = chi2;
    dataBase(nSub).pClinProp = p;
    
end

% FDR correction
pVals = [dataBase(:).pClinProp];
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

%% figure 2 
% make bar-plot with numbers of CCEPs in both SPES-clin and SPES-prop, 
% or in either of the two

% pre-allocation
tblClinPropComb = NaN(4,size(dataBase,2));

% combine all tables to enable making a bar-plot
for nSub = 1:size(dataBase,2)
    tbl_tmp = dataBase(nSub).tblClinProp;

    tblClinPropComb(:,nSub) = tbl_tmp(:);

end

h = figure(1);
b = bar(tblClinPropComb([4,2,3],:)'); % plot both SPES-clin and SPES-prop; only SPES-clin; only SPES-propofol

b(1).FaceColor = cmap(1,:);
b(2).FaceColor = cmap(3,:);
b(3).FaceColor = cmap(4,:);

% plot numbers on top of each bar
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xtips1 = b(2).XEndPoints;
ytips1 = b(2).YEndPoints;
labels1 = string(b(2).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xtips1 = b(3).XEndPoints;
ytips1 = b(3).YEndPoints;
labels1 = string(b(3).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

legend('in both SPES-protocols','only in SPES-clinical', ...
    'only in SPES-propofol','Location','northeast')
xlabel('Subjects');
ylabel('Number of CCEPs')
ylim([0 1.5*max(max(tblClinPropComb(2:4,:)))])

h.Units = 'normalized';
h.Position = [0.35 0.5 0.65 0.4];
ax = gca;
ax.XTickLabel = {dataBase(:).sub_label};

% save figure
figureName = 'figure2_distributionResponses';
print(h,'-vector','-depsc',fullfile(myDataPath.Figures,figureName))

fprintf('Figure is saved as .eps in \n %s \n', ...
    fullfile(myDataPath.Figures,figureName))

% housekeeping
clear ax b figureName h labels1 nSub tblClinPropComb tbl_tmp xtips1 ytips1

%% check to be sure that the legend in figure 2 is correct
nSub = 1;

numel(find(dataBase(nSub).ccepsClin == 1 & dataBase(nSub).ccepsProp == 1)) % both SPES
numel(find(dataBase(nSub).ccepsClin == 1 & dataBase(nSub).ccepsProp == 0)) % SPES clinical
numel(find(dataBase(nSub).ccepsClin == 0 & dataBase(nSub).ccepsProp == 1)) % SPES propofol

%% end of script