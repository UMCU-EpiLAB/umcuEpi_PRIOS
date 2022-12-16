%% displaying the mean latency and SEM for each subject in SPES propofol 
% and SPES clinical

% in this script, Figure 5A is made. 

% 1. set paths
% 2. select subjects
% 3. load CCEP data
% 4. load visually checked N1_peak_latencies
% 5. exclude CCEPs that are detected in only one of the SPES-protocols
% 6. calculate mean response and standard error of the mean for each subject
% 7. make figure 5A

clear
close all
clc

%% general variables

cmap = parula(8);

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

%% load data

foldername = fullfile(myDataPath.CCEPpath,'CCEPs');
dataBase = struct;

for nSub = 1:size(all_sublabels,2)
    filenameClin = [all_sublabels{nSub}, '_ses-1_task-SPESclin_CCEP.mat'];
    filenameProp = [all_sublabels{nSub}, '_ses-1_task-SPESprop_CCEP.mat'];

    tic
    tmpClin = load(fullfile(foldername,filenameClin));
    tmpProp = load(fullfile(foldername,filenameProp));
    toc

    fprintf('Data of %s is loaded\n',all_sublabels{nSub})

    dataBase(nSub).sub_label = tmpClin.sub_label;
    dataBase(nSub).ses_label = tmpClin.ses_label;
    dataBase(nSub).ch = tmpClin.ch;
    dataBase(nSub).cc_stimchans = tmpClin.cc_stimchans;
    dataBase(nSub).cc_stimsets = tmpClin.cc_stimsets;
    dataBase(nSub).spesClin = tmpClin;
    dataBase(nSub).spesProp = tmpProp;

    % housekeeping
    clear filenameClin filenameProp tmpClin tmpProp

end

disp('All data is loaded')

% housekeeping
clear foldername nSub

%% load visually checked n1-peaks

foldername = fullfile(myDataPath.CCEPpath,'checkedN1s');

for nSub = 1:size(all_sublabels,2)

    filenameClin = [all_sublabels{nSub}, '_ses-1_task-SPESclin_N1sChecked_comb.mat'];
    filenameProp = [all_sublabels{nSub}, '_ses-1_task-SPESprop_N1sChecked_comb.mat'];

    tic
    dataBase(nSub).spesClin.ccep = load(fullfile(foldername, filenameClin));
    dataBase(nSub).spesProp.ccep = load(fullfile(foldername, filenameProp));
    toc

    % housekeeping
    clear filenameClin filenameProp

    fprintf('Checked n1s of %s are loaded\n',all_sublabels{nSub})
end

% housekeeping
clear foldername nSub all_sublabels

%% find response electrodes in which a CCEP is evoked after stimulation
% a certain stimulus pair in both SPES-clinical and SPES-propofol

for nSub = 1:size(dataBase,2)

    AmatClin = ~isnan(dataBase(nSub).spesClin.ccep.n1_peak_sample_check);
    AmatProp = ~isnan(dataBase(nSub).spesProp.ccep.n1_peak_sample_check);

    [respCh,stimPair] = find(AmatClin == 1 & AmatProp == 1);

    ccepClin = NaN(size(dataBase(nSub).spesClin.cc_epoch_sorted_reref_avg));
    ccepProp = NaN(size(dataBase(nSub).spesClin.cc_epoch_sorted_reref_avg));
    AmatClinLat = NaN(size(dataBase(nSub).spesClin.ccep.n1_peak_sample_check));
    AmatPropLat = NaN(size(dataBase(nSub).spesProp.ccep.n1_peak_sample_check));

    for n = 1:size(respCh,1)
        ccepClin(respCh(n),stimPair(n),:) = dataBase(nSub).spesClin.cc_epoch_sorted_reref_avg(respCh(n),stimPair(n),:);
        ccepProp(respCh(n),stimPair(n),:) = dataBase(nSub).spesProp.cc_epoch_sorted_reref_avg(respCh(n),stimPair(n),:);
        AmatClinLat(respCh(n),stimPair(n)) = dataBase(nSub).spesClin.ccep.n1_peak_sample_check(respCh(n),stimPair(n));
        AmatPropLat(respCh(n),stimPair(n)) = dataBase(nSub).spesProp.ccep.n1_peak_sample_check(respCh(n),stimPair(n));

    end

    dataBase(nSub).ccepClin = ccepClin;
    dataBase(nSub).ccepProp = ccepProp;
    dataBase(nSub).AmatClinLat = AmatClinLat;
    dataBase(nSub).AmatPropLat = AmatPropLat;

    % housekeeping
    clear AmatClin AmatProp respCh stimPair ccepClin ccepProp AmatClinLat AmatPropLat n

end

% housekeeping
clear nSub

%% calculate mean and SEM of CCEP-signal and N1-latency

% pre-allocation
meanccepClin = NaN(size(dataBase,2),size(dataBase(1).ccepClin,3));
upperClin = NaN(size(dataBase,2),size(dataBase(1).ccepClin,3));
lowerClin = NaN(size(dataBase,2),size(dataBase(1).ccepClin,3));
meanccepProp = NaN(size(dataBase,2),size(dataBase(1).ccepProp,3));
upperProp = NaN(size(dataBase,2),size(dataBase(1).ccepProp,3));
lowerProp = NaN(size(dataBase,2),size(dataBase(1).ccepProp,3));

meanClinLat = NaN(size(dataBase,2),1);
upperClinLat = NaN(size(dataBase,2),1);
lowerClinLat = NaN(size(dataBase,2),1);
meanPropLat = NaN(size(dataBase,2),1);
upperPropLat = NaN(size(dataBase,2),1);
lowerPropLat = NaN(size(dataBase,2),1);

for nSub = 1:size(dataBase,2)

    ccepClin = dataBase(nSub).ccepClin;
    ccepClin = reshape(ccepClin,size(ccepClin,1)*size(ccepClin,2),size(ccepClin,3));
    ccepClin = ccepClin(~isnan(ccepClin(:,1)),:);
    ccepProp = dataBase(nSub).ccepProp;
    ccepProp = reshape(ccepProp,size(ccepProp,1)*size(ccepProp,2),size(ccepProp,3));
    ccepProp = ccepProp(~isnan(ccepProp(:,1)),:);
    clinLat = dataBase(nSub).AmatClinLat(:);
    propLat = dataBase(nSub).AmatPropLat(:);

    n = size(ccepClin,1);

    % calculate mean and standard error of the mean of SPES-clin
    meanccepClin(nSub,:) = mean(ccepClin,1);
    SEMccepClin = std(ccepClin,0,1)/sqrt(n);
    upperClin(nSub,:) = meanccepClin(nSub,:) + SEMccepClin;
    lowerClin(nSub,:) = meanccepClin(nSub,:) - SEMccepClin;

    % calculate mean and standard error of the mean of SPES-prop
    meanccepProp(nSub,:) = mean(ccepProp,1);
    SEMccepProp = std(ccepProp,0,1)/sqrt(n);
    upperProp(nSub,:) = meanccepProp(nSub,:) + SEMccepProp;
    lowerProp(nSub,:) = meanccepProp(nSub,:) - SEMccepProp;

    meanClinLat(nSub,:) = round(mean(clinLat,'omitnan'));
    SEMClinLat = round(std(clinLat,'omitnan'));
    upperClinLat(nSub,:) = meanClinLat(nSub,:) + SEMClinLat;
    lowerClinLat(nSub,:) = meanClinLat(nSub,:) - SEMClinLat;

    meanPropLat(nSub,:) = round(mean(propLat,'omitnan'));
    SEMPropLat = round(std(propLat,'omitnan'));
    upperPropLat(nSub,:) = meanPropLat(nSub,:) + SEMPropLat;
    lowerPropLat(nSub,:) = meanPropLat(nSub,:) - SEMPropLat;

    % housekeeping
    clear n ccepClin ccepProp clinLat propLat SEMccepClin SEMccepProp SEMClinLat SEMPropLat

end

% housekeeping
clear nSub

%% figure 5A

ySub = -1800;
yPropSup = -500;
ySEMclin = -1400;
ySEMprop = -1500;

tt = dataBase(1).spesClin.tt;

kk = figure;

for nSub = 1:size(dataBase,2)
  
    % plot SEM and mean of SPES-clinical and SPES-propofol
    fill([tt, fliplr(tt)],...
        (nSub * ySub) + [upperClin(nSub,:), fliplr(lowerClin(nSub,:))], ...
        cmap(nSub,:),'EdgeColor',cmap(nSub,:), ...
        'FaceAlpha',0.7,'EdgeAlpha',0.7)
    hold on
    h(nSub) = plot(tt, ...
        (nSub * ySub) + meanccepClin(nSub,:), ...
        'Color',cmap(nSub,:),'LineWidth',1);

    fill([tt, fliplr(tt)],...
        (nSub * ySub) + yPropSup + [upperProp(nSub,:), fliplr(lowerProp(nSub,:))], ...
        cmap(nSub,:),'EdgeColor',cmap(nSub,:), ...
        'FaceAlpha',0.3,'EdgeAlpha',0.3)
    plot(tt, ...
        (nSub * ySub) + yPropSup + meanccepProp(nSub,:), ...
        'Color',cmap(nSub,:),'LineWidth',1,'LineStyle',':')

    fill(tt([lowerClinLat(nSub) upperClinLat(nSub) upperClinLat(nSub) lowerClinLat(nSub)]),...
        (nSub * ySub) + ySEMclin + [1 1 100 100], ...
        cmap(nSub,:),'EdgeColor',cmap(nSub,:), ...
        'FaceAlpha',0.7,'EdgeAlpha',0.7)

    fill(tt([lowerPropLat(nSub) upperPropLat(nSub) upperPropLat(nSub) lowerPropLat(nSub)]),...
        (nSub * ySub) + ySEMprop + [1 1 100 100], ...
        cmap(nSub,:),'EdgeColor',cmap(nSub,:), ...
        'FaceAlpha',0.3,'EdgeAlpha',0.3)
    
end

fill([0 19/2048 19/2048 0],[-1500 -1500 -12400 -12400], ...
    [192/256 192/256 192/256],'FaceAlpha',0.2,'EdgeColor','none')

hold off

legend(h,{dataBase(:).sub_label},'location','bestoutside')
xlabel('Time (s)')
xlim([-0.05 0.1])
ylim([-12400 -1500])

kk.Units = 'normalized';
kk.Position = [0.3 0.1 0.35 0.8];

% save figure
figureName = 'figure5A_meanCcepClinProp';
print(kk,'-vector','-depsc',fullfile(myDataPath.Figures,figureName))

fprintf('Figure is saved as .eps in \n %s \n', ...
    fullfile(myDataPath.Figures,figureName))

%% end of script