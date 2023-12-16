%% comparing the median N1-peak-amplitude between the SPES protocols 

% in this script, figure 5C is made

% 1. set paths
% 2. select subjects
% 3. load visually checked N1_peak_latencies
% 5. exclude CCEPs that are detected in only one of the SPES-protocols
% 6. calculate wilcoxon signed rank test for the difference between the
%    N1-latency between SPES-clinical and SPES-protocol and apply FDR
%    correction
% 7. print statistics
% 8. make figure 5C

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

        error('You should first run PRIOS01_pipeline_preprocess.m by two observers, PRIOS02_visualCheckN1s.m, PRIOS03_pipeline_agreement.m, and PRIOS03_pipeline_n1peakAmplitude.m before you can run this script.')

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

        error('You should first run PRIOS01_pipeline_preprocess.m by two observers, PRIOS02_visualCheckN1s.m, PRIOS03_pipeline_agreement.m, and PRIOS03_pipeline_n1peakAmplitude.m before you can run this script.')

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

%% include only the N1-peak amplitudes of CCEPs that were detected in both
% find response electrodes in which a CCEP is evoked after stimulation
% a certain stimulus pair in both SPES-clinical and SPES-propofol

for nSub = 1:size(dataBase,2)

    AmatClin = ~isnan(dataBase(nSub).spesClin.n1_peak_amplitude_check);
    AmatProp = ~isnan(dataBase(nSub).spesProp.n1_peak_amplitude_check);

    resp = find(AmatClin == 1 & AmatProp == 1);

    AmatClinAmp = NaN(size(dataBase(nSub).spesClin.n1_peak_amplitude_check));
    AmatPropAmp = NaN(size(dataBase(nSub).spesProp.n1_peak_amplitude_check));

    AmatClinAmp(resp) = dataBase(nSub).spesClin.n1_peak_amplitude_check(resp);
    AmatPropAmp(resp) = dataBase(nSub).spesProp.n1_peak_amplitude_check(resp);

    dataBase(nSub).AmatClinAmp = AmatClinAmp(~isnan(AmatClinAmp));
    dataBase(nSub).AmatPropAmp = AmatPropAmp(~isnan(AmatPropAmp));

    % housekeeping
    clear AmatClin AmatProp resp AmatClinAmp AmatPropAmp

end

% housekeeping
clear nSub

%% calculate difference in latency between the two SPES-protocols 
% and apply FDR correction

for nSub = 1:size(dataBase,2)

    p(nSub) = signrank(dataBase(nSub).AmatClinAmp,dataBase(nSub).AmatPropAmp);

end

% FDR correction
pVals = p;
[pSort,pInd] = sort(pVals(:));

m = length(pVals);
thisVal = NaN(size(pSort));
for kk = 1:length(pSort)
    thisVal(kk) = (kk/m)*pFDR;
end

pSig = pVals;
pSig(pInd) = pSort < thisVal;

% housekeeping
% clear tbl p nSub chi2 ccepsClin ccepsProp kk pSort pInd pVals thisVal m

%% print statistics

clc
for nSub = 1:size(dataBase,2)
    fprintf('%s: SPES-clinical responses/stimulus pair (median, range): %3.1f (%3.1f-%3.1f) \n',...
        dataBase(nSub).sub_label, ...
        round(median(dataBase(nSub).AmatClinAmp)), ...
        min(dataBase(nSub).AmatClinAmp), ...
        max(dataBase(nSub).AmatClinAmp))
    fprintf('%s: SPES-propofol responses/stimulus pair (median, range): %3.1f (%3.1f-%3.1f), p = %2.3f \n\n',...
        dataBase(nSub).sub_label, ...
        round(median(dataBase(nSub).AmatPropAmp)), ...
        min(dataBase(nSub).AmatPropAmp), ...
        max(dataBase(nSub).AmatPropAmp),...
        p(nSub))

end

p_all = signrank(vertcat(dataBase(:).AmatClinAmp),vertcat(dataBase(:).AmatPropAmp));

fprintf('Combined: SPES-clinical responses/stimulus pair (median, range): %3.1f (%3.1f-%3.1f) \n',...
    round(median(vertcat(dataBase(:).AmatClinAmp))), ...
    min(vertcat(dataBase(:).AmatClinAmp)), ...
    max(vertcat(dataBase(:).AmatClinAmp)))

fprintf('Combined: SPES-propofol responses/stimulus pair (median, range): %3.1f (%3.1f-%3.1f), p = %2.3f \n\n',...
    round(median(vertcat(dataBase(:).AmatPropAmp))), ...
    min(vertcat(dataBase(:).AmatPropAmp)), ...
    max(vertcat(dataBase(:).AmatPropAmp)),...
    p_all)

%% supplementary figure 5:
% display the number of evoked CCEPs per stimulus pair during SPES-clinical
% and SPES-propofol for each subjects and connect both with a line.

close all
MkrSze = 10;

for nSub = 1:size(dataBase,2)

    ymin = min([dataBase(nSub).AmatClinAmp; dataBase(nSub).AmatPropAmp]);    
    ymax = max([dataBase(nSub).AmatClinAmp; dataBase(nSub).AmatPropAmp]);
    h = figure(nSub);
    hold on

    scatter((nSub-0.2)*ones(size(dataBase(nSub).AmatClinAmp)), ...
        dataBase(nSub).AmatClinAmp, ...
        MkrSze, cmap(nSub,:),'filled')
    scatter((nSub+0.2)*ones(size(dataBase(nSub).AmatPropAmp)), ...
        dataBase(nSub).AmatPropAmp, ...
        MkrSze, cmap(nSub,:),'filled')
    for n = 1:size(dataBase(nSub).AmatClinAmp,1)
        plot([nSub-0.2 nSub+0.2], ...
            [dataBase(nSub).AmatClinAmp(n), ...
            dataBase(nSub).AmatPropAmp(n)], ...
            'color',cmap(nSub,:),'LineWidth',0.5)
    end

    plot([nSub-0.2, nSub+0.2],...
        [round(median(dataBase(nSub).AmatClinAmp)), ...
        round(median(dataBase(nSub).AmatPropAmp))], ...
        'k','LineWidth',3)

    if p(nSub)<0.05 && pSig(nSub) == 1

        if p(nSub) < 0.001
            text(nSub,ymax+0.7*(abs(ymax)),'***','HorizontalAlignment','center')

        elseif p(nSub) < 0.01
            text(nSub,ymax+0.7*(abs(ymax)),'**','HorizontalAlignment','center')

        elseif p(nSub) < 0.05
            text(nSub,ymax+0.7*(abs(ymax)),'*','HorizontalAlignment','center')

        end
    end

    %     xlabel('Subject')
    ylabel('N1-peak-amplitude (\muV)')

    h.Units = 'normalized';
    h.Position = [0.35 0.5 0.65 0.4];

    ax = gca;
    ax.XTick = 1:6;
    ax.XTickLabel = {dataBase(:).sub_label};
    ylim([1.2*ymin ymax+(abs(ymax))])
    xlim([nSub-0.4, nSub+0.4])

    % save figure
    figureName = sprintf('supfig5_ClinPropAmplitude_%s',dataBase(nSub).sub_label);
    print(h,'-vector','-depsc',fullfile(myDataPath.Figures,figureName))

    fprintf('Figure is saved as .eps in \n %s \n', ...
        fullfile(myDataPath.Figures,figureName))

end

% housekeeping
clear ax figureName h MkrSze n nSub

%% end of script