%% responses to stimulation of subject
% in this script, CCEPs to stimulation are displayed in four figures:
% 1) all trials of a stimulus pair - response electrode combination in
% SPES-clinical, 2) all trials of the same stimulus pair - response
% electrode combination in SPES-propofol, 3) the averaged trial of
% SPES-clinical, 4) the averaged trial of SPES-propofol. 

% 1. set paths
% 2. select subject
% 3. load ccep data
% 4. load visually checked n1-peak data
% 5. figures with all trials
% 6. figures with averaged response.

close all
clear;
clc;

%% general variables

cmap = parula(5);

%% set paths

% add current path from folder which contains this script
rootPath = matlab.desktop.editor.getActiveFilename;
RepoPath = fileparts(rootPath);
matlabFolder = strfind(RepoPath,'matlab');
addpath(genpath(RepoPath(1:matlabFolder+6)));

myDataPath = PRIOS_setLocalDataPath(1);

% housekeeping
clear rootPath RepoPath matlabFolder

%% select subject

sub_label = 'sub-PRIOS09';

%% load data

foldername = fullfile(myDataPath.CCEPpath,'CCEPs');
filenameClin = [sub_label, '_ses-1_task-SPESclin_CCEP'];
filenameProp = [sub_label, '_ses-1_task-SPESprop_CCEP.mat'];

tic
dataBase(1) = load(fullfile(foldername,filenameClin));
dataBase(2) = load(fullfile(foldername,filenameProp));
toc

disp('Data is loaded')

% housekeeping
clear foldername filenameClin filenameProp

%% load visually checked n1-peaks & N1-peak amplitudes

foldername = fullfile(myDataPath.CCEPpath,'checkedN1s');
filenameClin = [sub_label, '_ses-1_task-SPESclin_N1sChecked_comb.mat'];
filenameProp = [sub_label, '_ses-1_task-SPESprop_N1sChecked_comb.mat'];

dataBase(1).ccep = load(fullfile(foldername, filenameClin));
dataBase(2).ccep = load(fullfile(foldername, filenameProp));

% housekeeping
clear foldername filenameClin filenameProp

%% maximal amplitude of N1-peak during SPES-clinical
clc

[respChan,stimChan] = find(dataBase(1).ccep.n1_peak_amplitude_check == min(dataBase(1).ccep.n1_peak_amplitude_check(:)));
tt = dataBase(1).tt;

figure,
plot(tt,squeeze(dataBase(1).cc_epoch_sorted_reref_avg(respChan,stimChan,:)),'Color',cmap(1,:)), 
hold on, 
plot(tt,squeeze(dataBase(2).cc_epoch_sorted_reref_avg(respChan,stimChan,:)),'Color',cmap(2,:)), 
plot(tt(dataBase(1).ccep.n1_peak_sample_check(respChan,stimChan)),dataBase(1).ccep.n1_peak_amplitude_check(respChan,stimChan),...
    'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(1,:))
plot(tt(dataBase(2).ccep.n1_peak_sample_check(respChan,stimChan)),dataBase(2).ccep.n1_peak_amplitude_check(respChan,stimChan),...
    'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:))
hold off

legend('SPES-clinical','SPES-propofol')
xlim([-0.1 0.5])

xlabel('Time (s)')
ylabel('Amplitude (\muV)')
title(sprintf('Maximal N1-amp SPES-clin: %s-%s, %s',dataBase(1).cc_stimchans{stimChan,:},dataBase(1).ch{respChan}))

fprintf('%s-%s, %s: Maximal N1-amp SPES-clin: %3f uV\n',...
    dataBase(1).cc_stimchans{stimChan,:},dataBase(1).ch{respChan},dataBase(1).ccep.n1_peak_amplitude_check(respChan,stimChan))
fprintf('%s-%s, %s: N1-amp SPES-prop: %3f uV\n',...
    dataBase(1).cc_stimchans{stimChan,:},dataBase(1).ch{respChan},dataBase(2).ccep.n1_peak_amplitude_check(respChan,stimChan))

% save figure
figureName = sprintf('%ssupfig4A_maxN1Amp_SPESclin_%s',...
    myDataPath.Figures,sub_label);

set(gcf,'PaperPositionMode','auto')
print(gcf,'-vector','-depsc',figureName)
print(gcf,'-dpng','-r300',figureName)

fprintf('Figure is saved as .eps and .png in \n %s \n',figureName)

% maximal N1 response when stimulating this pair during SPES-propofol
[respChan,stimChan] = find(dataBase(2).ccep.n1_peak_amplitude_check == min(dataBase(2).ccep.n1_peak_amplitude_check(:,stimChan)));
tt = dataBase(1).tt;

figure,
plot(tt,squeeze(dataBase(1).cc_epoch_sorted_reref_avg(respChan,stimChan,:)),'Color',cmap(1,:)), 
hold on, 
plot(tt,squeeze(dataBase(2).cc_epoch_sorted_reref_avg(respChan,stimChan,:)),'Color',cmap(2,:)), 
plot(tt(dataBase(1).ccep.n1_peak_sample_check(respChan,stimChan)),dataBase(1).ccep.n1_peak_amplitude_check(respChan,stimChan),...
    'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(1,:))
plot(tt(dataBase(2).ccep.n1_peak_sample_check(respChan,stimChan)),dataBase(2).ccep.n1_peak_amplitude_check(respChan,stimChan),...
    'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:))
hold off

legend('SPES-clinical','SPES-propofol')
xlim([-0.1 0.5])

xlabel('Time (s)')
ylabel('Amplitude (\muV)')
title(sprintf('Maximal N1-amp SPES-prop same stimpair: %s-%s, %s',dataBase(1).cc_stimchans{stimChan,:},dataBase(1).ch{respChan}))

fprintf('%s-%s, %s: Maximal N1-amp SPES-prop same stimpair: %3f uV\n',...
    dataBase(1).cc_stimchans{stimChan,:},dataBase(1).ch{respChan},dataBase(2).ccep.n1_peak_amplitude_check(respChan,stimChan))
fprintf('%s-%s, %s: N1-amp SPES-clin same stimpair: %3f uV\n',...
    dataBase(1).cc_stimchans{stimChan,:},dataBase(1).ch{respChan},dataBase(1).ccep.n1_peak_amplitude_check(respChan,stimChan))

% save figure
figureName = sprintf('%ssupfig4A_maxN1Amp_SPESclin_SPESprop_%s',...
    myDataPath.Figures,sub_label);

set(gcf,'PaperPositionMode','auto')
print(gcf,'-vector','-depsc',figureName)
print(gcf,'-dpng','-r300',figureName)

fprintf('Figure is saved as .eps and .png in \n %s \n',figureName)

%% maximal amplitude of N1-peak during SPES-propofol

[respChan,stimChan] = find(dataBase(2).ccep.n1_peak_amplitude_check == min(dataBase(2).ccep.n1_peak_amplitude_check(:)));
tt = dataBase(1).tt;

figure,
plot(tt,squeeze(dataBase(1).cc_epoch_sorted_reref_avg(respChan,stimChan,:)),'Color',cmap(1,:)), 
hold on, 
plot(tt,squeeze(dataBase(2).cc_epoch_sorted_reref_avg(respChan,stimChan,:)),'Color',cmap(2,:)), 
plot(tt(dataBase(1).ccep.n1_peak_sample_check(respChan,stimChan)),dataBase(1).ccep.n1_peak_amplitude_check(respChan,stimChan),...
    'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(1,:))
plot(tt(dataBase(2).ccep.n1_peak_sample_check(respChan,stimChan)),dataBase(2).ccep.n1_peak_amplitude_check(respChan,stimChan),...
    'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:))
hold off

legend('SPES-clinical','SPES-propofol')
xlim([-0.1 0.5])

xlabel('Time (s)')
ylabel('Amplitude (\muV)')
title(sprintf('Maximal N1-amp during SPES-prop %s-%s, %s',dataBase(1).cc_stimchans{stimChan,:},dataBase(1).ch{respChan}))

fprintf('%s-%s, %s: Maximal N1-amp SPES-prop: %3f uV\n',...
    dataBase(1).cc_stimchans{stimChan,:},dataBase(1).ch{respChan},dataBase(2).ccep.n1_peak_amplitude_check(respChan,stimChan))
fprintf('%s-%s, %s: N1-amp SPES-clin: %3f uV\n',...
    dataBase(1).cc_stimchans{stimChan,:},dataBase(1).ch{respChan},dataBase(1).ccep.n1_peak_amplitude_check(respChan,stimChan))

% save figure
figureName = sprintf('%ssupfig4A_maxN1Amp_SPESprop_%s',...
    myDataPath.Figures,sub_label);

set(gcf,'PaperPositionMode','auto')
print(gcf,'-vector','-depsc',figureName)
print(gcf,'-dpng','-r300',figureName)

fprintf('Figure is saved as .eps and .png in \n %s \n',figureName)

% maximal N1 response when stimulating this pair during SPES-clinical
[respChan,stimChan] = find(dataBase(1).ccep.n1_peak_amplitude_check == min(dataBase(1).ccep.n1_peak_amplitude_check(:,stimChan)));
tt = dataBase(1).tt;

figure,
plot(tt,squeeze(dataBase(1).cc_epoch_sorted_reref_avg(respChan,stimChan,:)),'Color',cmap(1,:)), 
hold on, 
plot(tt,squeeze(dataBase(2).cc_epoch_sorted_reref_avg(respChan,stimChan,:)),'Color',cmap(2,:)), 
plot(tt(dataBase(1).ccep.n1_peak_sample_check(respChan,stimChan)),dataBase(1).ccep.n1_peak_amplitude_check(respChan,stimChan),...
    'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(1,:))
plot(tt(dataBase(2).ccep.n1_peak_sample_check(respChan,stimChan)),dataBase(2).ccep.n1_peak_amplitude_check(respChan,stimChan),...
    'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:))
hold off

legend('SPES-clinical','SPES-propofol')
xlim([-0.1 0.5])

xlabel('Time (s)')
ylabel('Amplitude (\muV)')
title(sprintf('Maximal N1-amp SPES-clin same stimpair %s-%s, %s',dataBase(1).cc_stimchans{stimChan,:},dataBase(1).ch{respChan}))

fprintf('%s-%s, %s: Maximal N1-amp SPES-clin same stimpair: %3f uV\n',...
    dataBase(1).cc_stimchans{stimChan,:},dataBase(1).ch{respChan},dataBase(1).ccep.n1_peak_amplitude_check(respChan,stimChan))
fprintf('%s-%s, %s: N1-amp SPES-prop same stimpair: %3f uV\n',...
    dataBase(1).cc_stimchans{stimChan,:},dataBase(1).ch{respChan},dataBase(2).ccep.n1_peak_amplitude_check(respChan,stimChan))

% save figure
figureName = sprintf('%ssupfig4A_maxN1Amp_SPESprop_SPESclin_%s',...
    myDataPath.Figures,sub_label);

set(gcf,'PaperPositionMode','auto')
print(gcf,'-vector','-depsc',figureName)
print(gcf,'-dpng','-r300',figureName)

fprintf('Figure is saved as .eps and .png in \n %s \n',figureName)

%% end of script