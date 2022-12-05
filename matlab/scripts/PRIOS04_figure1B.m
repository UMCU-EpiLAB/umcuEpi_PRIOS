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

sub_label = 'sub-PRIOS03';

%% load data

foldername = fullfile(myDataPath.CCEPpath,'CCEPs');
filenameClin = [sub_label, '_ses-1_task-SPESclin_CCEP.mat'];
filenameProp = [sub_label, '_ses-1_task-SPESprop_CCEP.mat'];

tic
dataBase(1) = load(fullfile(foldername,filenameClin));
dataBase(2) = load(fullfile(foldername,filenameProp));
toc

disp('Data is loaded')

% housekeeping
clear foldername filenameClin filenameProp

%% load visually checked n1-peaks

foldername = fullfile(myDataPath.CCEPpath,'checkedN1s');
filenameClin = [sub_label, '_ses-1_task-SPESclin_N1sChecked_comb.mat'];
filenameProp = [sub_label, '_ses-1_task-SPESprop_N1sChecked_comb.mat'];

dataBase(1).ccep = load(fullfile(foldername, filenameClin));
dataBase(2).ccep = load(fullfile(foldername, filenameProp));

% housekeeping
clear foldername filenameClin filenameProp

%% find a stimulus pair-response electrode combination in which a CCEP is
% scored in both SPES-clinical and SPES-propofol

[respElecs,stimElecs] = find(~isnan(dataBase(1).ccep.n1_peak_sample_check) == 1 & ...
    ~isnan(dataBase(2).ccep.n1_peak_sample_check) ==1);

%% plot figure: all trials available for the stimulus pair - response electrode combination
% and a grey bar on the part where the stimulus artefact is present

close all

% select one of the stimulus pair - response electrode combinations
nElec = 4;
respElec = respElecs(nElec);
stimElec = stimElecs(nElec);

% check how many trials the stimulus pair is stimulated in SPES-clinical
% and SPES-propofol
numTrials(1) = sum(~isnan(dataBase(1).cc_epoch_sorted_reref(respElec,stimElec,:,1)));
numTrials(2) = sum(~isnan(dataBase(2).cc_epoch_sorted_reref(respElec,stimElec,:,1)));

ySup = -1200;
xmin = - 2;
xmax = 2;
ymax = 0;
ymin = (1+max(numTrials)) * ySup;

for protocol = 1:size(dataBase,2) % for SPES-clin and SPES-prop

    figure(protocol),

    hold on,
    for nTrial = 1:numTrials(protocol)  % for all trials that do not contain NaNs
        plot(dataBase(protocol).tt, ...
            nTrial*ySup + squeeze(dataBase(protocol).cc_epoch_sorted_reref(respElec,stimElec,nTrial,:)),'k')
    end

    patch([0 19/2048 19/2048 0],[ymin ymin ymax ymax], ...
        [211/256, 211/256, 211/256], ...
        'EdgeColor',[211/256, 211/256, 211/256])
    hold off

    xlabel('Time (s)')
    xlim([xmin xmax])
    ylim([ymin  ymax])
    title(sprintf('Stimulus pair: %s-%s, Response channel: %s',...
        dataBase(protocol).cc_stimchans{stimElec,:},dataBase(protocol).ch{respElec}))

    if ~exist(myDataPath.Figures,'dir')
        mkdir(myDataPath.Figures)
    end

    % save figure
    figureName = sprintf('figure1B_%s_%s',sub_label,dataBase(protocol).task_label);
    print(figure(protocol),'-vector','-depsc',fullfile(myDataPath.Figures,figureName))

    fprintf('Figure is saved as .eps in \n %s \n', ...
        fullfile(myDataPath.Figures,figureName))

end

%% plot figure: all trials available for the stimulus pair - response electrode combination
% and the averaged signal, including the N1-latency
% and a grey bar on top of the time window in which the stimulus artefact
% is present.

close all

xmin = - 0.1;
xmax = 0.1;
ymin = -1500;
ymax = 1500;

% for both SPES-clin and SPES-propofol
for protocol = 1:size(dataBase,2)

    figure(protocol),
    hold on,
    for nTrial = 1:size(dataBase(protocol).cc_epoch_sorted_reref,3)  % for all trials
        plot(dataBase(protocol).tt, ...
            squeeze(dataBase(protocol).cc_epoch_sorted_reref(respElec,stimElec,nTrial,:)),':k')
    end

    plot(dataBase(protocol).tt, ...
        squeeze(mean((dataBase(protocol).cc_epoch_sorted_reref(respElec,stimElec,:,:)),3,'omitnan')),'k','LineWidth',3)

    patch([0 19/2048 19/2048 0],[ymin ymin ymax ymax], ...
        [211/256, 211/256, 211/256], ...
        'EdgeColor',[211/256, 211/256, 211/256])

    hold off

    xlabel('Time (s)')
    xlim([xmin xmax])
    ylim([ymin  ymax])
    title(sprintf('Stimulus pair: %s-%s, Response channel: %s',...
        dataBase(protocol).cc_stimchans{stimElec,:},dataBase(protocol).ch{respElec}))

    % save figure
    figureName = sprintf('figure1B_%s_%s_avg',sub_label,dataBase(protocol).task_label);
    print(figure(protocol),'-vector','-depsc',fullfile(myDataPath.Figures,figureName))

    fprintf('Figure is saved as .eps in \n %s \n', ...
        fullfile(myDataPath.Figures,figureName))

    fprintf('N1-peak of %s is found at t = %2.1f ms \n',...
        dataBase(protocol).task_label, ...
        dataBase(protocol).tt(dataBase(protocol).ccep.n1_peak_sample_check(respElec,stimElec))*1000)

end

%% end of script