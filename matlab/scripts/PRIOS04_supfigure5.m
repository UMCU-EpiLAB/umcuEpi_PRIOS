%% effect of ten or two stimulus trials on latency and amplitude
% in this script, CCEPs to stimulation are displayed in four figures:
% 1) all trials of a stimulus pair - response electrode combination in
% SPES-clinical, 2) all trials of the same stimulus pair - response
% electrode combination in SPES-propofol, 3) the averaged trial of
% SPES-clinical, 4) the averaged trial of SPES-propofol. 

% 1. set paths
% 2. select subject
% 3. load ccep data
% 4. load visually checked n1-peak data
% 5. calculate average response with only two stimulus trials
% 5. figures with all trials
% 6. figures with averaged response.

close all
clear;
clc;

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

%% select subject

sub_label = 'sub-PRIOS09';

%% load data

foldername = fullfile(myDataPath.CCEPpath,'CCEPs');
filenameClin = [sub_label, '_ses-1_task-SPESclin_CCEP.mat'];

tic
dataBase(1) = load(fullfile(foldername,filenameClin));
toc

disp('Data is loaded')

% housekeeping
clear foldername filenameClin filenameProp

%% load visually checked n1-peaks

foldername = fullfile(myDataPath.CCEPpath,'checkedN1s');
filenameClin = [sub_label, '_ses-1_task-SPESclin_N1sChecked_comb'];

dataBase(1).ccep = load(fullfile(foldername, filenameClin));

% housekeeping
clear foldername filenameClin filenameProp

%% average the first and sixth stimulus (so first positive (1-2) and first negative (2-1))

cc_epoch_sorted_reref_avg2 = NaN(size(dataBase.cc_epoch_sorted_reref_avg));
tt = dataBase.tt;

for nChan = 1:size(dataBase.cc_epoch_sorted_reref,1)
    for nStim = 1:size(dataBase.cc_epoch_sorted_reref,2)

        signal = squeeze(mean(dataBase.cc_epoch_sorted_reref(nChan,nStim,[1,6],:)));

        tt_preStim = tt>-2 & tt<-.1;

        % subtract median baseline from signal: this is also performed with
        % N1-detection
        median_preStim = median(signal(tt_preStim),'omitnan');
        signal = signal - median_preStim;

        cc_epoch_sorted_reref_avg2(nChan,nStim,:) = signal;
    end
end

dataBase.cc_epoch_sorted_reref_avg2 = cc_epoch_sorted_reref_avg2;

%% check in visually checked N1s the N1-peak latency and amplitude

n1_peak_sample_check = dataBase.ccep.n1_peak_sample_check;

n1_peak_sample = dataBase.ccep.n1_peak_sample;
n1_peak_sample(isnan(n1_peak_sample_check)) = NaN;
n1_peak_sample(abs(n1_peak_sample_check - n1_peak_sample)>0) = NaN;

n1_peak_sample2 = NaN(size(n1_peak_sample_check));
n1_peak_amplitude2 = NaN(size(n1_peak_sample_check));
n1_peak_amplitude = NaN(size(n1_peak_sample_check));

tt = dataBase.tt;

for nChan = 1:size(n1_peak_sample,1)
    for nStim = 1:size(n1_peak_sample,2)
        if ~isnan(n1_peak_sample(nChan,nStim)) 

            [loc,amp] = ccep_peakfinder(squeeze(dataBase.cc_epoch_sorted_reref_avg2(nChan,nStim,tt>0.01 & tt<0.1)),5,[],-1,'false');
            loc = loc+find(tt>0.009,1,'first')+1;
            [minLoc,I] = min(abs(loc-n1_peak_sample(nChan,nStim)));
            
            if minLoc > 30 
%                 figure(1),
%                 plot(squeeze(dataBase.cc_epoch_sorted_reref_avg2(nChan,nStim,:)),'b'),
%                 hold on,
%                 plot(loc(I),amp(I),'*b')
%                 plot(squeeze(dataBase.cc_epoch_sorted_reref_avg(nChan,nStim,:)),'r'),
%                 plot(n1_peak_sample(nChan,nStim), dataBase.cc_epoch_sorted_reref_avg(nChan,nStim,n1_peak_sample(nChan,nStim)),'or')
%                 hold off
%                 pause
            end
            
            n1_peak_amplitude(nChan,nStim) = dataBase.cc_epoch_sorted_reref_avg(nChan,nStim,n1_peak_sample(nChan,nStim));
            n1_peak_sample2(nChan,nStim) = loc(I);
            n1_peak_amplitude2(nChan,nStim) = amp(I);

        end
    end
end

%% statistics
% latency
n1_peak_sample_vec = n1_peak_sample(:);
n1_peak_sample2_vec = n1_peak_sample2(:);

n1_peak_sample2_vec(isnan(n1_peak_sample_vec)) = [];
n1_peak_sample_vec(isnan(n1_peak_sample_vec)) = [];

n1_peak_sample_vec = (n1_peak_sample_vec/2048 - 2) *1000;
n1_peak_sample2_vec = (n1_peak_sample2_vec/2048 - 2) *1000;

lat2 = (mean(n1_peak_sample2_vec));
lat10 = (mean(n1_peak_sample_vec));

plat = signrank(n1_peak_sample2_vec,n1_peak_sample_vec);

fprintf('%s: latency ten trials: %3.1f ms, two trials: %3.1f ms, p = %1.3f \n', ...
    dataBase.sub_label, lat10, lat2, plat)

% amplitude
n1_peak_amplitude_vec = n1_peak_amplitude(:);
n1_peak_amplitude2_vec = n1_peak_amplitude2(:);

n1_peak_amplitude2_vec(isnan(n1_peak_amplitude_vec)) = [];
n1_peak_amplitude_vec(isnan(n1_peak_amplitude_vec)) = [];

amp2 = mean(n1_peak_amplitude2_vec);
amp10 = mean(n1_peak_amplitude_vec);

pamp = signrank(n1_peak_amplitude2_vec,n1_peak_amplitude_vec);

fprintf('%s: amplitude ten trials: %3.0f uV, two trials: %3.0f uV, p = %1.3f \n', ...
    dataBase.sub_label, amp10, amp2, pamp)


%% supplementary figure 5A:
% display the latencies of evoked CCEPs per stimulus pair during SPES-clinical
% for two and ten stimulus trials and connect both with a line.
nSub = find(strcmpi(dataBase.sub_label,{'sub-PRIOS01','sub-PRIOS02','sub-PRIOS03','sub-PRIOS04','sub-PRIOS05','sub-PRIOS09'}));

close all
MkrSze = 10;

h = figure(1);
hold on

scatter((0.8)*ones(size(n1_peak_sample_vec)), ...
    n1_peak_sample_vec, ...
    MkrSze, cmap(nSub,:),'filled')
scatter((1.2)*ones(size(n1_peak_sample2_vec)), ...
    n1_peak_sample2_vec, ...
    MkrSze, cmap(nSub,:),'filled')
for n = 1:size(n1_peak_sample_vec,1)
    plot([0.8 1.2], ...
        [n1_peak_sample_vec(n), ...
        n1_peak_sample2_vec(n)], ...
        'color',cmap(nSub,:),'LineWidth',0.5)
end

plot([0.8, 1.2],...
    [mean(n1_peak_sample_vec), ...
    mean(n1_peak_sample2_vec)], ...
    'k','LineWidth',3)

ylabel('N1 latency (ms)')

h.Units = 'normalized';
h.Position = [0.35 0.5 0.65 0.4];

ax = gca;
ax.XTick = 1;
ax.XTickLabel = {dataBase.sub_label};
ylim([0 100])
xlim([0.6, 1.4])

% save figure
figureName = sprintf('supfig5A_LatencyTwoTenTrials_%s',dataBase.sub_label);
print(h,'-vector','-depsc',fullfile(myDataPath.Figures,figureName))

fprintf('Figure is saved as .eps in \n %s \n', ...
    fullfile(myDataPath.Figures,figureName))

% housekeeping
clear ax figureName h MkrSze n nSub

%% supplementary figure 5B:
% display the amplitudes of evoked CCEPs per stimulus pair during SPES-clinical
% for two and ten stimulus trials and connect both with a line.

nSub = find(strcmpi(dataBase.sub_label,{'sub-PRIOS01','sub-PRIOS02','sub-PRIOS03','sub-PRIOS04','sub-PRIOS05','sub-PRIOS09'}));

MkrSze = 10;

ymin = min([n1_peak_amplitude_vec; n1_peak_amplitude2_vec]);
ymax = max([n1_peak_amplitude_vec; n1_peak_amplitude2_vec]);

h = figure(2);
hold on

scatter((0.8)*ones(size(n1_peak_amplitude_vec)), ...
    n1_peak_amplitude_vec, ...
    MkrSze, cmap(nSub,:),'filled')
scatter((1.2)*ones(size(n1_peak_amplitude2_vec)), ...
    n1_peak_amplitude2_vec, ...
    MkrSze, cmap(nSub,:),'filled')
for n = 1:size(n1_peak_amplitude_vec,1)
    plot([0.8 1.2], ...
        [n1_peak_amplitude_vec(n), ...
        n1_peak_amplitude2_vec(n)], ...
        'color',cmap(nSub,:),'LineWidth',0.5)
end

plot([0.8, 1.2],...
    [mean(n1_peak_amplitude_vec), ...
    mean(n1_peak_amplitude2_vec)], ...
    'k','LineWidth',3)

ylabel('N1 amplitude (uV)')

h.Units = 'normalized';
h.Position = [0.35 0.5 0.65 0.4];

ax = gca;
ax.XTick = 1;
ax.XTickLabel = {dataBase.sub_label};
xlim([0.6, 1.4])

% save figure
figureName = sprintf('supfig5B_AmplitudeTwoTenTrials_%s',dataBase.sub_label);
print(h,'-vector','-depsc',fullfile(myDataPath.Figures,figureName))

fprintf('Figure is saved as .eps in \n %s \n', ...
    fullfile(myDataPath.Figures,figureName))

% housekeeping
clear ax figureName h MkrSze n nSub

%% FDR correction

pFDR = 0.05;
pVals = [0.032 0.015 0.712 0.0001 0.012 0.001 0.506 0.0001 0.085 0.089 0.024 0.271];

[pSort,pInd] = sort(pVals(:));

m = length(pVals);
thisVal = NaN(size(pSort));
for kk = 1:length(pSort)
    thisVal(kk) = (kk/m)*pFDR;
end

pSig = pVals;
pSig(pInd) = pSort < thisVal;
