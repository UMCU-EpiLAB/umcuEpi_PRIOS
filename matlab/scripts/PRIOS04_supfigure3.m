%% plotting a figure with the brain and electrodes on the brain
% in this script, the electrodes in mni305 space are plotted on a
% freesurfer brain. This figure is also displayed in the article, in Figure
% 1A.

% 1. set paths
% 2. select subject
% 3. add electrodes.tsv (see BIDS structure)
% 4. calculate indegree, outdegree
% 5. load mni305 pial
% 6. plot subject and electrodes on the freesurfer brain

close all
clear
clc

%% general variables

cmap = parula(5);
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

%% select subject

sub_label = 'sub-PRIOS09';

%% add electrodes.tsv

elecsName = fullfile(myDataPath.dataPath, sub_label,'ses-1','ieeg',...
    [sub_label,'_ses-1_electrodes.tsv']);

tb_electrodes = readtable(elecsName,'FileType','text','Delimiter','\t');
idx_include = ~strcmp(tb_electrodes.group,'other');

tb_electrodes = tb_electrodes(idx_include,:);

ch = tb_electrodes.name;

% convert electrode positions into a matrix
if iscell(tb_electrodes.x)

    elecmatrix = NaN(size(tb_electrodes,1),3);
    for ll = 1:size(tb_electrodes,1)
        if ~isequal(tb_electrodes.x{ll},'n/a')
            elecmatrix(ll,:) = [str2double(tb_electrodes.x{ll}), ...
                str2double(tb_electrodes.y{ll}),...
                str2double(tb_electrodes.z{ll})];
        end
    end
else

    elecmatrix = [tb_electrodes.x,...
        tb_electrodes.y,...
        tb_electrodes.z];
end

disp('All electrodes positions are converted to a matrix.')

%% load mni305 pial
% Freesurfer subjects directory
FSsubjectsdir = fullfile(myDataPath.CCEPpath,'freesurfer');

% load mni305 pial
[Lmnipial_vert,Lmnipial_face] = read_surf(fullfile(FSsubjectsdir,'lh.pial'));
[Rmnipial_vert,Rmnipial_face] = read_surf(fullfile(FSsubjectsdir,'rh.pial'));

% housekeeping
clear FSsubjectsdir

%% load visually checked N1s

clear dataBase

% load visually checked N1s of SPES-clinical
if ~exist(fullfile(myDataPath.CCEPpath,'checkedN1s',...
        [sub_label,'_ses-1_task-SPESclin_N1sChecked_comb.mat']),'file')

    error('You should first run PRIOS01_pipeline_preprocess.m by two observers, PRIOS02_visualCheckN1s.m and PRIOS03_pipeline_agreement.m before you can run this script.')

else

    tmp = load(fullfile(myDataPath.CCEPpath,'checkedN1s',...
        [sub_label,'_ses-1_task-SPESclin_N1sChecked_comb.mat']));

    dataBase.sub_label = sub_label;
    dataBase.cc_stimchans = tmp.cc_stimchans;
    dataBase.cc_stimsets = tmp.cc_stimsets;
    dataBase.spesClin  = tmp;

end

% load visually checked N1s of SPES-propofol
if ~exist(fullfile(myDataPath.CCEPpath,'checkedN1s',...
        [sub_label,'_ses-1_task-SPESprop_N1sChecked_comb.mat']),'file')

    error('You should first run PRIOS01_pipeline_preprocess.m by two observers, PRIOS02_visualCheckN1s.m and PRIOS03_pipeline_agreement.m before you can run this script.')

else

    tmp = load(fullfile(myDataPath.CCEPpath,'checkedN1s',...
        [sub_label,'_ses-1_task-SPESprop_N1sChecked_comb.mat']));

    dataBase.sub_label = sub_label;
    dataBase.cc_stimchans = tmp.cc_stimchans;
    dataBase.cc_stimsets = tmp.cc_stimsets;
    dataBase.spesProp = tmp;
end

% housekeeping
clear all_sublabels tmp

%% number of evoked responses per stimulus pair

% find all non-nans and add this binary matrix in ccepsClin/ccepsProp
ccepsClin = ~isnan(dataBase.spesClin.n1_peak_sample_check);
ccepsProp = ~isnan(dataBase.spesProp.n1_peak_sample_check);

outdegreeClin = sum(ccepsClin,1);
outdegreeProp = sum(ccepsProp,1);

[~,rankOutdegreeClin] = ismember(outdegreeClin,sort(unique(outdegreeClin),'descend'));
[~,rankOutdegreeProp] = ismember(outdegreeProp,sort(unique(outdegreeProp),'descend'));

fprintf('Subject: %s \nhighest outdegreeClin: %2.0f --> outdegreeProp: %2.0f, ranking Prop: %2.0f \nhighest outdegreeProp: %2.0f --> outdegreeClin: %2.0f, ranking Clin: %2.0f \n',...
    dataBase.sub_label, outdegreeClin(find(rankOutdegreeClin == 1,1,'first')), outdegreeProp(find(rankOutdegreeClin == 1,1,'first')), rankOutdegreeProp(find(rankOutdegreeClin == 1,1,'first')),...
    outdegreeProp(find(rankOutdegreeProp == 1,1,'first')), outdegreeClin(find(rankOutdegreeProp == 1,1,'first')), rankOutdegreeClin(find(rankOutdegreeProp == 1,1,'first')))

% housekeeping
clear tbl p nSub chi2 ccepsClin ccepsProp kk pSort pInd pVals thisVal m

%% plot subject electrodes on mni brain

close all

mkrsz = 30;

% get hemisphere for each electrode
hemi = tb_electrodes.hemisphere;

% set the view for the correct hemisphere
if isequal(hemi{1},'L')
    g.faces = Lmnipial_face+1; % correct for zero index
    g.vertices = Lmnipial_vert;
    v_d = ([270 0]);
elseif isequal(hemi{1},'R')
    g.faces = Rmnipial_face+1; % correct for zero index
    g.vertices = Rmnipial_vert;
    v_d = ([90 0]);
end

% make the electrodes pop out of the brain cortex
a_offset = .2*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);

figure
ieeg_RenderGifti(g);

hold on

% all electrodes:
plot3(els(:,1), ...
    els(:,2), ...
    els(:,3), ...
    '.','Color', 'k','MarkerSize',mkrsz)

% electrodes with high outdegree:

stimChansClin = dataBase.cc_stimsets(find(rankOutdegreeClin == 1,1,'first'),:);
stimChansProp = dataBase.cc_stimsets(find(rankOutdegreeProp == 1,1,'first'),:);

plot3(els(stimChansClin,1), els(stimChansClin,2), els(stimChansClin,3), ...
    '.','Color', cmap(1,:),'MarkerSize',0.8*mkrsz)
plot3(els(stimChansProp,1), els(stimChansProp,2), els(stimChansProp,3), ...
    '.','Color', cmap(2,:),'MarkerSize',0.8*mkrsz)
hold off

ieeg_viewLight(v_d(1),v_d(2))

figureName = sprintf('%ssupfig3_outDegreeBrain_%s',...
    myDataPath.Figures,sub_label);

if ~exist(myDataPath.Figures,"dir")

    mkdir(myDataPath.Figures)
end

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)

%% end of script