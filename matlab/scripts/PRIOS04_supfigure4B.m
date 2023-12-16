%% plotting a figure with the brain and electrodes on the brain
% in this script, the electrodes in mni305 space are plotted on a
% freesurfer brain. This figure is also displayed in the article, in Figure
% 1A. 

% 1. set paths
% 2. select subject
% 3. add electrodes.tsv (see BIDS structure)
% 4. load mni305 pial
% 5. plot subject and electrodes on the freesurfer brain

close all
clear
clc

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

% sub_label = 'sub-PRIOS01'; 
% % stimChan = {'TL13','TL14'};
% % respChan = {'TL15'};
% stimChan = {'TH10','TH11'};
% respChan = {'TH03'};

% sub_label = 'sub-PRIOS02'; 
% stimChan = {'T02','T03'};
% respChan = {'T10'};
% stimChan = {'sT2','sT3'};
% respChan = {'sT1'};
% v_d = ([-90 -90]);

% sub_label = 'sub-PRIOS03'; 
% stimChan = {'F19','F20'};
% respChan = {'F12','F18'};
% stimChan = {'F05','F06'};
% respChan = {'F07'};
% 
% sub_label = 'sub-PRIOS04'; 
% % stimChan = {'IH11','IH12'};
% % respChan = {'IH03','IH04'};
% % v_d = ([90 0]);
% stimChan = {'C14','C15'};
% respChan = {'C22'};
% v_d = ([-90 45]);
% 
sub_label = 'sub-PRIOS05'; 
% stimChan = {'C36','C37'};
% respChan = {'C28'};
stimChan = {'C17','C18'};
respChan = {'C09','C10'};
% 
% sub_label = 'sub-PRIOS09'; 
% % stimChan = {'vT1','vT2'};
% % respChan = {'sT1'};
% % v_d = ([-90 -90]);
% stimChan = {'vT2','vT3'};
% respChan = {'vT1'};
% v_d = ([-90 -90]);


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

%% plot subject electrodes on mni brain

close all

mkrsz = 30;

% get hemisphere for each electrode
hemi = tb_electrodes.hemisphere;

% set the view for the correct hemisphere
if isequal(hemi{1},'L')
    g.faces = Lmnipial_face+1; % correct for zero index
    g.vertices = Lmnipial_vert;

    if ~exist('v_d','var')
        v_d = ([270 0]);
    end
elseif isequal(hemi{1},'R')
    g.faces = Rmnipial_face+1; % correct for zero index
    g.vertices = Rmnipial_vert;
    if ~exist('v_d','var')
        v_d = ([90 0]);
    end
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

% stimulated and response electrode:
stimelec = contains(tb_electrodes.name,stimChan);
respelec = find(contains(tb_electrodes.name,respChan));
plot3(els(stimelec==1,1), els(stimelec==1,2), els(stimelec==1,3), ...
    '.','Color', [184/256, 26/256, 93/256],'MarkerSize',0.8*mkrsz)
if size(respelec,1) == 2
    for nResp = 1:size(respelec,1)
        plot3(els(respelec(nResp),1), els(respelec(nResp),2), els(respelec(nResp),3), ...
            '.','Color', cmap(nResp,:),'MarkerSize',0.8*mkrsz)
    end
elseif size(respelec,1) == 1
    plot3(els(respelec,1), els(respelec,2), els(respelec,3), ...
        '.','Color', [0/256, 161/256, 154/256],'MarkerSize',0.8*mkrsz)
end

hold off

ieeg_viewLight(v_d(1),v_d(2))

if size(respChan,2) == 1
    figureName = sprintf('%ssupfig4B_brain_%s_%s-%s_%s',...
        myDataPath.Figures,sub_label,...
        stimChan{:},respChan{:});
elseif size(respChan,2) == 2
    figureName = sprintf('%ssupfig4B_brain_%s_%s-%s_%s_%s',...
        myDataPath.Figures,sub_label,...
        stimChan{:},respChan{:});
end

if ~exist(myDataPath.Figures,"dir")

    mkdir(myDataPath.Figures)
end

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)


%% end of script
