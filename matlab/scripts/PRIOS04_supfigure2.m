%% plotting a figure with the brain and electrodes on the brain
% in this script, the electrodes in mni305 space are plotted on a
% freesurfer brain. This figure is also displayed in the article, in Figure
% 1A. 

% 1. set paths
% 2. select subject
% 3. add electrodes.tsv (see BIDS structure)
% 4. find electrodes that are located on the endpoints of the arcuate
% fasciculus (electrodes with Destrieux labels: 12, 14, 15, 52, 53
% (Frontal), or 34, 36, 37, 38, 72, 73 (Temporal))
% 5. load mni305 pial
% 6. plot subject and electrodes on the freesurfer brain

close all
clear
clc

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

sub_label = 'sub-PRIOS01'; 

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

%% find electrodes in the Arcuate Fasciculus (AF)

if iscell(tb_electrodes.Destrieux_label)
    Destrieux_label = tb_electrodes.Destrieux_label;

    Destrieux_label = str2double(Destrieux_label);

else
    Destrieux_label = tb_electrodes.Destrieux_label;
end

AFChans_F = ismember(Destrieux_label, [12 14 15 52 53]);
AFChans_T = ismember(Destrieux_label, [34 36 37 38 72 73]);

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

% electrodes on AF:
    plot3(els(AFChans_T==1,1), els(AFChans_T==1,2), els(AFChans_T==1,3), ...
        '.','Color', [184/256, 26/256, 93/256],'MarkerSize',0.8*mkrsz)
    plot3(els(AFChans_F==1,1), els(AFChans_F==1,2), els(AFChans_F==1,3), ...
        '.','Color', [0/256, 156/256, 180/256],'MarkerSize',0.8*mkrsz)
hold off

ieeg_viewLight(v_d(1),v_d(2))

figureName = sprintf('%ssupfig2_brain_%s',...
    myDataPath.Figures,sub_label);

if ~exist(myDataPath.Figures,"dir")

    mkdir(myDataPath.Figures)
end

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)

%% end of script