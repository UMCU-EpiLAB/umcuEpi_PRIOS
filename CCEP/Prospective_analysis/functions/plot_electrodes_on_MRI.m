function plot_electrodes_on_MRI(myDataPath)

% This script can be used to create an MNI cortex (inflated) with
% electrodes in different colors for different locations for all patients
% used in this study. 
%
% Jaap van der Aar, Giulio Castegnaro, Dora Hermes, Dorien van Blooijs, 2019
%

%% Set paths
clc

% get a list of datasets
% theseSubs = ccep_getSubFilenameInfo(myDataPath);
participants_tsv = read_tsv(fullfile(myDataPath.dataPath,'participants.tsv'));

%% Get standardized electrodes through surface based registration or linear
% convert electrodes from patient's individual MRI to MNI305 space

% Freesurfer subjects directory
FSsubjectsdir = fullfile(myDataPath.dataPath,'derivatives','freesurfer');

elec_coords = struct();

for kk = 1:size(participants_tsv,1) 
    disp(['subj ' int2str(kk) ' of ' int2str(size(participants_tsv,1))])
    
    % subject freesurfer dir
    FSdir = fullfile(myDataPath.dataPath,'derivatives','freesurfer',participants_tsv.participant_id{kk},'ses-1',...
        [participants_tsv.participant_id{kk},'_ses-1','_T1w']);
    
    % get electrodes info
    elec_coords(kk).elecs_tsv = readtable(fullfile(myDataPath.dataPath,participants_tsv.participant_id{kk},'ses-1','ieeg',...
        [participants_tsv.participant_id{kk},'_ses-1_electrodes.tsv']),'FileType','text','Delimiter','\t');
    if iscell(elec_coords(kk).elecs_tsv.x)
        elecmatrix = NaN(size(elec_coords(kk).elecs_tsv,1),3);
        for ll = 1:size(elec_coords(kk).elecs_tsv,1)
            if ~isequal(elec_coords(kk).elecs_tsv.x{ll},'n/a')
                elecmatrix(ll,:) = [str2double(elec_coords(kk).elecs_tsv.x{ll}) str2double(elec_coords(kk).elecs_tsv.y{ll}) str2double(elec_coords(kk).elecs_tsv.z{ll})];
            end
        end
    else
        elecmatrix = [elec_coords(kk).elecs_tsv.x elec_coords(kk).elecs_tsv.y elec_coords(kk).elecs_tsv.z];
    end
    nElec = size(elecmatrix,1);
    
    % get hemisphere for each electrode
    these_json = dir(fullfile(myDataPath.dataPath, participants_tsv.participant_id{kk}, 'ses-1','ieeg',[participants_tsv.participant_id{kk},'_ses-1_task-SPESclin*_ieeg.json']));
    ieeg_json = jsonread(fullfile(these_json(1).folder,these_json(1).name));
    if isequal(ieeg_json.iEEGPlacementScheme,'left') || isequal(ieeg_json.iEEGPlacementScheme,'left;')
        hemi = num2cell(repmat('L',nElec,1));
    elseif isequal(ieeg_json.iEEGPlacementScheme,'right')|| isequal(ieeg_json.iEEGPlacementScheme,'right;')
        hemi = num2cell(repmat('R',nElec,1));
    elseif contains(ieeg_json.iEEGPlacementScheme,{'left','right'}) % check with kk=17
        hemi = cell(nElec,1);
        [hemi{:}] = deal('n/a');
        
        schemesplit = strsplit(ieeg_json.iEEGPlacementScheme,';');
        rightcell = find(contains(schemesplit,'right'));
        leftcell = find(contains(schemesplit,'left'));
        
        if rightcell < leftcell
            leftcells = extractAfter(ieeg_json.iEEGPlacementScheme,'left');
            rightcells = extractBetween(ieeg_json.iEEGPlacementScheme,'right','left');
            rightcells = rightcells{:};
        else
            rightcells = extractAfter(ieeg_json.iEEGPlacementScheme,'right');
            leftcells = extractBetween(ieeg_json.iEEGPlacementScheme,'left','right');
            leftcells = leftcells{:};
        end
        
        leftelec = strsplit(leftcells,';');
        leftelec =  leftelec(~cellfun('isempty',leftelec));
        rightelec = strsplit(rightcells,';');
        rightelec = rightelec(~cellfun('isempty',rightelec));
        
        for elec=1:size(leftelec,2)
           C = strsplit(leftelec{elec},{'[',']'});
           elecInd = find(contains(elec_coords(kk).elecs_tsv.name,C{1}));
           [hemi{elecInd}] = deal('L');
        end
        
        for elec=1:size(rightelec,2)
           C = strsplit(rightelec{elec},{'[',']'});
           elecInd = find(contains(elec_coords(kk).elecs_tsv.name,C{1}));
           [hemi{elecInd}] = deal('R');
        end
    end
    elec_coords(kk).hemi = hemi;
    % convert to MNI using surface
    elec_coords(kk).mni_coords = ccep_mni305ThroughFsSphere(elecmatrix,hemi,FSdir,FSsubjectsdir);
    % convert to MNI using linear transformations
    % elec_coords(kk).mni_coords = ccep_mni305linear(elecmatrix,FSdir);
    
end

save(fullfile(myDataPath.CCEPpath,'elec_coordinatesMNI305.mat'),'elec_coords')

%% Start here to make figures
% load MNI electrode positions (saved in previous section), MNI sphere, pial,
% and surface labels

% linear is not so nice...
% load(fullfile(myDataPath.CCEPpath,'derivatives','elec_coordinatesMNI305lin.mat'),'elec_coords')
% surface based is nice:
load(fullfile(myDataPath.CCEPpath,'elec_coordinatesMNI305.mat'),'elec_coords')

% Freesurfer subjects directory
FSsubjectsdir = fullfile(myDataPath.dataPath,'derivatives','freesurfer');

% load mni305 pial
[Lmnipial_vert,Lmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','lh.pial'));
[Rmnipial_vert,Rmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','rh.pial'));

% load mni305 inflated
[Lmniinfl_vert,Lmniinfl_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','lh.inflated'));
[Rmniinfl_vert,Rmniinfl_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','rh.inflated'));

% surface labels
[Lvertices, Llabel, Lcolortable] = read_annotation(fullfile(FSsubjectsdir,'fsaverage','label','lh.aparc.a2009s.annot'));
Lvert_label = Llabel; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
% mapping labels to colortable
for kk = 1:size(Lcolortable.table,1) % 76 are labels
    Lvert_label(Llabel==Lcolortable.table(kk,5)) = kk;
end
[Rvertices, Rlabel, Rcolortable] = read_annotation(fullfile(FSsubjectsdir,'fsaverage','label','rh.aparc.a2009s.annot'));
Rvert_label = Rlabel; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
% mapping labels to colortable
for kk = 1:size(Rcolortable.table,1) % 76 are labels
    Rvert_label(Rlabel==Rcolortable.table(kk,5)) = kk;
end

%% add all electrodes labels and left or right hemisphere into one 
% variable: allmni_coords and allmni_coords_infl

allmni_coords = [];
allmni_coords_infl = [];

allmni_labels = [];
all_hemi = [];
for kk = 1:length(elec_coords)
    Destrieux_label = elec_coords(kk).elecs_tsv.Destrieux_label;
    if iscell(Destrieux_label)
        for ll = 1:size(Destrieux_label,1)
            if ischar(Destrieux_label{ll})
                if isequal(Destrieux_label,'n/a') 
                    Destrieux_label{ll} = NaN;
                else
                    Destrieux_label{ll} = str2double(Destrieux_label{ll});
                end
            end
        end
        Destrieux_label = cell2mat(Destrieux_label);
    end
    
    allmni_coords = [allmni_coords; elec_coords(kk).mni_coords]; %#ok<AGROW>
    allmni_labels = [allmni_labels; Destrieux_label]; %#ok<AGROW>
    all_hemi = [all_hemi; elec_coords(kk).hemi]; %#ok<AGROW>
    
    % run through all coordinates and find the inflated points
    temp_inflated = NaN(size(elec_coords(kk).mni_coords));
    for ll = 1:size(Destrieux_label,1)
        if isequal(elec_coords(kk).hemi{ll},'L')
            [~,min_ind] = min(sqrt(sum((Lmnipial_vert-elec_coords(kk).mni_coords(ll,:)).^2,2)));
            temp_inflated(ll,:) = Lmniinfl_vert(min_ind,:);
        elseif isequal(elec_coords(kk).hemi{ll},'R')
            [~,min_ind] = min(sqrt(sum((Rmnipial_vert-elec_coords(kk).mni_coords(ll,:)).^2,2)));
            temp_inflated(ll,:) = Rmniinfl_vert(min_ind,:);
        end
    end
    
    allmni_coords_infl = [allmni_coords_infl; temp_inflated]; %#ok<AGROW>
end

%% labels for electrode areas we want to color

% categorize anatomical regions
ccep_categorizeAnatomicalRegions

lroi_label = Lvert_label;
lroi_label(ismember(lroi_label,roi_temporal+1)) = 200;
lroi_label(ismember(lroi_label,roi_central+1)) = 300;
lroi_label(lroi_label<100) = 0;
lroi_label = lroi_label/100;

%% Plot figure with left pial with electrodes in mni space

v_d = [270 0];

figure
gl.faces = Lmnipial_face+1;
gl.vertices = Lmnipial_vert;
gl = gifti(gl);
tH = ieeg_RenderGifti(gl); %#ok<NASGU>

% make sure electrodes pop out
a_offset = .1*max(abs(allmni_coords(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords+repmat(a_offset,size(allmni_coords,1),1);      
% els = allmni_coords;

ieeg_elAdd(els(ismember(all_hemi,'L') & ~ismember(allmni_labels,[roi_temporal roi_frontal roi_central roi_parietal]),:),'k',10)
% set(tH,'FaceAlpha',.5) % make transparent
ieeg_elAdd(els(ismember(all_hemi,'L') & ismember(allmni_labels,roi_temporal),:),[0 0 .8],15)
ieeg_elAdd(els(ismember(all_hemi,'L') & ismember(allmni_labels,roi_frontal),:),[1 .8 0],15)
ieeg_elAdd(els(ismember(all_hemi,'L') & ismember(allmni_labels,roi_central),:),[.8 .3 0],15)
ieeg_elAdd(els(ismember(all_hemi,'L') & ismember(allmni_labels,roi_parietal),:),[0 .5 0],15)
ieeg_viewLight(v_d(1),v_d(2))

figureName = fullfile(myDataPath.CCEPpath,'render','leftMNIpial');

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)


%% Plot figure with right pial with electrodes in mni space
v_d = [96 6];

figure
gr.faces = Rmnipial_face+1;
gr.vertices = Rmnipial_vert;
gr = gifti(gr);
tH = ieeg_RenderGifti(gr); %#ok<NASGU>

% make sure electrodes pop out
a_offset = .5*max(abs(allmni_coords(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords+repmat(a_offset,size(allmni_coords,1),1);      

ieeg_elAdd(els(ismember(all_hemi,'R') & ~ismember(allmni_labels,[roi_temporal roi_frontal roi_central roi_parietal]),:),'k',10)
% set(tH,'FaceAlpha',.5) % make transparent
ieeg_elAdd(els(ismember(all_hemi,'R') & ismember(allmni_labels,roi_temporal),:),[0 0 .8],15)
ieeg_elAdd(els(ismember(all_hemi,'R') & ismember(allmni_labels,roi_frontal),:),[1 .8 0],15)
ieeg_elAdd(els(ismember(all_hemi,'R') & ismember(allmni_labels,roi_central),:),[.8 .3 0],15)
ieeg_elAdd(els(ismember(all_hemi,'R') & ismember(allmni_labels,roi_parietal),:),[0 .5 0],15)
ieeg_viewLight(v_d(1),v_d(2))

figureName = fullfile(myDataPath.CCEPpath,'render','rightMNIpial'); 

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)

end