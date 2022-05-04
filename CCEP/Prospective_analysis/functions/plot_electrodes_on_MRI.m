function plot_electrodes_on_MRI(myDataPath, table_latency, dataBase, av_lat_elec)

% % This script can be used to create an MNI cortex (inflated) with
% % electrodes in different colors for different locations for all patients
% % used in this study. 
% %
% % Jaap van der Aar, Giulio Castegnaro, Dora Hermes, Dorien van Blooijs, 2019
% %
% 
% %% Set paths
% % clc
% 
% % get a list of datasets
% % theseSubs = ccep_getSubFilenameInfo(myDataPath);
participants_tsv = read_tsv(fullfile(myDataPath.dataPath,'participants.tsv'));
% participants_tsv(8:end,:) = [];
% 
% %% Get standardized electrodes through surface based registration or linear
% % convert electrodes from patient's individual MRI to MNI305 space
% 
% % Freesurfer subjects directory
% FSsubjectsdir = fullfile(myDataPath.dataPath,'derivatives','freesurfer');
% 
% elec_coords = struct();
% 
% for kk = 1:size(participants_tsv,1)
%     disp(['subj ' int2str(kk) ' of ' int2str(size(participants_tsv,1))])
%     
%     % subject freesurfer dir
%     FSdir = fullfile(myDataPath.dataPath,'derivatives','freesurfer',participants_tsv.participant_id{kk},'ses-1',...
%         [participants_tsv.participant_id{kk},'_ses-1','_T1w']);
%     
%     % get electrodes info
%     elec_coords(kk).elecs_tsv = readtable(fullfile(myDataPath.dataPath,participants_tsv.participant_id{kk},'ses-1','ieeg',...
%         [participants_tsv.participant_id{kk},'_ses-1_electrodes.tsv']),'FileType','text','Delimiter','\t');
%     if iscell(elec_coords(kk).elecs_tsv.x)
%         elecmatrix = NaN(size(elec_coords(kk).elecs_tsv,1),3);
%         for ll = 1:size(elec_coords(kk).elecs_tsv,1)
%             if ~isequal(elec_coords(kk).elecs_tsv.x{ll},'n/a')
%                 elecmatrix(ll,:) = [str2double(elec_coords(kk).elecs_tsv.x{ll}) str2double(elec_coords(kk).elecs_tsv.y{ll}) str2double(elec_coords(kk).elecs_tsv.z{ll})];
%             end
%         end
%     else
%         elecmatrix = [elec_coords(kk).elecs_tsv.x elec_coords(kk).elecs_tsv.y elec_coords(kk).elecs_tsv.z];
%     end
%     nElec = size(elecmatrix,1);
%     
%     % get hemisphere for each electrode
%     these_json = dir(fullfile(myDataPath.dataPath, participants_tsv.participant_id{kk}, 'ses-1','ieeg',[participants_tsv.participant_id{kk},'_ses-1_task-SPESclin*_ieeg.json']));
%     ieeg_json = jsonread(fullfile(these_json(1).folder,these_json(1).name));
%     if isequal(ieeg_json.iEEGPlacementScheme,'left') || isequal(ieeg_json.iEEGPlacementScheme,'left;')
%         hemi = num2cell(repmat('L',nElec,1));
%     elseif isequal(ieeg_json.iEEGPlacementScheme,'right')|| isequal(ieeg_json.iEEGPlacementScheme,'right;')
%         hemi = num2cell(repmat('R',nElec,1));
%     elseif contains(ieeg_json.iEEGPlacementScheme,{'left','right'}) % check with kk=17
%         hemi = cell(nElec,1);
%         [hemi{:}] = deal('n/a');
%         
%         schemesplit = strsplit(ieeg_json.iEEGPlacementScheme,';');
%         rightcell = find(contains(schemesplit,'right'));
%         leftcell = find(contains(schemesplit,'left'));
%         
%         if rightcell < leftcell
%             leftcells = extractAfter(ieeg_json.iEEGPlacementScheme,'left');
%             rightcells = extractBetween(ieeg_json.iEEGPlacementScheme,'right','left');
%             rightcells = rightcells{:};
%         else
%             rightcells = extractAfter(ieeg_json.iEEGPlacementScheme,'right');
%             leftcells = extractBetween(ieeg_json.iEEGPlacementScheme,'left','right');
%             leftcells = leftcells{:};
%         end
%         
%         leftelec = strsplit(leftcells,';');
%         leftelec =  leftelec(~cellfun('isempty',leftelec));
%         rightelec = strsplit(rightcells,';');
%         rightelec = rightelec(~cellfun('isempty',rightelec));
%         
%         for elec=1:size(leftelec,2)
%            C = strsplit(leftelec{elec},{'[',']'});
%            elecInd = find(contains(elec_coords(kk).elecs_tsv.name,C{1}));
%            [hemi{elecInd}] = deal('L');
%         end
%         
%         for elec=1:size(rightelec,2)
%            C = strsplit(rightelec{elec},{'[',']'});
%            elecInd = find(contains(elec_coords(kk).elecs_tsv.name,C{1}));
%            [hemi{elecInd}] = deal('R');
%         end
%     end
%     elec_coords(kk).hemi = hemi;
%     % convert to MNI using surface
%     elec_coords(kk).mni_coords = ccep_mni305ThroughFsSphere(elecmatrix,hemi,FSdir,FSsubjectsdir);
%     % convert to MNI using linear transformations
%     % elec_coords(kk).mni_coords = ccep_mni305linear(elecmatrix,FSdir);
%     
% end
% 
% save(fullfile(myDataPath.CCEPpath,'elec_coordinatesMNI305.mat'),'elec_coords')

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
Destrieux_label_pat = NaN(150,size(elec_coords,2)); % 150 is an extimation

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

    % Save all destrieux labes per patient
    % Required to determine median per lobe and min and max per patient
    Destrieux_label_pat(1:size(Destrieux_label,1),kk) = Destrieux_label;
    
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

%% plot electrodes on brain that have response in prop and not in clin

% prop_resp_MRI(dataBase, myDataPath)


%% Plot figure with left pial with electrodes in mni space

number_N1_MRI(Lmnipial_vert, Rmnipial_vert, Lmnipial_face, Rmnipial_face, allmni_coords, all_hemi, dataBase, myDataPath)


%% Electrodes on brain with N1-latency per electrode

latency_N1_MRI(av_lat_elec, Lmnipial_vert, Rmnipial_vert, Lmnipial_face, Rmnipial_face, allmni_coords, all_hemi, dataBase, myDataPath)


%% Determine the number of CCEPs evoked BY and ON electrodes located on a lobe
close all;
CCEPS_per_lobe_stimp(myDataPath, dataBase, Destrieux_label_pat, roi_central, roi_frontal, roi_occipital, roi_temporal, roi_parietal)
CCEPS_per_lobe_elec(myDataPath, dataBase, Destrieux_label_pat, roi_central, roi_frontal, roi_occipital, roi_temporal, roi_parietal)

%% Number of ERs per lobe
% Preallocation (NaN's are later removed)
% ERs_per_lobe_clin = NaN(100,size(dataBase,2));
% ERs_per_lobe_prop = NaN(100,size(dataBase,2));
% 
% % Create matrix with the nmber of ERs per electrode of all patients
% for pat = 1:size(dataBase,2)
%     ERs_per_lobe_clin(1:size(dataBase(pat).agreement_parameter.ERs_elecClin,2),pat) = dataBase(pat).agreement_parameter.ERs_elecClin' ; % number of ERs per electrode
%     ERs_per_lobe_prop(1:size(dataBase(pat).agreement_parameter.ERs_elecProp,2),pat) = dataBase(pat).agreement_parameter.ERs_elecProp' ; % number of ERs per electrode
% 
% end
% 
% % Concatenate to one array
% ERs_per_lobe_clin = ERs_per_lobe_clin(:);
% ERs_per_lobe_prop = ERs_per_lobe_prop(:);
% 
% % Remove NaN's
% ERs_per_lobe_clin(isnan(ERs_per_lobe_clin)) = [];
% ERs_per_lobe_prop(isnan(ERs_per_lobe_prop)) = [];
% 
% 
% %% Average N1 latency per lobe
% % Determine for clincal-SPES and propofol-SPES
% 
% allmni_labels = allmni_labels(~isnan(allmni_labels)); %destrieux labels per electrode
% 
% % make one array with averaged latencies of all electrodes of all patients
% % For clinical SPES
% pat = 1;
% for c = 1:2:size(av_lat_elec,2)
%     N1_latency_clin{c} = av_lat_elec(1:size(dataBase(pat).ccep_clin.ch,1), c);
%     pat = pat+1;
% end
% 
% n1_lat_clin = vertcat(N1_latency_clin{:});
% 
% % For propofol SPES
% pat = 1;
% for c = 2:2:size(av_lat_elec,2)
%     N1_latency_prop{c} = av_lat_elec(1:size(dataBase(pat).ccep_prop.ch,1), c);
%     pat = pat+1;
% end
% 
% n1_lat_prop = vertcat(N1_latency_prop{:});
% 
% mode = {'Temporal','Frontal','Parietal','Central'};
% N1_lobe_clin = nan(100,size(mode,2)); % 100 is an estimation
% N1_lobe_prop = nan(100,size(mode,2)); % 100 is an estimation
% 
% for m = 1:size(mode,2)
% 
%     if isequal(mode{m},'Temporal')
%         region = roi_temporal;
% 
%     elseif isequal(mode{m},'Frontal')
%          region = roi_frontal;
%          
%     elseif isequal(mode{m},'Parietal')
%          region = roi_parietal;
%            
%     elseif isequal(mode{m},'Central')
%          region = roi_central;
% 
%     end
% 
%     idx_lobe = ismember(allmni_labels, region);
%     
%     % N1 latency per lobe
%     N1_lobe_clin(1:sum(idx_lobe),m) = n1_lat_clin(idx_lobe);
%     N1_lobe_prop(1:sum(idx_lobe),m) = n1_lat_prop(idx_lobe);
% 
%     % Number of ERs per lobe
%     sum_ERs_per_lobe_clin(m) = sum(ERs_per_lobe_clin(idx_lobe));
%     median(ERs_per_lobe_clin(idx_lobe))
% 
%     min_ERs_clin = min(ERs_per_lobe_clin(idx_lobe));
%     max_ERs_clin = max(ERs_per_lobe_clin(idx_lobe));
%     
%     sum_ERs_per_lobe_prop(m) = sum(ERs_per_lobe_prop(idx_lobe));
%     median(ERs_per_lobe_prop(idx_lobe))
%     min_ERs_prop = min(ERs_per_lobe_prop(idx_lobe));
%     max_ERs_prop = max(ERs_per_lobe_prop(idx_lobe));
% 
%     fprintf('Clinical-SPES: Median latency in %s lobe = %1.1f ms, contained %1.0f electrodes \n', mode{m}, median(n1_lat_clin(idx_lobe),'omitnan'), sum(idx_lobe));
%     fprintf('Propofol-SPES: Median latency in %s lobe = %1.1f ms, \n', mode{m}, median(n1_lat_prop(idx_lobe),'omitnan'));
%     
%     fprintf('Clinical-SPES: Total ERs in %s lobe = %1.0f, min = %1.0f, max = %1.0f, \n', mode{m}, sum_ERs_per_lobe_clin(m), min_ERs_clin, max_ERs_clin);
%     fprintf('Propofol-SPES: Total ERs in %s lobe = %1.0f, min = %1.0f, max = %1.0f, \n', mode{m}, sum_ERs_per_lobe_prop(m), min_ERs_prop, max_ERs_prop);
% 
% 
%     % Significance
%     p(m) = signrank(n1_lat_clin(idx_lobe), n1_lat_prop(idx_lobe)) ; 
%     if p(m) < 0.05 
%         fprintf('The p-value between clin and prop for %s lobe = %1.4f. This means a significant difference \n', mode{m},p(m));
%     else
%         fprintf('The p-value between clin and prop for %s lobe = %1.4f. This means NO significant difference \n',mode{m}, p(m));
%     end
% 
% 
%     % Number of electrodes brain lobe
%     % Determine median number of electrodes per lobe and min and max per patient
%     total_lobe = sum(sum(ismember(Destrieux_label_pat, region)));
%     med_lobe = mean(sum(ismember(Destrieux_label_pat, region)));
%     min_lobe = min(sum(ismember(Destrieux_label_pat, region)));
%     max_lobe = max(sum(ismember(Destrieux_label_pat, region)));
% 
%     fprintf('%s lobe Total electrodes = %1.0f, mean = %1.1f, min/max = %1.0f/%1.0f \n', mode{m}, total_lobe, med_lobe, min_lobe, max_lobe);
% 
%     fprintf('****************** next lobe******************\n')


end


% %% Display the latency per lobe in a violin plot
% figure('Position',[205,424,1530,638]);
% N1_lobe_combined = zeros(100,size(mode,2)*2);
% N1_lobe_combined(:,1:2:size(N1_lobe_combined,2)) = N1_lobe_clin;
% N1_lobe_combined(:,2:2:size(N1_lobe_combined,2)) = N1_lobe_prop;
% 
% 
% violins = violinplot(N1_lobe_combined) ;
% for i = 1:2:size(mode,2)*2
% %     violins(1,i*2).ViolinColor = violins(1,i*2-1).ViolinColor ;
%     violins(1,i).ViolinColor = [1 0 0];
%     violins(1,i+1).ViolinColor = [0 0 1];
% 
% end
% 
% count = 1;
% ymax = max(max(N1_lobe_combined));
%  
% for m=1:size(mode,2)
%         if p(m)   < 0.001 
%             text(count+0.5,ymax-0.2,'**','FontSize',20,'FontWeight','bold')
%             plot(count+0.1:0.1:count+0.9, ymax-0.43*ones(9,1),'k','LineWidth',2)
%         
%         elseif p(m)  < 0.05 
%             text(count+0.5,ymax-0.2,'*','FontSize',20,'FontWeight','bold')
%             plot(count+0.1:0.1:count+0.9, ymax-0.43*ones(9,1),'k','LineWidth',2)
%                 
%         end
%         count = count+2;
% end
% 
% ax = gca;
% ax.XAxis.FontSize = 12;
% ax.YAxis.FontSize = 12;
% ax.XAxis.FontWeight = 'bold';
% ax.YAxis.FontWeight = 'bold';
% 
% % Set double xlabel
% ax.XTick = 1.5:2:size(N1_lobe_combined,2);
% ax.XTickLabel  = mode'; 
% 
% % Display medians on second row beneath the figure
% medians = median(N1_lobe_combined,'omitnan');
% 
% ymin = min(ylim);
% y_range = diff(ylim);
% x_as = 1:size(N1_lobe_combined,2);
% size_pat = size(N1_lobe_combined,2); 
% second_row_txt = cellstr(strsplit(num2str(medians,'%.1f '),' '));
% text([(x_as(1)-x_as(2))*0.5 x_as], ones(1,size_pat+1)*ymin-0.08*y_range, ['Median' second_row_txt],'HorizontalAlignment','center','FontSize', 12)
% 
% sum_ERs_text = [sum_ERs_per_lobe_clin(:) sum_ERs_per_lobe_prop(:)]';
% sum_ERs_text = sum_ERs_text(:)';
% third_row_txt = cellstr(strsplit(num2str(sum_ERs_text,'%1.0f '),' '));
% text([(x_as(1)-x_as(2))*0.5 x_as], ones(1,size_pat+1)*ymin-0.12*y_range, ['Total CCEPs' third_row_txt],'HorizontalAlignment','center','FontSize', 12)
% 
% 
% % Draw lines between patients 
% ymax = max(ylim); 
% for i = 1:2:size(x_as,2)
%     x1 = x_as(i)-0.5; 
%     if x1 > 0.5
%         hold on
%         line([x1,x1],[ymin,ymax],'color',[0.8 0.8 0.8]);
%     end
% end
% 
% title(sprintf('N1-Latency evoked on electrodes per lobe '),'FontSize', 15, 'FontWeight', 'bold')
% ylabel('Latency (milliseconds)','FontSize', 15, 'FontWeight', 'bold')
% legend([violins(1).ViolinPlot,violins(2).ViolinPlot], 'Clinical SPES','Propofol SPES','FontSize', 12, 'FontWeight', 'bold','Position',[0.78,0.84,0.12,0.07])
% 
% % Save figure
% outlabel=sprintf('Latency_violin_perLobe.png');
% path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/N1_compare/');
% if ~exist(path, 'dir')
%     mkdir(path);
% end
% saveas(gcf,[path,outlabel],'png')



end