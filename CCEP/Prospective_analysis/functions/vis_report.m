function vis_report(dataBase, myDataPath)
% Function used to group/sort all scripts only used for visualisation of
% results for the report

%% Remove patients with too low inter observer agreement
dataBase_remove = zeros(1,size(dataBase,2));
for s = 1:size(dataBase,2)
     if dataBase(s).ccep_clin.Ckappa <0.6 || dataBase(s).ccep_prop.Ckappa < 0.6
            % Skip because inter observer agreement is too low
        dataBase_remove(:,s) = 1;
     else
        dataBase_remove(:,s) = 0;
     end
end

loc_remove = find(dataBase_remove == 1);
for i = 1:size(loc_remove,2)    
    names_remove(i,:) = dataBase(loc_remove(i)).ccep_clin.sub_label;
end

% Remove patient from dataBase
dataBase(:,loc_remove) = [];


%% Visualise the number of ERs per SPES session per patient with bar graphs
figure('Position',[407,689,939,373])
ax1 = axes('Position',[0.074,0.11,0.9,0.82]);

% Pre-allocation
ERs_tot = zeros(size(dataBase,2)*2,1);
subjects = cell(1,size(dataBase,2));
combin = zeros(size(dataBase,2),3);

for subj = 1:size(dataBase,2)
   
    % prealloction of the column number   
    clin = 2*subj-1;                      
    prop = 2*subj; 
    
    % Place the total number of ERs per patient in an array.
    ERs_tot(clin,2) = sum(sum(~isnan(dataBase(subj).ccep_clin.n1_peak_sample)));
    ERs_tot(prop,2) = sum(sum(~isnan(dataBase(subj).ccep_prop.n1_peak_sample)));
    
    % Place the number of ERs only in clin of only in prop in the second
    % column
    ERs_tot(clin,1) = size(dataBase(subj).elec_in_clin,1);
    ERs_tot(prop,1) = size(dataBase(subj).elec_in_prop,1);


    %% Determine number of responses in both make that blue
    % Then only in clin in lighter blue
    % Only in prop in darker blue
    % So only one column per subject
    combin(subj,2) = ERs_tot(clin,2)- ERs_tot(clin,1);  % in both
    combin(subj,1)= ERs_tot(clin,1);                    % Only in clin
    combin(subj,3)= ERs_tot(prop,1);                    % Only in prop


    subjects{subj} = dataBase(subj).ccep_clin.sub_label;

end

%% Create bar graph
X = categorical(subjects);
b =   bar(ax1,X,combin,'stacked');         % Order: only in clin, in both, only in prop

b(1).FaceColor(:) = [194/255 228/255 255/255];            % only in clinical-SPES
b(2).FaceColor(:) =  [17/255 145/255 250/255];            % Both
b(3).FaceColor(:) = [0/255 66/255 133/255];               % only in propofol-SPES
 
% Place the Number of ERs next to the column
for i = 1:3%:size(dataBase,2)*2
    xtips1 = b(i).XEndPoints;
    ytips1 = b(i).YEndPoints;
    labels1 = string(b(i).YData);
    labels1(ismember(labels1, '0')) = NaN;
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',11) 
end

legend('Only in Clinical-SPES','Both sessions','Only in Propofol-SPES')
ylabel('Number of ERs');
title('Total number of CCEPs evoked per SPES session')

% Save figure
outlabel='CCEPs_per_session.png';
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')



%% Create Violinplot of the ranking of the number of ERs per stimulation pair. 

ERs_perStimp_violin(dataBase,myDataPath,subjects) 
    
    
end
