function vis_report(dataBase, myDataPath)
% Function used to group/sort all scripts only used for visualisation of
% results for the report

%% Visualise the number of ERs per SPES session per patient with bar graphs
% TODO: total number of ERs does not make much sense, because it depends on
% the number of stimulus pairs and number of channels. Normalizing it to
% ERs per stimulus pairs would make more sense. 
figure('Position',[407,689,939,373])
ax1 = axes('Position',[0.074,0.11,0.9,0.82]);

% Pre-allocation
ERs_tot = zeros(size(dataBase,2)*2,1);
subjects = cell(1,size(dataBase,2));

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

    
    subjects_1{subj} = dataBase(subj).ccep_clin.sub_label;
    subjects{clin} = [dataBase(subj).ccep_clin.sub_label,'_c'];
    subjects{prop} = [dataBase(subj).ccep_clin.sub_label,'_p'];

end
 
plot_tot(1:2:14,1:2) = ERs_tot(1:2:end,:);
plot_tot(1:2:14,3:4) = 0;

plot_tot(2:2:14,3:4 )= ERs_tot(2:2:end,:);
plot_tot(2:2:14,1:2) = 0;

% Create bar graph
X = categorical(subjects);
b =   bar(ax1,X,plot_tot,'stacked');         % bar(ax1,X,ERs_tot,1);

b(1).FaceColor(:) = [0.5 0 0.5];          %[0 0 1]
b(2).FaceColor(:) =  [0.9882 0.6157 0.5843];                        %[1 0 0]
b(3).FaceColor(:) = [0.5 0 0.5];                        %[1 0 0]
b(4).FaceColor(:) = [0.5843 0.8157 0.9882];                        %[1 0 0]
ax1.XTickLabel = [];
 
% Place the Number of ERs next to the column
for i = 1:4%:size(dataBase,2)*2
    xtips1 = b(i).XEndPoints;
    ytips1 = b(i).YEndPoints;
    labels1 = string(b(i).YData);
    labels1(ismember(labels1, '0')) = NaN;
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',11) 
end

legend('Only in one session','SPES-clin','','SPES-prop')
ylabel('Number of ERs');
title('Total number of ERs evoked per SPES session')

ymin = min(ylim);
y_range = diff(ylim)-450;
x_as = 1:2:size(subjects,2);
size_pat = size(subjects_1,2); 
second_row_txt = subjects_1; %cellstr(strsplit(num2str(medians,'%.2f '),' '));
text([-0.5 ;x_as(:)+0.5], ones(1,size_pat+1)*ymin-0.1*y_range, cellstr([' ',second_row_txt]),'HorizontalAlignment','center','FontSize', 10)



% Save figure
outlabel='ERs_per_stimp.png';
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')

%% Visualise the network characteristics in a heatmap to later plot on the MRI
% Create a heatmap of the network characteristics with the outlay of the
% electrodes from the matlabSjabloon in Excel.

% for subj = 1:size(dataBase,2)    
%     heat_map_grid(myDataPath, dataBase(subj).ccep_clin, dataBase(subj).agreement_parameter)
% end

%% Create Violinplot of the ranking of the number of ERs per stimulation pair. 

ERs_perStimp_violin(dataBase,myDataPath) 
    
    
end
