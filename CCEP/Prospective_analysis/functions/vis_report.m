function vis_report(dataBase, myDataPath)
% Function used to group/sort all scripts only used for visualisation of
% results for the report

%% Visualise the number of ERs per SPES session per patient with bar graphs
figure('Position',[407,689,939,373])
ax1 = axes('Position',[0.074,0.11,0.9,0.82]);

% Pre-allocation
ERs_tot = zeros(size(dataBase,2),(size(dataBase,2)*2));

for subj = 1:size(dataBase,2)
   
    % prealloction of the column number   
    clin_colm = 2*subj-1;                      
    prop_colm = 2*subj; 
    
    % Place the total number of ERs per patient in an array.
    ERs_tot(subj,clin_colm) = sum(sum(~isnan(dataBase(subj).ccep_clin.n1_peak_amplitude_check)));
    ERs_tot(subj,prop_colm) = sum(sum(~isnan(dataBase(subj).ccep_prop.n1_peak_amplitude_check)));
       
end

x = {'PRIOS01','PRIOS02','PRIOS03','PRIOS04','PRIOS05','PRIOS06'};
 
% Create bar graph
X = categorical(x);
b =   bar(ax1,X,ERs_tot,1);         % bar(ax1,X,ERs_tot,1);

for i = 1:2:size(ERs_tot,2)
        b(i).FaceColor(:) =  [0.5843 0.8157 0.9882];          %[0 0 1]
        b(i+1).FaceColor(:) = [0.9882 0.6157 0.5843];                        %[1 0 0]
end
 

% Place the Number of ERs next to the column
for i = 1:2:12
    xtips1 = b(i).XEndPoints;
    ytips1 = b(i).YEndPoints;
    labels1 = string(b(i).YData);
    labels1(ismember(labels1, '0')) = NaN;
    text(xtips1,ytips1,labels1,'HorizontalAlignment','right',...
    'VerticalAlignment','bottom','FontSize',11,'fontweight','Bold') 
end

for i = 2:2:12
    xtips1 = b(i).XEndPoints;
    ytips1 = b(i).YEndPoints;
    labels1 = string(b(i).YData);
    labels1(ismember(labels1, '0')) = NaN;
    text(xtips1,ytips1,labels1,'HorizontalAlignment','left',...
    'VerticalAlignment','bottom','FontSize',11,'fontweight','Bold') 
end
legend('SPES-clin','SPES-Prop')
ylabel('Number of ERs');
title('Total number of ERs evoked per SPES session')


% Save figure
outlabel='ERs_per_stimp.jpg';
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'jpg')

end
