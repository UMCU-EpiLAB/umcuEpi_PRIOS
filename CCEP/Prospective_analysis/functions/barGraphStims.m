function barGraphStims(dataBase,myDataPath)
%% Make horizontal bar graph for the in-degree, out-degree and betweenness centrality
% Combine all patients toghether

mode = {'In-degree','Out-degree','Betweenness Centrality'};
figure('Position',[719,2,1201,1055])

% Remove patients with too low interobserver agreement
skip_pat = zeros(size(dataBase,2),1);

% Remove patients that have an interobserver agreement lower than 0.6
for i = 1:size(dataBase,2)      
   if dataBase(i).ccep_clin.Ckappa < 0.6 && dataBase(i).ccep_prop.Ckappa <0.6
       skip_pat(i,:) = 1;
   end
end

dataBase(skip_pat==1) = [];

for m = 1:size(mode,2)
    
    start_row = 1;
    dataBarPlot = [];

    for subj = 1:size(dataBase,2)

        if isequal(mode{m},'In-degree')
            data_clin = dataBase(subj).agreement_parameter.indegreeN_Clin;
            data_prop = dataBase(subj).agreement_parameter.indegreeN_Prop;

        elseif isequal(mode{m},'Out-degree')
            data_clin = dataBase(subj).agreement_parameter.outdegreeN_Clin;
            data_prop = dataBase(subj).agreement_parameter.outdegreeN_Prop;
    
        elseif isequal(mode{m},'Betweenness Centrality')
            data_clin = dataBase(subj).agreement_parameter.BCN_Clin;
            data_prop = dataBase(subj).agreement_parameter.BCN_Prop;

        end
    
        % The data of all patients are concatenated in one column so it is necessary
        % to know what the last row of values is for each patient so the
        % values of the next patient can start right after. 
        lastrow = size(data_clin,2) + start_row-1;

        % Clin is second column, prop is first column (required to get
        % clin on the left side of the bar plot)
        dataBarPlot(start_row:lastrow, 2) = data_clin *(-1);   
        dataBarPlot(start_row:lastrow, 1) = data_prop;
        
        start_row = size(dataBarPlot,1)+1;
           
    end
       
    % Sort based on the clinical-SPES results in the second column. 
    dataBarPlot_sorted = sortrows(dataBarPlot,2,'descend');

    subplot(3,1,m)
    barh(dataBarPlot_sorted(:,2),'FaceColor',[17/255,145/255,250/255],'EdgeColor',[17/255,145/255,250/255]);     % Left --> Clinical
    
    hold on
    barh(dataBarPlot_sorted(:,1), 'FaceColor',[252/255,96/255,57/255],'EdgeColor',[252/255,96/255,57/255]);    % Right --> propofol
        
    if m == 2
        % Plot second subplot in same line as the first
        xlim([-0.379722518644821,0.620277481355179])

    end

    % Alter the ticks.
    xticks = get(gca, 'xtick');
  
    % Get the current labels
    labels = get(gca, 'xtickLabel'); 
    
    if ischar(labels)
        labels = cellstr(labels); 
    end
    
    % Figure out which ones we need to change
    toscale = xticks < 0;
    
    % Replace the text for the ones < 0
    labels(toscale) = arrayfun(@(x)sprintf('%1.2f', x), ...
                               abs(xticks(toscale) * -1 ), 'uniformoutput', false);
    
    % Update the tick locations and the labels
    set(gca, 'xtick', xticks, 'xticklabel', labels)       
    legend({'Clinical-SPES','Propofol-SPES'},'FontSize',8,'Location','southeast','Box','off')

    % To determine the significance for the ranking
    % Use Spearman 'rows','pairwise' to ensure that row's with NaN's in both columns are not considered in the analysis.
    [rho, pval] = corr(dataBarPlot(:,1), -dataBarPlot(:,2),'Type','Spearman','rows','pairwise');      
    
    if pval < 0.001
            title(sprintf('%s, p = <0.001, r_s = %1.3f', mode{m}, rho),'FontSize',12)               
    else
            title(sprintf('%s, p = %1.3f, r_s = %1.3f', mode{m}, pval, rho),'FontSize',12)
    end 

    set(gca,'ycolor','none','XGrid','on','Box','off')        
end



% Save figure
outlabel=sprintf('allsubs_barGraph.png');
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/BarGraph/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')

end