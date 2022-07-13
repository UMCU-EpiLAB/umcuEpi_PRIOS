function [statistics, rank] = statistical_agreement(myDataPath, agreement_parameter, ccep)
%% Wilcoxon signed rank non parametric test for two paired groups
% For the number of ERs detected in the SPESclin and SPESprop
% null hypothesis is that the two means are the same
 
% Use the Wilcoxon Signed rank test for all patients when the larger part
% of the population is not normally distributed.
% Wilcoxon tests the null hypothesis that data in x and y are samples from continuous distributions with equal medians
% So significant value indicates a significant difference between 2 parameters
p = signrank(agreement_parameter.ERs_stimpClin, agreement_parameter.ERs_stimpProp) ;           

if p<0.05
    fprintf('Wilcoxon signed rank test for the number of ERs evoked per stim-pair gives p-value = %1.4f. There is a significant difference for %s \n', p, ccep.sub_label);
else
    fprintf('Wilcoxon signed rank test for the number of ERs evoked per stim-pair gives p-value = %1.4f. There is NO significant difference for %s \n', p, ccep.sub_label);
end
     
%% Spearman test for ranking of stimpairs
% Rank the stimuli from stimulus pairs with most evoked ERs to stimulus
% pairs with least evoked ERs, NaN's are placed at the end of the ranking.
mode = {'SPES_clin','SPES_prop'};
rank = struct;

for i = 1:size(mode,2)
    
     if strcmp(mode{i},'SPES_clin')
        ERs = agreement_parameter.ERs_stimpClin;
     elseif strcmp(mode{i},'SPES_prop')
        ERs = agreement_parameter.ERs_stimpProp;
     end

    rank.(mode{i})(:,1:2) = ccep.stimsets_avg;
    rank.(mode{i})(:,3) = ERs;
    % Rank the number of ERs per stimulation pair from the moest to the lease with NaN's at the end of the ranking
    [~, order] = sort(rank.(mode{i})(:,3), 'descend','MissingPlacement','last');          
    rank.(['sort_' mode{i}]) = rank.(mode{i})(order, :);
        
    % If the next stimpair has the same number of ERs, give it the same order number
    rank.(['sort_' mode{i}])(1,4) = 1;
    for j = 2:size(rank.(['sort_' mode{i}]),1)
        if rank.(['sort_' mode{i}])(j,3) == rank.(['sort_' mode{i}])(j-1,3)          
            rank.(['sort_' mode{i}])(j,4) = rank.(['sort_' mode{i}])(j-1,4);
        else
            rank.(['sort_' mode{i}])(j,4) = j;
        end
    end
    
    rank.(['sort_names_' mode{i}]) = ccep.stimpnames_avg(order);
end

%% Create with lines drawn between the positions of the stimpairs in the two rankings.

% make labeling for figure
for i=1:size(mode,2)
    [~,~,ic] = unique(rank.(['sort_' mode{i}])(:,4));
    groups = splitapply(@(x){x},rank.(['sort_names_' mode{i}])',ic);
    n = cell(size(groups));
    [n{:}] = deal(' ');
    newgroups = reshape(horzcat(groups,n)',size(groups,1)*2,1);
    rank.(['fig_sort_names_' mode{i}]) = vertcat({' '}, newgroups{:});
end
 
% Add number of CCEPS behind each name
for i=1:size(mode,2)
    row = 1;
    for r = 1:size(rank.(['fig_sort_names_' mode{i}]),1)
        if isequal(rank.(['fig_sort_names_' mode{i}]){r,:}, ' ')
            rank.(['number_of_CCEPs_' mode{i}])(r,:) =  NaN;        
        else
            rank.(['number_of_CCEPs_' mode{i}])(r,:) =   rank.(['sort_' mode{i}])(row,3) ;
            row = row+1;
        end
    end
end

% COncatenate stimpair name and number of CCEPs
for i=1:size(mode,2)
    for c = 1: size(rank.(['fig_sort_names_' mode{i}]),1)
        if isnan(rank.(['number_of_CCEPs_' mode{i}])(c,1))
            rank.(['name_cceps_' mode{i}]){c,:}  = ' ';
            
        else
           rank.(['name_cceps_' mode{i}]){c,:} = [rank.(['fig_sort_names_' mode{i}]){c,1} ' (' num2str(rank.(['number_of_CCEPs_' mode{i}])(c,1)) ')'];
     
        end
    end
end

% Sorted matrix based on stimpair number, so ranking isunsorted 
for i=1:size(mode,2)
    [~,order] = sortrows(rank.(['sort_' mode{i}])(:,1:2));
    rank.(['unsort_' mode{i}]) = rank.(['sort_' mode{i}])(order,:);         
end

% If the order of the stimulation pairs is not equal, than the ranking
% cannot be compared
if ~isequal(rank.unsort_SPES_clin(:,1:2),rank.unsort_SPES_prop(:,1:2))
    error('Sorting stimulus pairs is incorrect and led to unequal matrices in SPES-clin and SPES-prop')
end


% Test the hypothesis of NO correlation
% When p <0.05, an rho is close to (-)1, approval of the hypothesis that
% there is (negative) correlation between the two columns 
[RHO_stmp,PVAL_stmp] = corr(rank.unsort_SPES_clin(:,3) , rank.unsort_SPES_prop(:,3) ,'Type','Spearman');            % Test the hypothesis that there is a correlation
fprintf('Spearman Corr between stimpair ranking of SPES-clin and SPES-prop gives, p-value = %1.4f, rho = %1.3f, for %s \n', PVAL_stmp, RHO_stmp, ccep.sub_label);

figure('Position',[1074,4,519,1052]);
cm = colormap(parula(max(rank.sort_SPES_clin(:,4))));
colororder({'k','k'})
set(gca,'YTick',(1:size(rank.fig_sort_names_SPES_clin,1)),'YTickLabel',rank.name_cceps_SPES_clin)
yyaxis left
set(gca, 'YDir', 'reverse');
set(gca,'TickLength',[0 0])
ylim([1, max([size(rank.fig_sort_names_SPES_clin,1) size(rank.fig_sort_names_SPES_prop,1)])])
ylabel('order SPES-clin')

yyaxis right
set(gca,'YTick',(1:size(rank.fig_sort_names_SPES_prop,1)),'YTickLabel',rank.name_cceps_SPES_prop)
set(gca, 'YDir', 'reverse');
ylim([1, max([size(rank.fig_sort_names_SPES_clin,1) size(rank.fig_sort_names_SPES_prop,1)])])
ylabel('order SPES-prop')


xlim([1, 2])
set(gca,'xtick',[])
if PVAL_stmp < 0.01
    str_main = sprintf('%s, (p = <0.01)',ccep.sub_label);  
elseif PVAL_stmp < 0.05
    str_main = sprintf('%s, (p = <0.05)',ccep.sub_label);  
else
    str_main = sprintf('%s, (p = %1.2f)',ccep.sub_label, PVAL_stmp); 
end
sgtitle(str_main)

n=1;
for k = 1:length(rank.fig_sort_names_SPES_clin)
    if ~strcmp(rank.fig_sort_names_SPES_clin{k},' ')
        [~,loc2] = ismember(rank.fig_sort_names_SPES_clin{k}, rank.fig_sort_names_SPES_prop);
        
        line([1, 2],[k, loc2],'Color',cm(rank.sort_SPES_clin(n,4),:), 'LineWidth',2)  ;
        n=n+1;
    end
end


% Save figure
outlabel=sprintf('sub-%s_ranking.png',ccep.sub_label);
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/Ranking/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')

close all

%% Make horizontal bar graph for the number of CCEPs per stimulation piar
% Take the ranking of the stimulation pairs based on the clinincal-SPES
% separate per patient
dataBarPlot(:,1:2) = rank.SPES_clin(:,1:2);
dataBarPlot(:,4) = rank.SPES_clin(:,3) *(-1);   
dataBarPlot(:,3) = rank.SPES_prop(:,3);
dataBarPlot_sorted = sortrows(dataBarPlot,4,'descend');

% % When value is 0, add a little bit so you do see a very small line next to
% % the zero center line
% dataBarPlot_sorted(dataBarPlot_sorted(:,3) == 0,3) = dataBarPlot_sorted(dataBarPlot_sorted(:,3) == 0, 3) +0.1;
% dataBarPlot_sorted(dataBarPlot_sorted(:,4) == 0,4) = dataBarPlot_sorted(dataBarPlot_sorted(:,4) == 0, 4) -0.1;


figure('Position',[680,227,1076,826])
b1 = barh(dataBarPlot_sorted(:,4), 'r');     % Left --> Clinical

hold on
b2 = barh(dataBarPlot_sorted(:,3), 'b');         % Right --> propofol

% legend('prop','clin','Location','southeast')


% Now alter the ticks.
xticks = get(gca, 'xtick');

% Get the current labels
labels = get(gca, 'xtickLabel'); 

if ischar(labels)
    labels = cellstr(labels); 
end

% Figure out which ones we need to change
toscale = xticks < 0;

% Replace the text for the ones < 0
labels(toscale) = arrayfun(@(x)sprintf('%1.0f', x), ...
                           abs(xticks(toscale) * -1 ), 'uniformoutput', false);

% Update the tick locations and the labels
set(gca, 'xtick', xticks, 'xticklabel', labels)

xmax_r = max(get(gca, 'xlim'));
xmax_l = min(get(gca,'xlim'));
label(1) = text(xmax_r/2, -2.5, 'Propofol SPES','FontSize',13,'FontWeight','bold');
label(2) = text(xmax_l/2, -2.5, 'Clinical SPES','FontSize',13,'FontWeight','bold','HorizontalAlignment','center');

title(sprintf('Number of CCEPs detected per stimulation pair for %s',ccep.sub_label))

set(gca,'ycolor','none','XGrid','on','Box','off')

% Make little txt note in each bar that indicates the number of CCEPs
xtips1 = b1.YEndPoints -1.3;
ytips1 = b1.XEndPoints;
labels1 = arrayfun(@(x)sprintf('%1.0f', x), ...
                           abs(b1.YData * -1 ), 'uniformoutput', false);
% labels1 = string(b1.YData);
text(xtips1,ytips1,labels1,'VerticalAlignment','middle')

xtips2 = b2.YEndPoints +0.3;
ytips2 = b2.XEndPoints;
labels2 = string(b2.YData);
text(xtips2,ytips2,labels2,'VerticalAlignment','middle')


% Save figure
outlabel=sprintf('sub-%s_barGraph.png',ccep.sub_label);
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/BarGraph/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')





% Bar graph with all stimulation pairs of all patients



%% Spearman correlation
% For the indegree, outdegree and betweenness centrality per electrode
measure = {'indegree','outdegree','BC'};

for n = 1:size(measure,2)
    
    for i = 1:size(mode,2)
        
        if strcmp(mode{i},'SPES_clin')
            ERs = agreement_parameter.([measure{n} 'N_Clin']);
        elseif strcmp(mode{i},'SPES_prop')
            ERs = agreement_parameter.([measure{n} 'N_Prop']);
        end

        rank.([measure{n}, mode{i}])(:,1) = 1:size(ccep.n1_peak_sample,1);
        rank.([measure{n}, mode{i}])(:,2) = ERs;
        [~, order] = sort(rank.([measure{n} mode{i}])(:,2), 'descend','MissingPlacement','last');       % Use 'MissingPlacement','last', to place NaN's at the end of the ranking
        rank.(['sort_' measure{n} mode{i}]) = rank.([measure{n} mode{i}])(order, :);
        
        % If the next stimpair has the same number of ERs, give it the same order number
        rank.(['sort_' measure{n} mode{i}])(1,3) = 1;
        for j = 2:size(rank.(['sort_' measure{n} mode{i}]),1)
            if rank.(['sort_' measure{n}  mode{i}])(j,2) == rank.(['sort_' measure{n} mode{i}])(j-1,2)
                rank.(['sort_' measure{n}  mode{i}])(j,3) = rank.(['sort_' measure{n} mode{i}])(j-1,3);
            else
                rank.(['sort_' measure{n} mode{i}])(j,3) = j;
            end
        end
        
        rank.(['sort_names_' measure{n} mode{i}]) = ccep.ch(order);
        
         % Sorted matrix based on stimpair number, so ranking is unsorted 
        [~,order3] = sort(rank.(['sort_' measure{n} mode{i}])(:,1),'MissingPlacement','last');
        rank.(['unsort_' measure{n} mode{i}]) = rank.(['sort_' measure{n} mode{i}])(order3,:);        
    end
    
    % To determine the significance for the ranking
    % Use Spearman 'rows','pairwise' to ensure that row's with NaN's in both columns are not considered in the analysis.
    [RHO.(measure{n}), PVAL.(measure{n})] = corr(rank.(['unsort_' measure{n} mode{1}])(:,2), rank.(['unsort_', measure{n} mode{2}])(:,2),'Type','Spearman','rows','pairwise');
    fprintf('Spearman Corr (ranking) between %s per electrode of SPES-clin and SPES-prop gives, p-value = %1.4f, rho = %1.3f, for %s \n', measure{n}, PVAL.(measure{n}), RHO.(measure{n}), ccep.sub_label);

end

fprintf('------------ NEXT PATIENT ------------\n')

% Write to struct
statistics.p_BC = PVAL.BC;
statistics.rho_BC = RHO.BC;
statistics.p_indegree= PVAL.indegree;
statistics.rho_indegree = RHO.indegree;
statistics.p_outdegree = PVAL.outdegree;
statistics.rho_outdegree = RHO.outdegree;
statistics.p_stimp = PVAL_stmp;
statistics.rho_stimp = RHO_stmp;
statistics.p_ERsperStimp = p;  

clear PVAL_stmp
        
end       

