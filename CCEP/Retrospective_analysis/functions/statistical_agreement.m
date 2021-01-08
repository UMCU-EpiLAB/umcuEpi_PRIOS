function [statistics,rank] = statistical_agreement(myDataPath, agreement_parameter,ccep10)
%% Wilcoxon signed rank non parametric test for two paired groups
% For the number of ERs evoked in 2 and 10 stimuli
% null hypothesis is that the two means/medians are the same
SubjectName = extractBetween(ccep10.dataName,'ieeg/sub-','_ses-1_');

ER_stimp10 = agreement_parameter.ERs_stimp10;
ER_stimp2 = agreement_parameter.ERs_stimp2;

NorDis10 = lillietest(ER_stimp10);                  % null hypothesis that x is normally distributed, results in 1 when the null hypothesis is rejected
NorDis2 = lillietest(ER_stimp2);

% Check for monotonic relation
figure(2)
scatter(ER_stimp10,ER_stimp2)
refline

% Check for normal distribution
figure(1)
subplot(2,1,1)
normplot(ER_stimp10)                                % normal distribution is datapoints along the reference line
subplot(2,1,2)
normplot(ER_stimp2)

p = signrank(ER_stimp10, ER_stimp2) ;           % tests the null hypothesis that data in x and y are samples from continuous distributions with equal medians

if p<0.05
    fprintf('Test between the number of ERs evoked in 10 and in 2 stims per stimulation pair gives p-value = %1.4f. This means that there is a significant difference between stimulating 10 or 2 times for %s \n', p, SubjectName{1});
else
    fprintf('Test between the number of ERs evoked in 10 and in 2 stims per stimulation pair gives p-value = %1.4f. This means that there is NO significant difference between stimulating 10 or 2 times for %s \n', p, SubjectName{1});
end

%% Spearman test for ranking of stimpairs
% Rank the stimuli from stimulus pairs with most evoked ERs to stimulus
% pairs with least evoked ERs
mode = {'stims10','stims2'};
rank = struct;

for i=1:size(mode,2)
    
    if strcmp(mode{i},'stims10')
        ERs = ER_stimp10;
    elseif strcmp(mode{i},'stims2')
        ERs = ER_stimp2;
    end
    
    rank.(mode{i})(:,1:2) = ccep10.stimsets_avg;
    rank.(mode{i})(:,3) = ERs;
    [~, order] = sort(rank.(mode{i})(:,3), 'descend','MissingPlacement','last');          % most ER evoking stimpairs first, place NaN's at the end of the ranking)
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
    
    rank.(['sort_names_' mode{i}]) = ccep10.stimpnames_avg(order);
    
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

% Sorted matrix based on stimpair number, so ranking isunsorted 
for i=1:size(mode,2)
    [~,order] = sortrows(rank.(['sort_' mode{i}])(:,1:2));
    rank.(['unsort_' mode{i}]) = rank.(['sort_' mode{i}])(order,:);         
end

% If the order of the stimulation pairs is not equal, than the ranking
% cannot be compared
if ~isequal(rank.unsort_stims10(:,1:2),rank.unsort_stims2(:,1:2))
    error('Sorting stimulus pairs is incorrect and led to unequal matrices in 10 and 2 stims')
end

% Test the hypothesis of NO correlation
% When p <0.05, an rho is close to (-)1, rejection of the hypothesis that
% no correlation exists between the two columns i.e. there is a
% correlation
[RHO_stmp,PVAL_stmp] = corr(rank.unsort_stims10(:,4) , rank.unsort_stims2(:,4) ,'Type','Spearman');            % Test the hypothesis that the correlation is NOT 0
fprintf('Spearman Corr between stimpair ranking of 10 stims and 2 stims gives, p-value = %1.4f, rho = %1.3f, for %s \n', PVAL_stmp, RHO_stmp, SubjectName{1});

figure('Position',[1074,4,519,1052]);
cm = colormap(parula(max(rank.sort_stims10(:,4))));
colororder({'k','k'})
set(gca,'YTick',(1:size(rank.fig_sort_names_stims10,1)),'YTickLabel',rank.fig_sort_names_stims10)
yyaxis left
set(gca, 'YDir', 'reverse');
set(gca,'TickLength',[0 0])
ylim([1, max([size(rank.fig_sort_names_stims10,1) size(rank.fig_sort_names_stims2,1)])])
ylabel('order 10 stims')

yyaxis right
set(gca,'YTick',(1:size(rank.fig_sort_names_stims2,1)),'YTickLabel',rank.fig_sort_names_stims2)
set(gca, 'YDir', 'reverse');
ylim([1, max([size(rank.fig_sort_names_stims10,1) size(rank.fig_sort_names_stims2,1)])])
ylabel('order 2 stims')

xlim([1, 2])
set(gca,'xtick',[])
if PVAL_stmp(:) <0.01
       str_main = sprintf('%s, p < 0.01',SubjectName{1});
    elseif PVAL_stmp(:) <0.05
        str_main = sprintf('%s, p < 0.05',SubjectName{1});
    else
        str_main = sprintf('sub-%s, p = %1.2f', SubjectName{1},PVAL_stmp);
    end
sgtitle(str_main)

n=1;
for k = 1:length(rank.fig_sort_names_stims10)
    if ~strcmp(rank.fig_sort_names_stims10{k},' ')
        [~,loc2] = ismember(rank.fig_sort_names_stims10{k}, rank.fig_sort_names_stims2);
        
        line([1, 2],[k, loc2],'Color',cm(rank.sort_stims10(n,4),:), 'LineWidth',2)  ;
        n=n+1;
    end
end

% Save figure
outlabel=sprintf('sub-%s_ranking.jpg',SubjectName{1});
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/Ranking/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'jpg')


%% Spearman correlation
% For the indegree, outdegree and betweenness centrality per electrode
measure = {'indegree','outdegree','BC'};
for n=1:size(measure,2)
    
    for i=1:size(mode,2)
        
        if strcmp(mode{i},'stims10')
            ERs = agreement_parameter.([measure{n} 'N_10']);
        elseif strcmp(mode{i},'stims2')
            ERs = agreement_parameter.([measure{n} 'N_2']);
        end
        
        rank.([measure{n}, mode{i}])(:,1) = 1:size(ccep10.ch,1);
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
        
        rank.(['sort_names_' measure{n} mode{i}]) = ccep10.ch(order);
        
        % Sorted matrix based on stimpair number, so ranking isunsorted 
        [~,order] = sort(rank.(['sort_' measure{n} mode{i}])(:,1),'MissingPlacement','last');
        rank.(['unsort_' measure{n} mode{i}]) = rank.(['sort_' measure{n} mode{i}])(order,:);        
    end
    
    % To determine the significance for the ranking
    [RHO.(measure{n}), PVAL.(measure{n})] = corr(rank.(['unsort_' measure{n} mode{1}])(:,3), rank.(['unsort_', measure{n} mode{2}])(:,3),'Type','Spearman','rows','pairwise');
    
    % To determine the significance for the absolute value.
    pval.(measure{n}) = signrank(rank.(['unsort_' measure{n} mode{1}])(:,2), rank.(['unsort_', measure{n} mode{2}])(:,2));

    fprintf('Spearman Corr (ranking) between %s per electrode of 10 and 2 stimuli gives, p-value = %1.4f, rho = %1.3f, for %s \n', measure{n}, PVAL.(measure{n}), RHO.(measure{n}), SubjectName{1});
    fprintf('Wilcoxon (abs values) between %s per electrode of SPES-clin and SPES-prop gives, p-value = %1.4f, for %s \n', measure{n}, pval.(measure{n}), SubjectName{1});

end

fprintf('------------ NEXT PATIENT ------------\n')

 % Write to variable
statistics.p_BC = PVAL.BC;
statistics.rho_BC = RHO.BC;
statistics.p_indegree= PVAL.indegree;
statistics.rho_indegree = RHO.indegree;
statistics.p_outdegree = PVAL.outdegree;
statistics.rho_outdegree = RHO.outdegree;
statistics.p_stimp = PVAL_stmp;
statistics.rho_stimp = RHO_stmp;
statistics.p_ERsperStimp = p;
% statistics.ranking2stimp = sort_rank2;
% statistics.ranking10stimp = sort_rank10;

statistics.p_indegree_abs = pval.indegree;
statistics.p_outdegree_abs = pval.outdegree;
statistics.p_BC_abs = pval.BC;

end

