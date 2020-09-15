function statistics = statistical_agreement(myDataPath, agreement_parameter,ccep10)
%% Paired t-test or Wilcoxon signed rank non parametric test for two paired groups
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

if NorDis10 == 1 && NorDis2 ==1
    p = signrank(ER_stimp10, ER_stimp2) ;           % tests the null hypothesis that data in x and y are samples from continuous distributions with equal medians
    fprintf('The detected ERs per stimulation pair is NOT normally distributed, Wilcoxon Signed Rank test is used.\n')
else
    fprintf('The detected ERs per stimulation pair is normally distributed, Paired T-test is used.\n')
    [~,p] = ttest(ER_stimp10, ER_stimp2);          % alpha default = 0.05
end

if p<0.05
    fprintf('Test between the number of ERs evoked in 10 and in 2 stims per stimulation pair gives p-value = %1.4f. This means that there is a significant difference between stimulating 10 or 2 times for %s \n', p, SubjectName{1});
else
    fprintf('Test between the number of ERs evoked in 10 and in 2 stims per stimulation pair gives p-value = %1.4f. This means that there is NO significant difference between stimulating 10 or 2 times for %s \n', p, SubjectName{1});
end

%% Spearman test for ranking of stimpairs
% % pre-allocation
% stimsets = cell(1,length(ccep10.stimsets_avg));
% for stimp = 1:size(ER_stimp10,1)
%   stimsets{stimp} = unique([ccep10.stimsets_avg(stimp,1), ccep10.stimsets_avg(stimp,2)]);
% end
% stimsets = stimsets';
%
% for i = 1:length(stimsets)
%     stimnames(i,:) = {[ccep10.ch{stimsets{i,:}}]} ;
% end
%

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
    [~, order] = sort(rank.(mode{i})(:,3), 'descend');          % most ER evoking stimpairs first
    rank.(['sort_' mode{i}]) = rank.(mode{i})(order, :);
    
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

% rank10(:,1) = stimsets;
% rank10(:,2) = num2cell(ER_stimp10);
% [sortedValue_10, order10] = sort([rank10{:,2}], 'descend');
% sort_rank10 = rank10(order10, :);
% [~, ~, ic10] = unique(sortedValue_10);
% Ranked10 = num2cell(max(ic10)-ic10+1);
% sort_rank10(:,3) = Ranked10;
% order10 = cell2mat(Ranked10);
% Unique_rank10 = num2cell(unique(sort(ic10(:),'ascend')));
%
% for i = 1:length(sort_rank10)
%     names_ranked10(i,:) ={[ccep10.ch{sort_rank10{i,1}}]}  ;
% end

% for i = 1:length(Unique_rank10)
%
%     rankingnum10 = Unique_rank10{i};
%     rowsWithRanknum10 = find(order10 == rankingnum10);
%     stimpWithRank10 = sort_rank10(rowsWithRanknum10,1)';
%     nmbrOfVal10 = numel(stimpWithRank10);
%     stim_per_rank10(i,1) = {[stimpWithRank10{1, 1:nmbrOfVal10};]};
%     stimnames_rank10(i,1) = {[ccep10.ch{stim_per_rank10{i,1}}]};
% end
% stim_per_rank10(:,2) = Unique_rank10;
%
% For 2 stimuli
% rank2(:,1) = stimsets;
% rank2(:,2) = num2cell(ER_stimp2);
% [sortedValue_2, order2] = sort([rank2{:,2}], 'descend');
% sort_rank2 = rank2(order2, :);
% [~, ~, ic2] = unique(sortedValue_2);
% Ranked2 = num2cell(max(ic2)-ic2+1);
% sort_rank2(:,3) = Ranked2;
% order2 = cell2mat(Ranked2);
% Unique_rank2 = num2cell(unique(sort(ic2(:),'ascend')));
%
% for i = 1:length(sort_rank2)
%     names_ranked2(i,:) ={[ccep10.ch{sort_rank2{i,1}}]}  ;
% end
%
% for i = 1:length(Unique_rank2)
%
%     rankingnum2 = Unique_rank2{i};
%     rowsWithRanknum2 = find(order2 == rankingnum2);
%     stimpWithRank2 = sort_rank2(rowsWithRanknum2,1)';
%     nmbrOfVal2 = numel(stimpWithRank2);
%     stim_per_rank2(i,1) = {[stimpWithRank2{1, 1:nmbrOfVal2};]};
%     stimnames_rank2(i,1) = {[ccep10.ch{stim_per_rank2{i,1}}]};
%
% end
%     stim_per_rank2(:,2) = Unique_rank2;
%
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

cm = colormap(parula(max(rank.sort_stims10(:,4))));
figure('Position',[1074,4,519,1052]);
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
str_main = sprintf('sub-%s', SubjectName{1});
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

% Plot the same figure as above though now the stimpairs with the same
% number of ERs are on the same row instead of alphabetical order.
% figure('Position',[694,4,1224,1052]);
% colororder({'k','k'})
% ylim([1, numel(stimnames_rank10)])
%
% if numel(stimnames_rank10) > numel(stimnames_rank2)
%    yyaxis right
%    ylim([1, numel(stimnames_rank10)])
%    cm = colormap(lines(max(ic)));
% else %numel(stimnames_rank2) => numel(stimnames_rank10)
%     yyaxis left
%     ylim([0.5, numel(stimnames_rank2)])
%     yyaxis right
%     ylim([0.5, numel(stimnames_rank2)])
%     cm = colormap(lines(max(ic2)));
% end
%
% yyaxis left
% set(gca,'YTick',(1:numel(stimnames_rank10)),'YTickLabel',stimnames_rank10)
% set(gca, 'YDir', 'reverse');
% ylabel('order 10 stims')
%
% yyaxis right
% set(gca,'YTick',(1:numel(stimnames_rank2)),'YTickLabel',stimnames_rank2)
% set(gca, 'YDir', 'reverse');
%
% xlim([1, 2])
% set(gca,'xtick',[])
% str_main = sprintf('sub-%s', SubjectName{1});
% sgtitle(str_main)
% ylabel('order 2 stims')
%
% for k = 1:length(sort_rank10)
%     name10 = names_ranked10{k,1} ;
%     [~,loc2] = ismember(name10, names_ranked2(:,1));
%
%     line([1, 2],[sort_rank10{k, 3}, sort_rank2{loc2, 3}],'Color',cm(ic2(k),:), 'LineWidth',2)  ;                 % [x,x],[rank10, rank2]
%     hold on
% end
%
% % Save figure
% outlabel=sprintf('sub-%s_rankOrdered.jpg',subj{1});
% path = [fullfile(myDataPath.CCEPpath,'Visualise_agreement/Ranking/')];
% if ~exist(path, 'dir')
%    mkdir(path);
% end
% saveas(gca,[path,outlabel],'jpg')


%% Concatenate/merge electrodes of the stimulation pair to be able to sort them
% for i = 1: length(sort_rank10)
%     stimpair_concat(i,:) = sscanf(sprintf('%d%d,',[sort_rank10{i, 1}(1).';sort_rank10{i, 1}(2)]),'%d,'); %horzcat(sort_rank10{i, 1}(1),sort_rank10{i, 1}(2));
%     stimpair_concat2(i,:) = sscanf(sprintf('%d%d,',[sort_rank2{i, 1}(1).';sort_rank2{i, 1}(2)]),'%d,'); %horzcat(sort_rank2{i, 1}(1),sort_rank2{i, 1}(2));
% end
%
% stimpair_concat(:,2) = cell2mat(sort_rank10(:,3));
% Stimpair_sorted = sortrows([stimpair_concat], 'ascend');
%
% stimpair_concat2(:,2) = cell2mat(sort_rank2(:,3));
% Stimpair_sorted2 = sortrows([stimpair_concat2], 'ascend');
%
% Statistic_mat(:,1:2) = Stimpair_sorted(:,1:2);
% Statistic_mat(:,3) = Stimpair_sorted2(:,2);
% colNames = {'Stimpair', 'Ranking10', 'Ranking2'};
% T = array2table(Statistic_mat,'VariableNames', colNames);    % able to check

for i=1:size(mode,2)
    [~,order] = sortrows(rank.(['sort_' mode{i}])(:,1:2));
    rank.(['unsort_' mode{i}]) = rank.(['sort_' mode{i}])(order,:);
end

if ~isequal(rank.unsort_stims10(:,1:2),rank.unsort_stims2(:,1:2))
    error('Sorting stimulus pairs is incorrect and led to unequal matrices in stims10 and stims2')
end

% Test the hypothesis of NO correlation
% When p <0.05, an rho is close to (-)1, rejection of the hypothesis that no correlation exists between the two columns
% [RHO_stmp,PVAL_stmp] = corr(Statistic_mat(:,2) , Statistic_mat(:,3) ,'Type','Spearman');            % Test the hypothesis that the correlation is NOT 0
[RHO_stmp,PVAL_stmp] = corr(rank.unsort_stims10(:,4) , rank.unsort_stims2(:,4) ,'Type','Spearman');            % Test the hypothesis that the correlation is NOT 0
fprintf('Spearman Corr between stimpair ranking of 10 and 2 stimuli gives, p-value = %1.4f, rho = %1.3f, for %s \n', PVAL_stmp, RHO_stmp, SubjectName{1});

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
        [~, order] = sort(rank.([measure{n} mode{i}])(:,2), 'descend');
        rank.(['sort_' measure{n} mode{i}]) = rank.([measure{n} mode{i}])(order, :);
        
        rank.(['sort_' measure{n} mode{i}])(1,3) = 1;
        for j = 2:size(rank.(['sort_' measure{n} mode{i}]),1)
            if rank.(['sort_' measure{n}  mode{i}])(j,2) == rank.(['sort_' measure{n} mode{i}])(j-1,2)
                rank.(['sort_' measure{n}  mode{i}])(j,3) = rank.(['sort_' measure{n} mode{i}])(j-1,3);
            else
                rank.(['sort_' measure{n} mode{i}])(j,3) = j;
            end
        end
        
        rank.(['sort_names_' measure{n} mode{i}]) = ccep10.ch(order);
        
        [~,order] = sort(rank.(['sort_' measure{n} mode{i}])(:,1));
        rank.(['unsort_' measure{n} mode{i}]) = rank.(['sort_' measure{n} mode{i}])(order,:);
        
    end
    
    
    
    % % For 10 stimuli
    % indegree10_elec(:,1) = ccep10.ch;
    % indegree10_elec(:,2) = num2cell(agreement_parameter.indegreeN_10') ;
    %
    % [sortedInd10, orderInd10] = sort([indegree10_elec{:,2}], 'descend');            % rank based on highst value first
    % Ind10 = indegree10_elec(orderInd10, :);
    % [~, ~, ic_ind10] = unique(sortedInd10);
    % Ranked_ind10 = num2cell(max(ic_ind10)-ic_ind10+1);
    % Ind10(:,3) = Ranked_ind10;
    % sort_Ind10 = sortrows([Ind10], 'ascend');
    % rank_sort_Ind10 = [sort_Ind10{:,3}]';
    %
    % % For 2 stimuli
    % indegree2_elec(:,1) = ccep10.ch;
    % indegree2_elec(:,2) = num2cell(agreement_parameter.indegreeN_2') ;
    %
    % [sortedInd2, orderInd2] = sort([indegree2_elec{:,2}], 'descend');            % rank based on highst value first
    % Ind2 = indegree2_elec(orderInd2, :);
    % [~, ~, ic_ind2] = unique(sortedInd2);
    % Ranked_ind2 = num2cell(max(ic_ind2)-ic_ind2+1);
    % Ind2(:,3) = Ranked_ind2;
    % sort_Ind2 = sortrows([Ind2], 'ascend');
    % rank_sort_Ind2 = [sort_Ind2{:,3}]';
    %
    %
    % % Check
    % if ~ismember(sort_Ind2(:,1),sort_Ind10(:,1))
    %     fprintf('The stimulation pair order of the two situations is not equal. \n');
    % end
    
    % [RHO_ind, PVAL_ind] = corr(rank_sort_Ind10 , rank_sort_Ind2 ,'Type','Spearman');
    [RHO.(measure{n}), PVAL.(measure{n})] = corr(rank.(['unsort_' measure{n} mode{1}])(:,3), rank.(['unsort_', measure{n} mode{2}])(:,3),'Type','Spearman','rows','pairwise');
    
    fprintf('Spearman Corr between %s per electrode of 10 and 2 stimuli gives, p-value = %1.4f, rho = %1.3f, for %s \n', measure{n}, PVAL.(measure{n}), RHO.(measure{n}), SubjectName{1});
end

fprintf('------------ NEXT PATIENT ------------\n')

%% Rank biserial correlation
% For the outdegree per electrode
% For 10 stimuli
% outdegree10_elec(:,1) = ccep10.ch;
% outdegree10_elec(:,2) = num2cell(agreement_parameter.outdegreeN_10') ;
%
% [sortedOutd10, orderOutd10] = sort([outdegree10_elec{:,2}], 'descend');            % rank based on highst value first
% Outd10 = outdegree10_elec(orderOutd10, :);
% [~, ~, ic_outd10] = unique(sortedOutd10);
% Ranked_outd10 = num2cell(max(ic_outd10)-ic_outd10+1);
% Outd10(:,3) = Ranked_outd10;
% sort_Outd10 = sortrows([Outd10], 'ascend') ;
% rank_sort_Outd10 = [sort_Outd10{:,3}]';
%
% % For 2 stimuli
% outdegree2_elec(:,1) = ccep10.ch;
% outdegree2_elec(:,2) = num2cell(agreement_parameter.outdegreeN_2') ;
%
% [sortedOutd2, orderOutd2] = sort([outdegree2_elec{:,2}], 'descend');            % rank based on highst value first
% Outd2 = outdegree2_elec(orderOutd2, :);
% [~, ~, ic_outd2] = unique(sortedOutd2);
% Ranked_outd2 = num2cell(max(ic_outd2)-ic_outd2+1);
% Outd2(:,3) = Ranked_outd2;
% sort_Outd2 = sortrows([Outd2], 'ascend');
% rank_sort_Outd2 = [sort_Outd2{:,3}]';
%
%
% % Check
% if ~ismember(sort_Outd2(:,1),sort_Outd10(:,1))
%     fprintf('The stimulation pair order of the two situations is not equal. \n');
% end
%
% [RHO_outd, PVAL_outd] = corr(rank_sort_Outd10 , rank_sort_Outd2 ,'Type','Spearman');
%
% fprintf('Spearman Corr between outdegree per electrode of 10 and 2 stimuli gives, p-value = %1.4f, rho = %1.3f, for %s \n', PVAL_outd, RHO_outd, SubjectName{1});
%
%
%% Spearman correlation
% For the Betweenness correlation
% For 10 stimuli
% BC10(:,1) = ccep10.ch;
% BC10(:,2) = num2cell(agreement_parameter.BCN_all_10') ;
%
% % remove NaN's
% remove_NaN = cellfun(@isnan,BC10(:,2),'UniformOutput',false);
% tot_NaN = sum([remove_NaN{:}]);
% for i = 1:length(remove_NaN)
%     if ismember(remove_NaN{i,1}, 1)
%         BC10{i,2} = 0;
%     end
% end
% fprintf('%d NaNs are replaced with 0 for 10 stims, %s \n', tot_NaN, SubjectName{1});
%
%
% [sortedBC10, orderBC10] = sort([BC10{:,2}], 'descend');            % rank based on highst value first
% BC10 = BC10(orderBC10, :);
% [~, ~, ic_BC10] = unique(sortedBC10);
% Ranked_BC10 = num2cell(max(ic_BC10)-ic_BC10+1);
% BC10(:,3) = Ranked_BC10;
% sort_BC10 = sortrows([BC10], 'ascend') ;
% rank_sort_BC10 = [sort_BC10{:,3}]';
%
%
% % For 2 stimuli
% BC2(:,1) = ccep10.ch;
% BC2(:,2) = num2cell(agreement_parameter.BCN_2') ;
%
%
% % remove NaN's
% remove_NaN_2 = cellfun(@isnan,BC2(:,2),'UniformOutput',false);
% tot_NaN_2 = sum([remove_NaN_2{:}]);
% for i = 1:length(remove_NaN_2)
%     if ismember(remove_NaN_2{i,1}, 1)
%         BC2{i,2} = 0;
%     end
% end
% fprintf('%d NaNs are replaced with 0 for 2 stims, %s \n', tot_NaN_2, SubjectName{1});
%
%
% [sortedBC2, orderBC2] = sort([BC2{:,2}], 'descend');            % rank based on highst value first
% BC2 = BC2(orderBC2, :);
% [~, ~, ic_BC2] = unique(sortedBC2);
% Ranked_BC2 = num2cell(max(ic_BC2)-ic_BC2+1);
% BC2(:,3) = Ranked_BC2;
% sort_BC2 = sortrows([BC2], 'ascend') ;
% rank_sort_BC2 = [sort_BC2{:,3}]';
%
%
% % Check
% if ~ismember(sort_BC10(:,1),sort_BC2(:,1))
%     fprintf('The stimulation pair order of the two situations is not equal. \n');
% end
%
% [RHO_BC, PVAL_BC] = corr(rank_sort_BC10 , rank_sort_BC2 ,'Type','Spearman');
%
% fprintf('Spearman Corr between outdegree per electrode of 10 and 2 stimuli gives, p-value = %1.4f, rho = %1.3f, for %s \n\n\n', PVAL_BC, RHO_BC, SubjectName{1});
%
% % Write to variable
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

end

