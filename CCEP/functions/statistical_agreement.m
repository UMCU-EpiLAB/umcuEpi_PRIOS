function statistics = statistical_agreement(agreement_parameter,ccep10)

% Mann-Whitney-Wilcoxon non parametric test for two paired /unpaired groups
% For the overall, positive and negative agreement





%% Spearman correlation for the number of ERs detected in the 2 stimuli and the 10 stimuli
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

% Statistics for NUMBER of detected ERS
[rho,p] = corr(ER_stimp10, ER_stimp2, 'Type', 'Spearman');          % not normally distributed
%[rho,p] = corr(ER_stimp10, ER_stimp2, 'Type', 'Pearson');          % normally distributed



%% Spearman test for ranking of stimpairs

for stimp = 1:size(ER_stimp10,1)
  stimsets{stimp} = unique([ccep10.stimsets_avg(stimp,1), ccep10.stimsets_avg(stimp,2)]);
end
stimsets = stimsets';

% For all stimuli
rank10(:,1) = stimsets;
rank10(:,2) = num2cell(ER_stimp10);
[sortedValue, order] = sort([rank10{:,2}], 'descend');          % most ER evoking stimpairs first
sort_rank10 = rank10(order, :) ;
[~, ~, ic] = unique(sortedValue);
Ranked = num2cell(max(ic)-ic+1);
sort_rank10(:,3) = Ranked;                                      % Assign label based on the ranking (1 is highest)

% For 2 stimuli
rank2(:,1) = stimsets;
rank2(:,2) = num2cell(ER_stimp2);
[sortedValue_2, order2] = sort([rank2{:,2}], 'descend');
sort_rank2 = rank2(order2, :);
[~, ~, ic2] = unique(sortedValue_2);
Ranked2 = num2cell(max(ic2)-ic2+1);
sort_rank2(:,3) = Ranked2;

% Concatenate/merge electrodes of the stimulation pair to be able to sort them
for i = 1: length(sort_rank10)
    stimpair_concat(i,:) = sscanf(sprintf('%d%d,',[sort_rank10{i, 1}(1).';sort_rank10{i, 1}(2)]),'%d,'); %horzcat(sort_rank10{i, 1}(1),sort_rank10{i, 1}(2));
    stimpair_concat2(i,:) = sscanf(sprintf('%d%d,',[sort_rank2{i, 1}(1).';sort_rank2{i, 1}(2)]),'%d,'); %horzcat(sort_rank2{i, 1}(1),sort_rank2{i, 1}(2));
end

stimpair_concat(:,2) = cell2mat(sort_rank10(:,3));
Stimpair_sorted = sortrows([stimpair_concat], 'ascend');

stimpair_concat2(:,2) = cell2mat(sort_rank2(:,3));
Stimpair_sorted2 = sortrows([stimpair_concat2], 'ascend');

Statistic_mat(:,1:2) = Stimpair_sorted(:,1:2);
Statistic_mat(:,3) = Stimpair_sorted2(:,2);
% colNames = {'Stimpair', 'Ranking10', 'Ranking2'};
% T = array2table(Statistic_mat,'VariableNames', colNames);    % able to check           

[RHO,PVAL] = corr(Statistic_mat(:,2) , Statistic_mat(:,3) ,'Type','Spearman');


SubjectName = extractBetween(ccep10.dataName,'ieeg/sub-','_ses-1_');
fprintf('Spearman Corr between stimpair ranking of 10 and 2 stimuli gives, p-value = %1.4f, rho = %1.3f, for %s \n', PVAL, RHO, SubjectName{1});



%% Rank biserial correlation 
% For the indegree and outdegree per electrode

% For 10 stimuli
indegree10_elec(:,1) = ccep10.ch;
indegree10_elec(:,2) = num2cell(agreement_parameter.indegreeN_10') ;

[sortedInd10, orderInd10] = sort([indegree10_elec{:,2}], 'descend');            % rank based on highst value first
Ind10 = indegree10_elec(orderInd10, :);
[~, ~, ic_ind10] = unique(sortedInd10);
Ranked_ind10 = num2cell(max(ic_ind10)-ic_ind10+1);
Ind10(:,3) = Ranked_ind10;
sort_Ind10 = sortrows([Ind10], 'ascend');
rank_sort_Ind10 = [sort_Ind10{:,3}]';

% For 2 stimuli
indegree2_elec(:,1) = ccep10.ch;
indegree2_elec(:,2) = num2cell(agreement_parameter.indegreeN_2') ;

[sortedInd2, orderInd2] = sort([indegree2_elec{:,2}], 'descend');            % rank based on highst value first
Ind2 = indegree2_elec(orderInd2, :);
[~, ~, ic_ind2] = unique(sortedInd2);
Ranked_ind2 = num2cell(max(ic_ind2)-ic_ind2+1);
Ind2(:,3) = Ranked_ind2;
sort_Ind2 = sortrows([Ind2], 'ascend');
rank_sort_Ind2 = [sort_Ind2{:,3}]';

[RHO_ind, PVAL_ind] = corr(rank_sort_Ind10 , rank_sort_Ind2 ,'Type','Spearman');

fprintf('Spearman Corr between stimpair ranking of 10 and 2 stimuli gives, p-value = %1.4f, rho = %1.3f, for %s \n', PVAL_ind, RHO_ind, SubjectName{1});





%% Rank biserial correlation 
% For the Betweenness correlation





end

