function statistics = statistical_agreement(agreement_parameter,ccep10)

% Mann-Whitney-Wilcoxon non parametric test for two paired /unpaired groups
% For the overall, positive and negative agreement





% Mann-Whitney-Wilcoxon non parametric test for two paired /unpaired groups
% For the number of ERs detected in the 2 stimuli and the 10 stimuli
% signals.
ER_stimp10 = agreement_parameter.ERs_stimp10;
ER_stimp10_doub = num2cell(ER_stimp10);
ER_stimp2 = agreement_parameter.ERs_stimp2;
ER_stimp2_doub = num2cell(ER_stimp2);

NorDis10 = lillietest(ER_stimp10);                  % null hypothesis that x is normally distributed, results in 1 when the null hypothesis is rejected 
NorDis2 = lillietest(ER_stimp2);

figure(1)
subplot(2,1,1)
normplot(ER_stimp10)
subplot(2,1,2)
normplot(ER_stimp2)

% Statistics for NUMBER of detected ERS
[rho,p] = corr(ER_stimp10, ER_stimp2, 'Type', 'Spearman');


for stimp = 1:size(ER_stimp10,1)
  stimsets{stimp} = unique([ccep10.stimsets_avg(stimp,1), ccep10.stimsets_avg(stimp,2)]);
end
stimsets = stimsets';

rank10(:,1) = stimsets;
rank10(:,2) = ER_stimp10_doub;
[~, order] = sort([rank10{:,2}], 'descend');
sort_rank10 = rank10(order, :);
sorted_ranking10_SP = sort_rank10(:,1);

rank2(:,1) = stimsets;
rank2(:,2) = ER_stimp2_doub;
[~, order2] = sort([rank2{:,2}], 'descend');
sort_rank2 = rank2(order2, :);
sorted_ranking2_SP = [sort_rank2(:,1)];




% Statistics for STIMPAIR ranking
[rho_sp,p_sp] = corr(sorted_ranking10_SP, sorted_ranking2_SP, 'Type', 'Spearman');


% Rank biserial correlation 
% For the indegree and outdegree







% Rank biserial correlation 
% For the Betweenness correlation





end

