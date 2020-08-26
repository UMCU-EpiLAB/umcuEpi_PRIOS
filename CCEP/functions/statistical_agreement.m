function statistics = statistical_agreement(agreement_parameter,ccep10)
%% Mann-Whitney-Wilcoxon non parametric test for two paired /unpaired groups
% For the number of ERs detected in the 2 stimuli and the 10 stimuli
% null hypothesis is that the two means are the same
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
    [p,h,stats] = ranksum(ER_stimp10, ER_stimp2) ;           % tests the null hypothesis that data in x and y are samples from continuous distributions with equal medians
else
    fprintf('The detected ERs per stimulation pair is normally distributed, Paired T-test is used')
    [h,p,ci,stats] = ttest(ER_stimp10, ER_stimp2);          % alpha default = 0.05
     
end

fprintf('Mann Whitney test between the ERs in 10 and ERs in 2 stims per stimulation pair gives p-value = %1.4f, reject the null-hypothesis? = %d, for %s \n', p, h, SubjectName{1});

%% Spearman test for ranking of stimpairs
stimsets = cell(1,length(ccep10.stimsets_avg));
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

[RHO_stmp,PVAL_stmp] = corr(Statistic_mat(:,2) , Statistic_mat(:,3) ,'Type','Spearman');
fprintf('Spearman Corr between stimpair ranking of 10 and 2 stimuli gives, p-value = %1.4f, rho = %1.3f, for %s \n', PVAL_stmp, RHO_stmp, SubjectName{1});

clearvars -except Statistic_mat RHO_stmp PVAL_stmp ccep10 agreement_parameter SubjectName


%% Spearman correlation 
% For the indegree per electrode

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


% Check
if ~ismember(sort_Ind2(:,1),sort_Ind10(:,1))
    fprintf('The stimulation pair order of the two situations is not equal. \n');
end

[RHO_ind, PVAL_ind] = corr(rank_sort_Ind10 , rank_sort_Ind2 ,'Type','Spearman');

fprintf('Spearman Corr between indegree per electrode of 10 and 2 stimuli gives, p-value = %1.4f, rho = %1.3f, for %s \n', PVAL_ind, RHO_ind, SubjectName{1});


clearvars -except Statistic_mat RHO_stmp PVAL_stmp ccep10 agreement_parameter SubjectName sort_Ind10 sort_Ind2 RHO_ind PVal_ind

%% Rank biserial correlation 
% For the outdegree per electrode
% For 10 stimuli
outdegree10_elec(:,1) = ccep10.ch;
outdegree10_elec(:,2) = num2cell(agreement_parameter.outdegreeN_10') ;

[sortedOutd10, orderOutd10] = sort([outdegree10_elec{:,2}], 'descend');            % rank based on highst value first
Outd10 = outdegree10_elec(orderOutd10, :);
[~, ~, ic_outd10] = unique(sortedOutd10);
Ranked_outd10 = num2cell(max(ic_outd10)-ic_outd10+1);
Outd10(:,3) = Ranked_outd10;
sort_Outd10 = sortrows([Outd10], 'ascend') ;
rank_sort_Outd10 = [sort_Outd10{:,3}]';

% For 2 stimuli
outdegree2_elec(:,1) = ccep10.ch;
outdegree2_elec(:,2) = num2cell(agreement_parameter.outdegreeN_2') ;

[sortedOutd2, orderOutd2] = sort([outdegree2_elec{:,2}], 'descend');            % rank based on highst value first
Outd2 = outdegree2_elec(orderOutd2, :);
[~, ~, ic_outd2] = unique(sortedOutd2);
Ranked_outd2 = num2cell(max(ic_outd2)-ic_outd2+1);
Outd2(:,3) = Ranked_outd2;
sort_Outd2 = sortrows([Outd2], 'ascend');
rank_sort_Outd2 = [sort_Outd2{:,3}]';


% Check
if ~ismember(sort_Outd2(:,1),sort_Outd10(:,1))
    fprintf('The stimulation pair order of the two situations is not equal. \n');
end

[RHO_outd, PVAL_outd] = corr(rank_sort_Outd10 , rank_sort_Outd2 ,'Type','Spearman');

fprintf('Spearman Corr between outdegree per electrode of 10 and 2 stimuli gives, p-value = %1.4f, rho = %1.3f, for %s \n', PVAL_outd, RHO_outd, SubjectName{1});

clearvars -except Statistic_mat RHO_stmp PVAL_stmp ccep10 agreement_parameter SubjectName sort_Ind10 sort_Ind2 RHO_ind PVal_ind sort_Outd10 sort_Outd2 RHO_outd PVal_outd


%% Spearman correlation 
% For the Betweenness correlation
% For 10 stimuli
BC10(:,1) = ccep10.ch;
BC10(:,2) = num2cell(agreement_parameter.BCN_all_10') ;

% remove NaN's
remove_NaN = cellfun(@isnan,BC10(:,2),'UniformOutput',false);
tot_NaN = sum([remove_NaN{:}]);
for i = 1:length(remove_NaN)
    if ismember(remove_NaN{i,1}, 1) 
        BC10{i,2} = 0;
    end
end
fprintf('%d NaNs are replaced with 0 for 10 stims, %s \n', tot_NaN, SubjectName{1});


[sortedBC10, orderBC10] = sort([BC10{:,2}], 'descend');            % rank based on highst value first
BC10 = BC10(orderBC10, :);
[~, ~, ic_BC10] = unique(sortedBC10);
Ranked_BC10 = num2cell(max(ic_BC10)-ic_BC10+1);
BC10(:,3) = Ranked_BC10;
sort_BC10 = sortrows([BC10], 'ascend') ;
rank_sort_BC10 = [sort_BC10{:,3}]';


% For 2 stimuli
BC2(:,1) = ccep10.ch;
BC2(:,2) = num2cell(agreement_parameter.BCN_2') ;


% remove NaN's
remove_NaN_2 = cellfun(@isnan,BC2(:,2),'UniformOutput',false);
tot_NaN_2 = sum([remove_NaN_2{:}]);
for i = 1:length(remove_NaN_2)
    if ismember(remove_NaN_2{i,1}, 1) 
        BC2{i,2} = 0;
    end
end
fprintf('%d NaNs are replaced with 0 for 2 stims, %s \n', tot_NaN_2, SubjectName{1});


[sortedBC2, orderBC2] = sort([BC2{:,2}], 'descend');            % rank based on highst value first
BC2 = BC2(orderBC2, :);
[~, ~, ic_BC2] = unique(sortedBC2);
Ranked_BC2 = num2cell(max(ic_BC2)-ic_BC2+1);
BC2(:,3) = Ranked_BC2;
sort_BC2 = sortrows([BC2], 'ascend') ;
rank_sort_BC2 = [sort_BC2{:,3}]';


% Check
if ~ismember(sort_BC10(:,1),sort_BC2(:,1))
    fprintf('The stimulation pair order of the two situations is not equal. \n');
end

[RHO_BC, PVAL_BC] = corr(rank_sort_BC10 , rank_sort_BC2 ,'Type','Spearman');

fprintf('Spearman Corr between outdegree per electrode of 10 and 2 stimuli gives, p-value = %1.4f, rho = %1.3f, for %s \n', PVAL_BC, RHO_BC, SubjectName{1});

end

