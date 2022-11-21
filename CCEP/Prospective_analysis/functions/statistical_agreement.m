function [statistics, rank] = statistical_agreement(agreement_parameter, ccep)
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
statistics.p_ERsperStimp = p;

clear PVAL_stmp
        
end       

