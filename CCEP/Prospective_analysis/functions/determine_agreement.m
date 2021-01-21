function [agreement] = determine_agreement(runs)

% Use the function below when >2 SPES sessions have to be compared, be
% aware that the values for W, Z and XandY have to be changed.
% files = dir(fullfile(myDataPath.CCEPpath, cfg.sub_labels{1}, 'ses-*' ,'run-*',...
%    [cfg.sub_labels{1} '_ses-*_' cfg.task_label '_*'  '_CCEP_*.mat']));

%% two different runs
% loop through all runs to find 10stims and 10stims of two different runs, to
% compare adjacency matrices and calculate agreement
countnum = 1;
agreement_stim = struct;
if size(runs,2) >1
    for i=1:size(runs,2)-1
        for j=1:size(runs,2)-1
            
            OA = NaN; PA = NaN; NA = NaN;
            
            if ~strcmp(extractBetween(runs(i).name,'_CCEP_','.mat'), extractBetween(runs(i+j).name,'_CCEP_','.mat'))        % if stim size is NOT equal
                if ~strcmp(extractBetween(runs(i).name,'_run-','_CCEP'), extractBetween(runs(i+j).name,'_run-','_CCEP'))    % if run label is not equal
                    if all(size(runs(i).ccep.n1_peak_amplitude_check) == size(runs(i+j).ccep.n1_peak_amplitude_check))                  % if size of adjacency matrix is equal
                        
                        titleclin = extractBetween(runs(i).name,'_run-','_CCEP');
                        titleprop = extractBetween(runs(i+j).name,'_run-','_CCEP');
                        run_labelClin = extractBetween(runs(i).name,'_CCEP_','.mat');
                        run_labelProp = extractBetween(runs(i+j).name,'_CCEP_','.mat');

                        AmatClin = runs(i).ccep.n1_peak_amplitude_check;
                        AmatProp = runs(i+j).ccep.n1_peak_amplitude_check;
                        
                        AmatClin(~isnan(AmatClin)) = 1; % all non-NaNs are amplitudes, so N1s --> 1
                        AmatClin(isnan(AmatClin)) = 0; % all NaNs are no N1s --> 0
                        
                        AmatProp(~isnan(AmatProp)) = 1;
                        AmatProp(isnan(AmatProp)) = 0;
                        
                        compare_mat = AmatClin + AmatProp;
                        truetrue = sum(compare_mat(:) == 2);
                        truefalse = sum(compare_mat(:) == 1);
                        falsefalse = sum(compare_mat(:) == 0);
                        
                        total = truetrue + truefalse + falsefalse;
                        
                        if total ~= size(compare_mat,1) * size(compare_mat,2)
                            error('Summation of all values is not equal to size of adjacency matrix')
                            
                        else
                            
                            % Overall, positive and negative agreement between the matrices.
                            OA = (truetrue + falsefalse) / (truetrue+ truefalse + falsefalse);
                            PA = (2 * truetrue) / ((2 * truetrue) + truefalse);
                            NA = (2 * falsefalse) / ((2 * falsefalse) + truefalse);
                            
                        end
                        
                    else
                        error('Size of adjacency matrices is not equal!')
                    end
                    
                    figure('Position',[680 137 1036 794]),
                    subplot(2,3,1), imagesc(AmatClin), xlabel('Stimpairs'), ylabel('Channels'), title(sprintf('%s in run %s',titleclin{:},run_labelClin{:}))
                    subplot(2,3,2), imagesc(AmatProp),  xlabel('Stimpairs'), ylabel('Channels'),title(sprintf('%s in run %s',titleprop{:},run_labelProp{:}))
                    subplot(2,3,3), imagesc(compare_mat),  xlabel('Stimpairs'), ylabel('Channels'),title(sprintf('Comparison of %s vs %s',titleclin{:},titleprop{:}))
                    subplot(2,3,4), imagesc((AmatClin-AmatProp)==1), xlabel('Stimpairs'), ylabel('Channels'), title(sprintf('only in %s',titleclin{:}))
                    subplot(2,3,5), imagesc((AmatProp-AmatClin)==1), xlabel('Stimpairs'), ylabel('Channels'), title(sprintf('only in %s',titleprop{:}))
                    
                    
                end
            end
        end
        agreement_stim(countnum).OA = OA;
        agreement_stim(countnum).PA = PA;
        agreement_stim(countnum).NA = NA;
        agreement_stim(countnum).compare_mat = compare_mat;
        
        countnum = countnum +1;
    end
end


agreement.agreement_stim = agreement_stim;
agreement.compare_mat = compare_mat;
agreement.AmatClin = AmatClin;
agreement.AmatProp = AmatProp;
end
