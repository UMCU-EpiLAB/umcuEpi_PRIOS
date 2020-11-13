function [agreement] = determine_agreement(runs)

% Same run: 10stims and 2stims
% loop through all runs to find 10stims and 2stims of the same run, to
% compare adjacency matrices and calculate agreement

countnum = 1;
agreement_run = struct;
for i=1:size(runs,2)-1
    for j=1:size(runs,2)-1
        
        OA = NaN; PA = NaN; NA = NaN;                   % overall agreement, positive agreement, negative agreement
        
        if strcmp(extractBetween(runs(i).name,'_run-','_CCEP'), extractBetween(runs(i+j).name,'_run-','_CCEP')) % if run label is equal
            if ~strcmp( extractBetween(runs(i).name,'_CCEP_','.mat'), extractBetween(runs(i+j).name,'_CCEP_','.mat')) % if stim num is not equal
                if all(size(runs(i).ccep.n1_peak_amplitude) == size(runs(i+j).ccep.n1_peak_amplitude)) % if size of adjacency matrix is equal
                    
                    title10 = extractBetween(runs(i).name,'_CCEP_','.mat');     % 10 stims
                    title2 = extractBetween(runs(i+j).name,'_CCEP_','.mat');    % 2 stims
                    run_label = extractBetween(runs(i).name,'_run-','_CCEP');
                    Amat10 = runs(i).ccep.n1_peak_amplitude;
                    Amat2 = runs(i+j).ccep.n1_peak_amplitude;
                    
                    Amat10(~isnan(Amat10)) = 1;                   % all non-NaNs are amplitudes, so N1s --> 1
                    Amat10(isnan(Amat10)) = 0;                    % all NaNs are no N1s --> 0
                    
                    Amat2(~isnan(Amat2)) = 1;
                    Amat2(isnan(Amat2)) = 0;
                    
                    compare_mat = Amat10 + Amat2;
                    dif_mat = Amat10-Amat2;                      % 1 = ER in 10stims, -1 = ER in 2stims
                    truetrue = sum(compare_mat(:) == 2);         % both are ER
                    truefalse = sum(compare_mat(:) == 1);        % one is ER and other is non-ER
                    falsefalse = sum(compare_mat(:) == 0);       % both are non-ER
                    
                    total = truetrue + truefalse + falsefalse;
                    
                    % total number of ones in the compare matrix (1 = ER versus non-ER)
                    TotOnesStim = zeros(size(compare_mat,2),1);
                    for n = 1:size(compare_mat,2)
                        TotOnesStim(n,1) = sum(compare_mat(:,n) == 1) ;         % Number of ones detected per stimulation pair
                    end
                    
                    % number of ERs detected with 10 stims
                    Tot10 = zeros(size(Amat10,2),1);
                    for n = 1:size(Amat10,2)
                        Tot10(n,1) = sum(Amat10(:,n) == 1) ;      % ERs detected per stimulation pair
                    end
                    TotERs10 = sum(Tot10);
                    
                    % number of ERs detected with 2 stims
                    Tot2 = zeros(size(Amat2,2),1);
                    for n = 1:size(Amat2,2)
                        Tot2(n,1) = sum(Amat2(:,n) == 1) ;       % ERs detected per stimulation pair
                    end
                    TotERs2 = sum(Tot2);
                    
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
                
            end
        end
    end
    
    agreement_run(countnum).OA = OA;
    agreement_run(countnum).PA = PA;
    agreement_run(countnum).NA = NA;
    agreement_run(countnum).compare_mat = compare_mat;
    
    figure('Name',runs(1).sub_label,'Position',[680 137 1036 794]),
    subplot(2,3,1), imagesc(Amat10), xlabel('Stimpairs'), ylabel('Channels'), title(sprintf('%s in run %s',title10{:},run_label{:}))
    subplot(2,3,2), imagesc(Amat2),  xlabel('Stimpairs'), ylabel('Channels'),title(sprintf('%s in run %s',title2{:},run_label{:}))
    subplot(2,3,3), imagesc(compare_mat),  xlabel('Stimpairs'), ylabel('Channels'),title(sprintf('Comparison of %s vs %s',title10{:},title2{:}))
    subplot(2,3,4), imagesc((Amat10-Amat2)==1), xlabel('Stimpairs'), ylabel('Channels'), title(sprintf('only in %s',title10{:}))
    subplot(2,3,5), imagesc((Amat2-Amat10)==1), xlabel('Stimpairs'), ylabel('Channels'), title(sprintf('only in %s',title2{:}))
    
    countnum = countnum +1;
end

agreement.agreement_run =agreement_run ;
agreement.compare_mat = compare_mat;
agreement.dif_mat = dif_mat;
agreement.TotERs10 = TotERs10;
agreement.TotERs2 = TotERs2;
agreement.TotOnesStim = TotOnesStim;
agreement.Amat10 = Amat10;
agreement.Amat2 = Amat2;
end