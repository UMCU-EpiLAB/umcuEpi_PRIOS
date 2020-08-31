function [agreement_run] = determine_agreement_checked(myDataPath,cfg)

files = dir(fullfile(myDataPath.CCEPpath, cfg.sub_labels{1}, 'ses-*' ,cfg.run_label{1},...
    [cfg.sub_labels{1} '*_checked*.mat']));

% Use the function below when >2 SPES sessions have to be compared, be
% aware that the values for W, Z and XandY have to be changed.
% files = dir(fullfile(myDataPath.CCEPpath, cfg.sub_labels{1}, 'ses-*' ,'run-*',...
%    [cfg.sub_labels{1} '_ses-*_' cfg.task_label '_*'  '_CCEP_*.mat']));

runs = struct;
if isempty(files)
    fprintf('WARNING: No runs are found');
else
    for i = 1:size(files,1)
        names = fullfile(files(i).folder, files(i).name);
        ccep = load(names);
        runs(i).name = names;
        runs(i).ccep = ccep.ccep;
    end
end

%% same run: 10stims and 2stims checked
% loop through all runs to find 10stims and 2stims of the same run, to
% compare adjacency matrices and calculate agreement
countnum = 1;
agreement_run = struct;
if size(runs,2) >1
    for i=1:size(runs,2)-1
        for j=1:size(runs,2)-1
            
            OA = NaN; PA = NaN; NA = NaN;
            
                if ~strcmp( extractBetween(runs(i).name,'CCEP_','_checked'), extractBetween(runs(i+j).name,'CCEP_','_checked')) % if stim num is not equal
                    if all(size(runs(i).ccep.n1_peak_amplitude) == size(runs(i+j).ccep.n1_peak_amplitude)) % if size of adjacency matrix is equal
                        
                        title1 = extractBetween(runs(i).name,'CCEP_','_checked');      % 10 stims
                        title2 = extractBetween(runs(i+j).name,'CCEP_','_checked');    % 2 stims    
                        run_label = extractBetween(runs(i).name,'_run-','_CCEP');
                        Amat1 = runs(i).ccep.n1_peak_amplitude_check;     % visually checked CCEPs
                        Amat2 = runs(i+j).ccep.n1_peak_amplitude_check;
                        
                        Amat1(~isnan(Amat1)) = 1;                   % all non-NaNs are amplitudes, so N1s --> 1
                        Amat1(isnan(Amat1)) = 0;                    % all NaNs are no N1s --> 0
            
                        Amat2(~isnan(Amat2)) = 1;
                        Amat2(isnan(Amat2)) = 0;
                        
                        compare_mat = Amat1 + Amat2;
                        dif_mat = Amat1-Amat2;                      % 1 = ER in 10stims, -1 = ER in 2stims
                        truetrue = sum(compare_mat(:) == 2);
                        truefalse = sum(compare_mat(:) == 1);
                        falsefalse = sum(compare_mat(:) == 0);
                        
                        total = truetrue + truefalse + falsefalse;
                        
                        % total number of ones in the compare matrix (1 = ER versus non-ER)
                        for i = 1:size(compare_mat,2)
                            TotOnesStim(i,1) = sum(compare_mat(:,i) == 1) ;         % Number of ones detected per stimulation pair
                        end

                        % number of ERs detected with 10 stims
                        for i = 1:size(Amat1,2)
                            Tot10(i,1) = sum(Amat1(:,i) == 1) ;      % ERs detected per stimulation pair    
                        end
                        TotERs10 = sum(Tot10);

                        % number of ERs detected with 2 stims
                        for i = 1:size(Amat2,2)
                            Tot2(i,1) = sum(Amat2(:,i) == 1) ;       % ERs detected per stimulation pair   
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
        agreement_run(countnum).OA = OA;
        agreement_run(countnum).PA = PA;
        agreement_run(countnum).NA = NA;
        agreement_run(countnum).compare_mat = compare_mat;

        
        figure('Position',[680 137 1036 794]),
        subplot(2,3,1), imagesc(Amat1), xlabel('Stimpairs'), ylabel('Channels'), title(sprintf('%s in run %s',title1{:},run_label{:}))
        subplot(2,3,2), imagesc(Amat2),  xlabel('Stimpairs'), ylabel('Channels'),title(sprintf('%s in run %s',title2{:},run_label{:}))
        subplot(2,3,3), imagesc(compare_mat),  xlabel('Stimpairs'), ylabel('Channels'),title(sprintf('Comparison of %s vs %s',title1{:},title2{:}))
        subplot(2,3,4), imagesc((Amat1-Amat2)==1), xlabel('Stimpairs'), ylabel('Channels'), title(sprintf('only in %s',title1{:}))
        subplot(2,3,5), imagesc((Amat2-Amat1)==1), xlabel('Stimpairs'), ylabel('Channels'), title(sprintf('only in %s',title2{:}))
        
        countnum = countnum +1;
    end
end
end
    
