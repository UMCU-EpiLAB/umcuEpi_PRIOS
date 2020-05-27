function [agreement_run, agreement_stim, compare_mat] = determine_agreement(myDataPath,cfg)

files = dir(fullfile(myDataPath.CCEPpath, cfg.sub_labels{1}, 'ses-*' ,cfg.run_label{1},...
    [cfg.sub_labels{1} '_ses-*_' cfg.task_label '_*'  '_CCEP_*.mat']));

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

%% same run: 10stims and 2stims
% loop through all runs to find 10stims and 2stims of the same run, to
% compare adjacency matrices and calculate agreement
countnum = 1;
agreement_run = struct;
if size(runs,2) >1
    for i=1:size(runs,2)-1
        for j=1:size(runs,2)-1
            
            OA = NaN; PA = NaN; NA = NaN;
            
            if strcmp(extractBetween(runs(i).name,'_run-','_CCEP'), extractBetween(runs(i+j).name,'_run-','_CCEP')) % if run label is equal
                if ~strcmp( extractBetween(runs(i).name,'_CCEP_','.mat'), extractBetween(runs(i+j).name,'_CCEP_','.mat')) % if stim num is not equal
                    if all(size(runs(i).ccep.n1_peak_amplitude) == size(runs(i+j).ccep.n1_peak_amplitude)) % if size of adjacency matrix is equal
                        
                        title1 = extractBetween(runs(i).name,'_CCEP_','.mat');
                        title2 = extractBetween(runs(i+j).name,'_CCEP_','.mat');
                        run_label = extractBetween(runs(i).name,'_run-','_CCEP');
                        Amat1 = runs(i).ccep.n1_peak_amplitude;
                        Amat2 = runs(i+j).ccep.n1_peak_amplitude;
                        
                        Amat1(~isnan(Amat1)) = 1; % all non-NaNs are amplitudes, so N1s --> 1
                        Amat1(isnan(Amat1)) = 0; % all NaNs are no N1s --> 0
                        
                        Amat2(~isnan(Amat2)) = 1;
                        Amat2(isnan(Amat2)) = 0;
                        
                        compare_mat = Amat1 + Amat2;
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
    
%% two different runs
% loop through all runs to find 10stims and 10stims of two different runs, to
% compare adjacency matrices and calculate agreement
countnum = 1;
agreement_stim = struct;
if size(runs,2) >1
    for i=1:size(runs,2)-1
        for j=1:size(runs,2)-1
            
            OA = NaN; PA = NaN; NA = NaN;
            
            if strcmp(extractBetween(runs(i).name,'_CCEP_','.mat'), extractBetween(runs(i+j).name,'_CCEP_','.mat')) % if stim size is equal
                if ~strcmp(extractBetween(runs(i).name,'_run-','_CCEP'),extractBetween(runs(i+j).name,'_run-','_CCEP')) % if run label is not equal
                if all(size(runs(i).ccep.n1_peak_amplitude) == size(runs(i+j).ccep.n1_peak_amplitude)) % if size of adjacency matrix is equal
                    
                    title1 = extractBetween(runs(i).name,'_run-','_CCEP');
                    title2 = extractBetween(runs(i+j).name,'_run-','_CCEP');
                    run_label = extractBetween(runs(i).name,'_CCEP_','.mat');
                    Amat1 = runs(i).ccep.n1_peak_amplitude;
                    Amat2 = runs(i+j).ccep.n1_peak_amplitude;
                    
                    Amat1(~isnan(Amat1)) = 1; % all non-NaNs are amplitudes, so N1s --> 1
                    Amat1(isnan(Amat1)) = 0; % all NaNs are no N1s --> 0
                    
                    Amat2(~isnan(Amat2)) = 1;
                    Amat2(isnan(Amat2)) = 0;
                    
                    compare_mat = Amat1 + Amat2;
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
                
            end
            end
        end
        agreement_stim(countnum).OA = OA;
        agreement_stim(countnum).PA = PA;
        agreement_stim(countnum).NA = NA;
        agreement_stim(countnum).compare_mat = compare_mat;
        
        figure('Position',[680 137 1036 794]), 
        subplot(2,3,1), imagesc(Amat1), xlabel('Stimpairs'), ylabel('Channels'), title(sprintf('%s in run %s',title1{:},run_label{:}))
        subplot(2,3,2), imagesc(Amat2),  xlabel('Stimpairs'), ylabel('Channels'),title(sprintf('%s in run %s',title2{:},run_label{:}))
        subplot(2,3,3), imagesc(compare_mat),  xlabel('Stimpairs'), ylabel('Channels'),title(sprintf('Comparison of %s vs %s',title1{:},title2{:}))
        subplot(2,3,4), imagesc((Amat1-Amat2)==1), xlabel('Stimpairs'), ylabel('Channels'), title(sprintf('only in %s',title1{:}))
        subplot(2,3,5), imagesc((Amat2-Amat1)==1), xlabel('Stimpairs'), ylabel('Channels'), title(sprintf('only in %s',title2{:}))
        
        countnum = countnum +1;
    end
end

    
% 
% %%
% 
% % Compare the electrodes which show a N1 in different networks.
% if length(runs) > 1
%     % for the number of files to compare
%     for n = 1:(length(runs)-1)
%         % Find the extra stimulation pairs
%         C = setdiff(runs(n+1).ccep.stimpnames,runs(n).ccep.stimpnames);
%         same_stimpair1 = setdiff(runs(n+1).ccep.stimpnames, C);
%         
%         D = setdiff(runs(n).ccep.stimpnames, runs(n+1).ccep.stimpnames);
%         same_stimpair2 = setdiff(runs(n).ccep.stimpnames, D);
%         
%         if length(same_stimpair1) ~= length(same_stimpair2)
%             fprintf('WARNING: Extra stimpairs of %s are NOT removed.\n', cfg.sub_labels{:});
%         end
%         
%         % Find the extra electrodes
%         E = setdiff(runs(n+1).ccep.ch', runs(n).ccep.ch');
%         same_chan1 = setdiff(runs(n+1).ccep.ch', E);
%         
%         F = setdiff(runs(n).ccep.ch', runs(n+1).ccep.ch');
%         same_chan2 = setdiff(runs(n).ccep.ch', F);
%         
%         if length(same_chan1) ~= length(same_chan2)
%             fprintf('WARNING: Extra channels of %s are NOT removed.\n', cfg.sub_labels{:});
%         end
%         
%         for Q = 1:length(runs)
%             correct_stimpairs{Q,1} = find(ismember(runs(Q).ccep.stimpnames, same_stimpair2));
%             correct_chans{Q,1} = find(ismember(runs(Q).ccep.ch', same_chan2));
%             % does not matter whether same_chan1 or 2, should be the same
%             
%             % Make matrix with the channels and stimpairs which are
%             % present in all runs
%             matrix{Q,1} = runs(Q).ccep.n1_peak_sample(correct_chans{Q}, correct_stimpairs{Q});
%             adjacency_matrix{Q,1} = isnan(matrix{Q,1});  % NAN = 1, value = 0
%         end
%         % For two runs, they can be added and the values will be 2, 1
%         % or 0. When multiple runs are compared, new values have to be
%         % filled in.
%         if length(runs) == 2
%             compare_mat = adjacency_matrix{1} + adjacency_matrix{2};
%             Z = sum(compare_mat(:)==2);         % Both no-ER
%             XandY = sum(compare_mat(:)==1);     % One NAN one ER
%             W = sum(compare_mat(:)==0);         % Two ERs
%             total = W + XandY + Z ;
%             check = size(compare_mat,1) .* size(compare_mat,2);
%         else
%             fprintf('WARNING: The matrices of >2 runs cannot be added up');
%         end
%     end
%     %end
%     
%     % Overall, positive and negative agreement between the matrices.
%     OA = (W + Z) / (W + XandY + Z);
%     PA = (2 * W) / (2 * W + XandY);
%     NA = (2 * Z) / (2 * Z + XandY);
%     
% end