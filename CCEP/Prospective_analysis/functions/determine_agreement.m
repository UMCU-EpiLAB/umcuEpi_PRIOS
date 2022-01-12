function [agreement] = determine_agreement(runs, myDataPath)
% Use the function below when 2 SPES sessions have to be compared.
if size(runs,2) == 2
           
    % Preallocation
    OA = NaN; PA = NaN; NA = NaN;
    
    if all(size(runs(1).ccep.n1_peak_amplitude) == size(runs(2).ccep.n1_peak_amplitude))   % Extra check whether size of adjacency matrix is equal

        AmatClin = ~isnan(runs(1).ccep.n1_peak_amplitude); % all non-NaNs are amplitudes, so N1s --> 1
        AmatProp = ~isnan(runs(2).ccep.n1_peak_amplitude);
        
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
        
        % Visualise the results of the agreement in scaled colours. 
        figure('Position',[680 137 1036 794]),
        subplot(2,3,1), imagesc(AmatClin), xlabel('Stimpairs'), ylabel('Channels'), title(sprintf('%s during %s',runs(1).ccep.sub_label, runs(1).ccep.task_label))
        subplot(2,3,2), imagesc(AmatProp),  xlabel('Stimpairs'), ylabel('Channels'),title(sprintf('%s during %s',runs(2).ccep.sub_label, runs(2).ccep.task_label))
        subplot(2,3,3), imagesc(compare_mat),  xlabel('Stimpairs'), ylabel('Channels'),title(sprintf('Comparison of %s vs %s',runs(1).ccep.task_label,runs(2).ccep.task_label))
        subplot(2,3,4), imagesc((AmatClin-AmatProp)==1), xlabel('Stimpairs'), ylabel('Channels'), title(sprintf('only in %s',runs(1).ccep.task_label))
        subplot(2,3,5), imagesc((AmatProp-AmatClin)==1), xlabel('Stimpairs'), ylabel('Channels'), title(sprintf('only in %s',runs(2).ccep.task_label))
            
        % Save the overall, positive and negative agreement in a struct
        agreement.OA = OA;
        agreement.PA = PA;
        agreement.NA = NA;
        agreement.compare_mat = compare_mat;
        agreement.AmatProp = AmatProp;
        agreement.AmatClin = AmatClin;
    else
        error('Size of adjacency matrices is not equal!')
    end

    % Save agreement
   
    filename = [runs(1).ccep.sub_label,'_agreementPar.mat'];
    filename_fig = [runs(1).ccep.sub_label,'_agreementFig.png'];

    filefolder = fullfile(myDataPath.CCEPpath, 'Visualise_agreement','/','Agreement','/');
    if ~exist(filefolder,'dir')
        mkdir(filefolder)
    end

    % Save file during scoring in case of error
    save(fullfile(filefolder,filename),'-struct','agreement');
    saveas(gcf, fullfile(filefolder,filename_fig));

    close all;


else
    warning('More that two runs are loaded')
end

end
