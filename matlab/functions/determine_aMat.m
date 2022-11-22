function [Amat] = determine_aMat(runs)
% Use the function below when 2 SPES sessions have to be compared.
if size(runs,2) == 2    
    if all(size(runs(1).ccep.n1_peak_sample) == size(runs(2).ccep.n1_peak_sample))   % Extra check whether size of adjacency matrix is equal

        AmatClin = ~isnan(runs(1).ccep.n1_peak_sample); % all non-NaNs are samples, so N1s --> 1
        AmatProp = ~isnan(runs(2).ccep.n1_peak_sample);
        
        % Save value to struct to use later
        Amat.AmatProp = AmatProp;
        Amat.AmatClin = AmatClin;
    else
        error('Size of adjacency matrices is not equal!')
    end

else
    warning('More that two runs are loaded')
end

end
