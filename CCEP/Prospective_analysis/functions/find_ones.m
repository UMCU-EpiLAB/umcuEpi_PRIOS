function LocOnes = find_ones(dataBase,agreement_run)

    stimchans =  dataBase_clin.stimpnames_avg';

    % [electrodes,stimpairs]. 
    [x(:,2),x(:,1)] = find(compare_mat ==1);
    %x_ColN = {'Stimpairs','Electrodes'};
    %x = array2table(x,'VariableNames',x_ColN);

    FindOnes = compare_mat;
    rowNames = dataBase.ch;
    colNames = stimchans;
    FindOnes = array2table(FindOnes,'RowNames',rowNames,'VariableNames',colNames);

    % de stimulatie paren
    for i = 1:size(x(:,2))                            % For the number of ones detected
        LocOnes(i,1) = stimchans(x(i,1))';
        LocOnes(i,2) = dataBase(1).ch(x(i,2));
        
    end

    ColN = {'Stimpairs','Electrodes'};
    LocOnes = array2table(LocOnes,'VariableNames',ColN);
     
end
