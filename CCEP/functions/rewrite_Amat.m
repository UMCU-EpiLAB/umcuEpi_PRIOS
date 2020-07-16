function dataBase = rewrite_Amat(dataBase,Amat1)
    stimsets = dataBase.cc_stimsets_avg  ;
    col1 = stimsets(:,1);
    col2 = stimsets(:,2);
    ERs_col = Amat1';                                                      % Electrodes in columns, stimpairs in rows
    
    
    new_mat = zeros(size(ERs_col,2),size(ERs_col,2));
    for i = 1:size(ERs_col,2)
        elecname = dataBase.ch{i};
        if ismember(i,col1,'rows') &&  ismember(i,col2,'rows')              % True --> electrode is part of two stimpairs   
           [~,loc1] =  ismember(i,col1,'rows') ;                            % Find row with the electrode of the first stimair 
           [~,loc2] = ismember(i,col2,'rows') ;                             % Find row with the electrode of the second stimair 
           new_mat(i,:) = ERs_col(loc1,:) + ERs_col(loc2,:);                  % Number of ERs detected per electrode present in two stimpairs
        elseif ismember(i,col1,'rows')                                      % electrode is only in one stimpair
            [~,loc1] =  ismember(i,col1,'rows') ; 
            new_mat(i,:) = ERs_col(loc1,:); 
        elseif ismember(i,col2,'rows')
            [~,loc2] =  ismember(i,col2,'rows') ; 
            new_mat(i,:) = ERs_col(loc2,:); 
        end
    end
    
    ColN = dataBase.ch';
    rowNames = dataBase.ch;
    indegree = array2table(new_mat,'RowNames',rowNames,'VariableNames',ColN);   %  
    dataBase.elec_Amat = new_mat;
end
