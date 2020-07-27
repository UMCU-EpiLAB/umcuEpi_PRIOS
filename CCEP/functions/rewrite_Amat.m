function dataBase = rewrite_Amat(dataBase,Amat1)
stimsets = dataBase.cc_stimsets_avg  ;
stim1 = stimsets(:,1);
stim2 = stimsets(:,2);
ERs_col = Amat1';                                                      % Electrodes in columns, stimpairs in rows

elec_mat = zeros(size(ERs_col,2),size(ERs_col,2)); % adjacency matrix with electrodes to electrodes
for chan = 1:size(ERs_col,2)
    %         elecname = dataBase.ch{i};
    if ismember(chan,stim1,'rows') &&  ismember(chan,stim2,'rows')              % True --> electrode is part of two stimpairs
        
        [~,loc1] =  ismember(chan,stim1,'rows') ;                            % Find row with the electrode of the first stimair
        [~,loc2] = ismember(chan,stim2,'rows') ;                             % Find row with the electrode of the second stimair
        elec_mat(chan,:) = ERs_col(loc1,:) + ERs_col(loc2,:);                  % Number of ERs detected per electrode present in two stimpairs
        
    elseif ismember(chan,stim1,'rows')                                      % electrode is only in one stimpair
        
        [~,loc1] =  ismember(chan,stim1,'rows') ;
        elec_mat(chan,:) = ERs_col(loc1,:);
        
    elseif ismember(chan,stim2,'rows')
        
        [~,loc2] =  ismember(chan,stim2,'rows') ;
        elec_mat(chan,:) = ERs_col(loc2,:);
    end
end

%     ColN = dataBase.ch';
%     rowNames = dataBase.ch;
%     indegree = array2table(new_mat,'RowNames',rowNames,'VariableNames',ColN);   %
dataBase.elec_Amat = elec_mat;
end
