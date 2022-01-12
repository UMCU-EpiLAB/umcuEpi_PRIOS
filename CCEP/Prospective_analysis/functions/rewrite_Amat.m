function dataBase = rewrite_Amat(dataBase,Amat)
% Rewrite the adjacency matrices with stimluation pairs in the columns and the electrodes in the rows
% to an adjacency matrix with electrodes in the columns and rows. 
% This way the network characteristics can be determined per electrode
% instead of per stimulation pair.

% Find electrodes in stimulation pair
stimsets = dataBase.stimsets_avg  ;
stim1 = stimsets(:,1);                              % First electrode of a stimulation pair
stim2 = stimsets(:,2);                              % Second electrode of a stimulation pair

% Inverse original adjacency matrix to electrodes in columns, stimpairs in rows
ERs_col = Amat';                                    

% Pre-allocation
% Elec_mat is used to determine the number of responses per electrode. 
elec_mat = zeros(size(ERs_col,2),size(ERs_col,2));                    % Adjacency matrix with electrodes to electrodes

% First determine how many stimulation pairs an electrode is part of
for chan = 1:size(ERs_col,2)

    if ismember(chan,stim1,'rows') &&  ismember(chan,stim2,'rows')           % True --> electrode is part of two stimpairs
        
        [~,loc1] =  ismember(chan,stim1,'rows') ;                            % Find row with the electrode of the first stimair
        [~,loc2] = ismember(chan,stim2,'rows') ;                             % Find row with the electrode of the second stimair
        elec_mat(chan,:) = ERs_col(loc1,:) + ERs_col(loc2,:);                % Number of ERs detected per electrode present in two stimpairs
        
    elseif ismember(chan,stim1,'rows')                                       % electrode is only in one stimpair
        
        [~,loc1] =  ismember(chan,stim1,'rows') ;
        elec_mat(chan,:) = ERs_col(loc1,:);                                  % Fill the elec_mat with an 1 when a ER is detected on this electrode which is in only 1 stimpair
        
    elseif ismember(chan,stim2,'rows')
        
        [~,loc2] =  ismember(chan,stim2,'rows') ;
        elec_mat(chan,:) = ERs_col(loc2,:);
    end
end
dataBase.elec_Amat = elec_mat;
end
