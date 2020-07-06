function [indegree, outdegree, BC] = agreement_parameters(Amat2)

    wantedAmat = Amat2';                                    % Stimparen are the source nodes (rows) elektrodes are the target nodes (columns)                 
    Amat = zeros(size(wantedAmat,1)+size(wantedAmat,2)); % Matrix with 4 quadrants, left under must be filled with the adjecency matrix            
    rowStart = size(wantedAmat,2)+1;
    columnEnd = size(wantedAmat,2);
    
    Amat((rowStart:end),(1:columnEnd)) = wantedAmat;
    
    G = digraph(Amat);                         % Matrix must have same length and width
    plot(G,'Layout','force')

        % NETWERKMATEN
        BC = centrality(G,'betweenness','Cost',G.Edges.Weight);

        for stimPair = size(wantedAmat,2)+1:size(Amat,1)
            outdegree(stimPair-60,1) = sum(Amat(stimPair,:));                % total number of ERs evoked by a stimulation pair 
            %T = minspantree(minnGraph)                                %minimum graph thats connects all of the nodes.
            % [tree,d] =shortstpathtree(minnGraph,1561)
        end
        
        for elec = 1:size(wantedAmat,2)
            indegree(1,elec) = sum(Amat(:,elec));                            % total number of ERs evoked in an electrode (el) 
        end
    
end


     
    
% % stimulation pairs/sets
%     sort_stimsets = sort(ccep.ccep.stimsets,2);                         % all stimsets negative and positive toghether
%     stimsets_dir = unique(sort_stimsets,'rows');                        % every unique stimset
%     stimchan = [ccep.ccep.ch(stimsets_dir)]  ;                          % find all names matching all unique stimsets  
%     trialtot = size(stimchan,1);                                        % total number of trials (verschillende mogelijkheden tot stimuleren)
% % Electrodes  
%     stimelek = ccep.ccep.ch;                                            % all stimulated electrodes 
%     stimelektot = size(stimelek,1);                                     % total number of stimulated electrodes         
%        
     
% % Number of times an elektrode is stimulated
% % Find how often a elektrode number of found in the stimsets_dir
%         for el=1:size(stimelek,1)
%             trialelek(el,1) = size(find(stimsets_dir==el),1);
%         end
% 
% % Total number of possible connections 
%     for el=1:size(stimelek,1)                                               % number of electrodes
%         % OUT
%         % the number of trials in which the specific electrode is
%         % stimulated (trialelek(el)) multiplied by the total number of
%         % possible responses electrodes (stimelektot) minus 2. 
%         n_outtot(el) = trialelek(el)*(stimelektot-2);
% 
%         % IN
%         % 2 (the number of electrodes in a stimulation pair) multiplied by
%         % the total number of trials (trialtot) minus the number of trials
%         % in which the specific electrode was stimulated (trialelek)
%         n_intot(el) = 2*(trialtot - trialelek(el));
%     end
    