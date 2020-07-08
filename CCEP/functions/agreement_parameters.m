function [indegree, outdegree, BC, rank_stimp, rank_elec] = agreement_parameters(Amat2, dataBase, stimchans)

    wantedAmat = Amat2';                                    % Stimparen are the source nodes (rows) elektrodes are the target nodes (columns)                 
    Amat = zeros(size(wantedAmat,1)+size(wantedAmat,2));    % Matrix with 4 quadrants, left under must be filled with the adjecency matrix            
    rowStart = size(wantedAmat,2)+1;
    columnEnd = size(wantedAmat,2);
    ElecName = dataBase.ch;
    StimName = stimchans;
       
    Amat((rowStart:end),(1:columnEnd)) = wantedAmat;
    
    G = digraph(Amat);                         % Matrix must have same length and width
    plot(G,'Layout','force')
    % plot(G,'NodeLabel',{'alpha','beta','gamma','delta','epsilon','zeta'})

        % NETWERKMATEN
        % BC IS ALWAYS ZERO BECAUSE ALL THE WEIGHTS ARE THE SAME = 1
        % THEREFORE NO ELECTRODE OR STIMPAIR IS MORE IMPORTANT. 
        
        BC = centrality(G,'betweenness','Cost',G.Edges.Weight);    % Weighted betweenness centrality
        n = numnodes(G);                                        
        p.NodeCData = 2*BC./((n-2)*(n-1));                         % normalised WBC
        colormap(flip(autumn,1));                                   % Plot
        title('Betweenness Centrality Scores - Weighted')


        for stimPair = size(wantedAmat,2)+1:size(Amat,1)
            outdegree((stimPair-(rowStart-1)),1) = sum(Amat(stimPair,:));                % total number of ERs evoked by a stimulation pair 
            %T = minspantree(minnGraph)                                %minimum graph thats connects all of the nodes.
            % [tree,d] =shortstpathtree(minnGraph,1561)
        end
        
        ColN = {'number of ERs'};
        rowNames = StimName;
        rank_stimp = array2table(outdegree,'RowNames',rowNames,'VariableNames',ColN);   % Outdegree
        rank_stimp = sortrows(rank_stimp,1,'descend');
        
        for elec = 1:size(wantedAmat,2)
            indegree(1,elec) = sum(Amat(:,elec));                            % total number of ERs evoked in an electrode (el) 
        end
        
        indegree = indegree';
        RowNames = ElecName;
        rank_elec = array2table(indegree,'RowNames',RowNames,'VariableNames',ColN);   % Indegree
        rank_elec = sortrows(rank_elec,1,'descend');
      
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
    