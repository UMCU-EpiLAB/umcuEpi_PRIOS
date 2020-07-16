function [indegree, outdegree, BC, rank_stimp, rank_elec] = agreement_parameters(Amat2, dataBase, stimchans)

    wantedAmat = Amat2';                                    % Stimparen are the source nodes (rows) elektrodes are the target nodes (columns)                 
    Amat = zeros(size(wantedAmat,1)+size(wantedAmat,2));    % Matrix with 4 quadrants, left under must be filled with the adjecency matrix            
    rowStart = size(wantedAmat,2)+1;
    columnEnd = size(wantedAmat,2);
    ElecName = dataBase.ch;
    StimName = stimchans;
       
    Amat((rowStart:end),(1:columnEnd)) = wantedAmat;
    
    G = digraph(Amat);                         % Matrix must have same length and width
    % plot(G,'NodeLabel',{'alpha','beta','gamma','delta','epsilon','zeta'})

        % NETWERKMATEN
        % BC IS ALWAYS ZERO BECAUSE ALL THE WEIGHTS ARE THE SAME = 1
        % THEREFORE NO ELECTRODE OR STIMPAIR IS MORE IMPORTANT. 
        
        BC = centrality(G,'betweenness','Cost',G.Edges.Weight);    % Weighted betweenness centrality
       
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
      
%% Work with the electrodes in the columns and rows 
    elec_mat = dataBase.elec_Amat;
    G = digraph(elec_mat);
    p = plot(G,'Layout','force','NodeLabel',dataBase.ch,'NodeColor','r','MarkerSize',bins);
    %%% Misschien mag iets als linesize of linecolour ook wel gebaseerd op
    %%% indegree en outdegree.

    wbc = centrality(G,'betweenness','Cost',G.Edges.Weight);        % 'Cost' edge weights specify the length of the edges.
    hub_ranks = centrality(G,'hubs');
    Indegree = centrality(G,'indegree','Importance',G.Edges.Weight);
    Outdegree = centrality(G,'outdegree','Importance',G.Edges.Weight);
    
    edges = linspace(min(Indegree),max(Indegree),7);            % create 7 equally-spaced bins based on their indegree score
    bins = discretize(Indegree,edges);
     
    
    
    
%     G.Nodes.Hubs = hub_ranks;
%     n = numnodes(G);
%     p.NodeCData = wbc./((n-2)*(n-1));
%     colormap(flip(autumn,1));
%     title('Betweenness Centrality Scores - Weighted')
    
    
     % total number of electrodes
    stimelektot = size(dataBase.ch,1);
    % total number of stimpairs
    stimptot = size(dataBase.stimpnames_avg,2);
    
    % Number of times an electrode is used in a stimulation pair  
    for el=1:size(elec_mat,1)
        trialelek(el) = size(find(dataBase.cc_stimsets_avg==el),1);
    end
    
    % totaal number of possible connections
    for el=1:size(elec_mat,1)
        n_outtot(el) = trialelek(el)*(stimelektot-2);
        n_intot(el) = 2*(stimptot - trialelek(el)); 
    end
    
    % NETWERKMATEN         
    for el=1:size(elec_mat,1)
        % indegree en outdegree = som van 1tjes
        indegreenorm(el) = indegree(el)/n_intot(el);
        outdegreenorm(el) = outdegree(el)/n_outtot(el);
        BCnorm(el) = wbc(el)/(n_intot(el)*n_outtot(el));
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
    