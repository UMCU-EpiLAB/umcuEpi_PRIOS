function [indegreenorm, outdegreenorm, BCnorm] = agreement_parameters(Amat2,ccep)

    wanted = zeros(length(Amat2));
    wanted(1:numel(Amat2)) = Amat2;
    Amat2 = wanted;

    G = digraph(Amat2);
    
% stimulation pairs/sets
    sort_stimsets = sort(ccep.ccep.stimsets,2);                         % all stimsets negative and positive toghether
    stimsets_dir = unique(sort_stimsets,'rows');                        % every unique stimset
    stimchan = [ccep.ccep.ch(stimsets_dir)]  ;                          % find all names matching all unique stimsets  
    trialtot = size(stimchan,1);                                        % total number of trials (verschillende mogelijkheden tot stimuleren)
% Electrodes  
    stimelek = ccep.ccep.ch;                                            % all stimulated electrodes 
    stimelektot = size(stimelek,1);                                     % total number of stimulated electrodes         
       
     
% Number of times an elektrode is stimulated
% Find how often a elektrode number of found in the stimsets_dir
        for el=1:size(stimelek,1)
            trialelek(el,1) = size(find(stimsets_dir==el),1);
        end

% Total number of possible connections 
    for el=1:size(stimelek,1)                                               % number of electrodes
        % OUT
        % the number of trials in which the specific electrode is
        % stimulated (trialelek(el)) multiplied by the total number of
        % possible responses electrodes (stimelektot) minus 2. 
        n_outtot(el) = trialelek(el)*(stimelektot-2);

        % IN
        % 2 (the number of electrodes in a stimulation pair) multiplied by
        % the total number of trials (trialtot) minus the number of trials
        % in which the specific electrode was stimulated (trialelek)
        n_intot(el) = 2*(trialtot - trialelek(el));

    end

        % NETWERKMATEN
        BC = centrality(G,'betweenness','Cost',G.Edges.Weight);
        %minnGraph = graph(Problem.A);

        for el=1:size(Amat2,1)
            indegree(el) = sum(Amat2(:,el));                            % total number of ERs evoked in electrode (el) per stimulation pair
            indegreenorm(el) = indegree(el)/n_intot(el);
            outdegree(el) = sum(Amat2(el,:));                           % total number of ERs evoked a stimulation pair 
            outdegreenorm(el) = outdegree(el)/n_outtot(el);
            BCnorm(el) = BC(el)/(n_intot(el).*n_outtot(el));
            %T = minspantree(minnGraph)                     %minimum graph thats connects all of the nodes.
            % [tree,d] =shortstpathtree(minnGraph,1561)
        end
    
end
    