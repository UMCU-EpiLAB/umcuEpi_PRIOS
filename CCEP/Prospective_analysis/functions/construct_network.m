function [AmatERsvis2, stimelektot, trialtot, trialelek, n_intot, n_outtot,...
    indegreenorm, outdegreenorm, BCnorm] = ...
    construct_network(AmatERsvis,chanoverlap,stimchan,stimelek)

% NETWERK OPSTELLEN EN PARAMETERS ERUIT HALEN
     AmatERsvis2 = AmatERsvis;
%     % zorgen dat alle -1tjes (stimulatieparen) 0 worden (dus eigen elektrodes=0)
     AmatERsvis2(AmatERsvis ==-1)=0;
%     % zorgen dat elektrodes die niet gestimuleerd zijn, ook 0 worden
     nietstim = setdiff(1:size(AmatERsvis,2),stimelek);
     AmatERsvis2(nietstim,:)=0;
     AmatERsvis2(:,nietstim)=0;
    
    % dubbel check dat overlappende elektrodes ook 0 worden en niet meedoen in metingen
     AmatERsvis2(:,chanoverlap)=0;
     AmatERsvis2(chanoverlap,:)=0;
    
    % graph maken zodat ik de BC kan gaan uitrekenen
    G = digraph(AmatERsvis2);
    
    % totaal aantal gestimuleerde elektroden
    stimelektot = size(ccep.ch,1);
    % totaal aantal trials
    trialtot = size(stimchan,1);
    
    % hoe vaak elke elektrode wordt gestimuleerd
    for el=1:size(AmatERsvis2,1)
        trialelek(el) = size(find(stimchan==el),1);
    end
    
    % totaal aantal mogelijke verbindingen
    for el=1:size(AmatERsvis2,1)
        % UIT
        % = trials waarin de elektrode zit * (totaal aantal elektroden - 2 (elektroden in stimpaar)
        n_outtot(el) = trialelek(el)*(stimelektot-2);
        
        % IN
        % = 2 (elektroden in een stimpaar) * (totaal aantal trials - aantal trials waarin de elektrode zit)
        n_intot(el) = 2*(trialtot - trialelek(el));
        
    end
    
    % NETWERKMATEN
    BC = centrality(G,'betweenness','Cost',G.Edges.Weight);
    
    for el=1:size(AmatERsvis2,1)
        % indegree en outdegree = som van 1tjes
        indegree(el) = sum(AmatERsvis2(:,el));
        indegreenorm(el) = indegree(el)/n_intot(el);
        outdegree(el) = sum(AmatERsvis2(el,:));
        outdegreenorm(el) = outdegree(el)/n_outtot(el);
        BCnorm(el) = BC(el)/(n_intot(el)*n_outtot(el));
    end
    
end