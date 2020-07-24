function agreement_parameter = agreement_parameters(agreement, dataBase10, dataBase2,stimchans)
    agreement_parameter = struct;    
    wantedAmat10 = agreement.Amat10'; 
    wantedAmat2 = agreement.Amat2';
    Amat_10 = zeros(size(wantedAmat10,1)+size(wantedAmat10,2));     % Matrix with 4 quadrants, left under must be filled with the adjacency matrix            
    Amat_2 = zeros(size(wantedAmat2,1)+size(wantedAmat2,2));        % Matrix with 4 quadrants, left under must be filled with the adjacency matrix            

    rowStart = size(wantedAmat10,2)+1;      % Same for 2 stims
    columnEnd = size(wantedAmat10,2);
    ElecName = dataBase10.ch;
       
    Amat_10((rowStart:end),(1:columnEnd)) = wantedAmat10;
    Amat_2((rowStart:end),(1:columnEnd)) = wantedAmat2;

    for stimPair = size(wantedAmat10,2)+1:size(Amat_10,1)
        ERs_stimp10((stimPair-(rowStart-1)),1) = sum(Amat_10(stimPair,:));       % total number of ERs evoked by a stimulation pair 
        ERs_stimp2((stimPair-(rowStart-1)),1) = sum(Amat_2(stimPair,:));
    end
       
    ColN = {'number of ERs'};
    rowNames = stimchans;
    rank_stimp10 = array2table(ERs_stimp10,'RowNames',rowNames,'VariableNames',ColN);   % Number of outgoing stimulations (ERs evoked per stimpair)
    rank_stimp10 = sortrows(rank_stimp10,1,'descend');
    rank_stimp2 = array2table(ERs_stimp2,'RowNames',rowNames,'VariableNames',ColN);   
    rank_stimp2 = sortrows(rank_stimp2,1,'descend');
    
    for elec = 1:size(wantedAmat10,2)
        ERs_elec10(1,elec) = sum(Amat_10(:,elec));                            % total number of ERs evoked in an electrode (el) 
        ERs_elec2(1,elec) = sum(Amat_2(:,elec)); 
    end

    ERs_elec10 = ERs_elec10';
    ERs_elec2 = ERs_elec2';
    RowNames = ElecName;
    rank_elec10 = array2table(ERs_elec10,'RowNames',RowNames,'VariableNames',ColN);   % Indegree
    rank_elec10 = sortrows(rank_elec10,1,'descend');
    rank_elec2 = array2table(ERs_elec2,'RowNames',RowNames,'VariableNames',ColN);   % Indegree
    rank_elec2 = sortrows(rank_elec2,1,'descend');
      
%% Work with the electrodes in the columns and rows 
% All stimulations (10)    
    elec_mat10 = dataBase10.elec_Amat;
    G_10 = digraph(elec_mat10);
    wbc_10 = centrality(G_10,'betweenness','Cost',G_10.Edges.Weight);        % 'Cost' edge weights specify the length of the edges.
    hub_ranks_10 = centrality(G_10,'hubs');
    Indegree_10 = centrality(G_10,'indegree','Importance',G_10.Edges.Weight);
    Outdegree_10 = centrality(G_10,'outdegree','Importance',G_10.Edges.Weight);
    
    edges_10 = linspace(min(Indegree_10),max(Indegree_10),7);            % create 7 equally-spaced bins based on their indegree score (number of detected ERs)
    bins_10 = discretize(Indegree_10,edges_10);
    
% Only 2 stimulations, first pos and first neg direction
    elec_mat2 = dataBase2.elec_Amat;
    G_2 = digraph(elec_mat2);
    wbc_2 = centrality(G_2,'betweenness','Cost',G_2.Edges.Weight);        % 'Cost' edge weights specify the length of the edges.
    hub_ranks_2 = centrality(G_2,'hubs');
    Indegree_2 = centrality(G_2,'indegree','Importance',G_2.Edges.Weight);
    Outdegree_2 = centrality(G_2,'outdegree','Importance',G_2.Edges.Weight);
    
    edges_2 = linspace(min(Indegree_2),max(Indegree_2),7);            % create 7 equally-spaced bins based on their indegree score
    bins_2 = discretize(Indegree_2,edges_2);
     
    figure('Position',[1,2,1600,756])
    subplot(1,2,1, 'Position',[0.018,0.11,0.48,0.82])
    p_in10 = plot(G_10,'Layout','force','NodeLabel',dataBase10.ch,'NodeColor','r','MarkerSize',bins_10*2,'NodeFontSize',10);
    high_bins_10 = max(bins_10);
    highest_ind_10 = find(bins_10 == high_bins_10);
    highlight(p_in10,highest_ind_10,'NodeColor','g')
    title({'{\bf\fontsize{13} Betweenness centrality all, highest indegree electrodes}'; 'Size of markers indicates the indegree ranking'},'FontWeight','Normal')
    
    subplot(1,2,2,'Position',[0.52,0.11,0.48,0.81])
    p_in2 = plot(G_2,'Layout','force','NodeLabel',dataBase2.ch,'NodeColor','r','MarkerSize',bins_2*2,'NodeFontSize',10);
    high_bins_2 = max(bins_2);
    highest_ind_2 = find(bins_2 == high_bins_2);
    highlight(p_in2,highest_ind_2,'NodeColor','g')
    title({'{\bf\fontsize{13} Betweenness centrality 2 stims, highest indegree electrodes}'; 'Size of markers indicates the indegree ranking'},'FontWeight','Normal')
    
    
    figure('Position',[1,2,1600,756])
    subplot(1,2,1, 'Position',[0.018,0.11,0.48,0.82])
    edges_o_10 = linspace(min(Outdegree_10),max(Outdegree_10),7);            % create 7 equally-spaced bins based on their indegree score
    bins_o_10 = discretize(Outdegree_10,edges_o_10);
    p_o_10 = plot(G_10,'Layout','force','NodeLabel',dataBase10.ch,'NodeColor','r','MarkerSize',bins_o_10*2,'NodeFontSize',10);
    high_bins_o_10 = max(bins_o_10);
    highest_ind_o_10 = find(bins_o_10 ==high_bins_o_10);
    highlight(p_o_10,highest_ind_o_10,'NodeColor','g')
    title({'{\bf\fontsize{13} Betweenness centrality all stims, highest outdegree electrodes}'; 'Size of markers indicates the outdegree ranking'},'FontWeight','Normal')
    
    subplot(1,2,2,'Position',[0.52,0.11,0.48,0.81])
    edges_o_2 = linspace(min(Outdegree_2),max(Outdegree_2),7);            % create 7 equally-spaced bins based on their indegree score
    bins_o_2 = discretize(Outdegree_2,edges_o_2);
    p_o_2 = plot(G_2,'Layout','force','NodeLabel',dataBase2.ch,'NodeColor','r','MarkerSize',bins_o_2*2,'NodeFontSize',10);
    high_bins_o_2 = max(bins_o_2);
    highest_ind_o_2 = find(bins_o_2 ==high_bins_o_2);
    highlight(p_o_2,highest_ind_o_2,'NodeColor','g')
    title({'{\bf\fontsize{13} Betweenness centrality 2 stims, highest outdegree electrodes}'; 'Size of markers indicates the outdegree ranking'},'FontWeight','Normal')
 
    figure('Position',[1,2,1600,756])
    subplot(1,2,1, 'Position',[0.018,0.11,0.48,0.82])
    p10 = plot(G_10,'Layout','force','NodeLabel',dataBase10.ch,'NodeColor','r','NodeFontSize',10);
    title({'{\bf\fontsize{13}Betweenness centrality all stims}'; 'In green the electrode with a high ranking in indegree and outdegree'},'FontWeight','Normal');
    for i = 1:length(dataBase10.ch)
        if ismember(i, highest_ind_10, 'rows') && ismember(i,highest_ind_o_10,'rows')               % When the electrode is highest ranked in the indegree and in outdegree
            highlight(p10,i,'NodeColor','g','MarkerSize',10)
        end
    end
    
    subplot(1,2,2,'Position',[0.52,0.11,0.48,0.81])
    p2 = plot(G_2,'Layout','force','NodeLabel',dataBase2.ch,'NodeColor','r','NodeFontSize',10);
    title({'{\bf\fontsize{13}Betweenness centrality of 2 stims}'; 'In green the electrode with a high ranking in indegree and outdegree'},'FontWeight','Normal');
    for i = 1:length(dataBase2.ch)
        if ismember(i, highest_ind_2, 'rows') && ismember(i,highest_ind_o_2,'rows')               % When the electrode is highest ranked in the indegree and in outdegree
            highlight(p2,i,'NodeColor','g','MarkerSize',10)
        end
    end

%%% Normalisation %%%
    % total number of electrodes
    stimelektot = size(dataBase10.ch,1);
    % total number of stimpairs
    stimptot = size(dataBase10.stimpnames_avg,2);
    
    % Number of times an electrode is used in a stimulation pair  
    trialelek = zeros(1,stimelektot);
    for el=1:size(elec_mat10,1)
        trialelek(el) = size(find(dataBase10.cc_stimsets_avg==el),1);
    end
    
    % totaal number of possible connections
    n_outtot = zeros(1,stimelektot);
    n_intot = zeros(1,stimelektot);
    for el=1:size(elec_mat10,1)
        n_outtot(el) = trialelek(el)*(stimelektot-2);
        n_intot(el) = 2*(stimptot - trialelek(el)); 
    end
    
    % NETWERKMATEN   
    outdegreenorm10 = zeros(1,stimelektot);
    indegreenorm10 = zeros(1,stimelektot);
    BCnorm10 = zeros(1,stimelektot); 
    
    outdegreenorm2 = zeros(1,stimelektot);
    indegreenorm2 = zeros(1,stimelektot);
    BCnorm2 = zeros(1,stimelektot); 
    
    for el=1:size(elec_mat10,1)
        outdegreenorm10(el) = Outdegree_10(el)/n_outtot(el);
        indegreenorm10(el) = Indegree_10(el)/n_intot(el);
        BCnorm10(el) = wbc_10(el)/(n_intot(el)*n_outtot(el));
        
        outdegreenorm2(el) = Outdegree_2(el)/n_outtot(el);
        indegreenorm2(el) = Indegree_2(el)/n_intot(el);
        BCnorm2(el) = wbc_2(el)/(n_intot(el)*n_outtot(el));
    end
    
  agreement_parameter.indegreeN_10 = indegreenorm10;
  agreement_parameter.outdegreeN_10 = outdegreenorm10;
  agreement_parameter.BCN_all_10 = BCnorm10;
  agreement_parameter.rank_stimp_10 = rank_stimp10;
  agreement_parameter.rank_elec_10 = rank_elec10;
   
  agreement_parameter.indegreeN_2 = indegreenorm2;
  agreement_parameter.outdegreeN_2 = outdegreenorm2;
  agreement_parameter.BCN_2 = BCnorm2;
  agreement_parameter.rank_stimp_2 = rank_stimp2;
  agreement_parameter.rank_elec_2 = rank_elec2;
  
  
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
    