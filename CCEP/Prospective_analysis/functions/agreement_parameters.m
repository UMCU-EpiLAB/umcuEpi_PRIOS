function agreement_parameter = agreement_parameters(agreement, dataBase_clin, dataBase_prop, myDataPath)
    agreement_parameter = struct;    
    wantedAmatClin = agreement.AmatClin'; 
    wantedAmatProp = agreement.AmatProp';
    Amat_clin = zeros(size(wantedAmatClin,1)+size(wantedAmatClin,2));     % Matrix with 4 quadrants, left under must be filled with the adjecency matrix            
    Amat_prop = zeros(size(wantedAmatProp,1)+size(wantedAmatProp,2));        % Matrix with 4 quadrants, left under must be filled with the adjecency matrix            

    rowStart = size(wantedAmatClin,2)+1;      % Same for 2 stims
    columnEnd = size(wantedAmatClin,2);
    ElecName = dataBase_clin.ch;
       
    Amat_clin((rowStart:end),(1:columnEnd)) = wantedAmatClin;
    Amat_prop((rowStart:end),(1:columnEnd)) = wantedAmatProp;

    for stimPair = size(wantedAmatClin,2)+1:size(Amat_clin,1)
        ERs_stimpClin((stimPair-(rowStart-1)),1) = sum(Amat_clin(stimPair,:));       % total number of ERs evoked by a stimulation pair 
        ERs_stimpProp((stimPair-(rowStart-1)),1) = sum(Amat_prop(stimPair,:));
    end
       
    
    for elec = 1:size(wantedAmatClin,2)
        ERs_elecClin(1,elec) = sum(Amat_clin(:,elec));                            % total number of ERs evoked in an electrode (el) 
        ERs_elecProp(1,elec) = sum(Amat_prop(:,elec)); 
    end
      
%% Work with the electrodes in the columns and rows 
% All stimulations (10)    
    elec_matClin = dataBase_clin.elec_Amat;
    G_Clin = digraph(elec_matClin);
    wbc_Clin = centrality(G_Clin,'betweenness','Cost',G_Clin.Edges.Weight);        % 'Cost' edge weights specify the length of the edges.
%     hub_ranks_10 = centrality(G_Clin,'hubs');
    Indegree_Clin = centrality(G_Clin,'indegree','Importance',G_Clin.Edges.Weight);
    Outdegree_Clin = centrality(G_Clin,'outdegree','Importance',G_Clin.Edges.Weight);
    
    edges_Clin = linspace(min(Indegree_Clin),max(Indegree_Clin),7);            % create 7 equally-spaced bins based on their indegree score (number of detected ERs)
    bins_Clin = discretize(Indegree_Clin,edges_Clin);
    
% Only 2 stimulations, first pos and first neg direction
    elec_matProp = dataBase_prop.elec_Amat;
    G_Prop = digraph(elec_matProp);
    wbc_Prop = centrality(G_Prop,'betweenness','Cost',G_Prop.Edges.Weight);        % 'Cost' edge weights specify the length of the edges.
%     hub_ranks_2 = centrality(G_Prop,'hubs');
    Indegree_Prop = centrality(G_Prop,'indegree','Importance',G_Prop.Edges.Weight);
    Outdegree_Prop = centrality(G_Prop,'outdegree','Importance',G_Prop.Edges.Weight);
    
    edges_Prop = linspace(min(Indegree_Prop),max(Indegree_Prop),7);            % create 7 equally-spaced bins based on their indegree score
    bins_Prop = discretize(Indegree_Prop,edges_Prop);
     
%     figure('Position',[1,2,1600,756])
%     subplot(1,2,1, 'Position',[0.018,0.11,0.48,0.82])
%     p_in10 = plot(G_10,'Layout','force','NodeLabel',dataBase10.ch,'NodeColor','r','MarkerSize',bins_10*2,'NodeFontSize',10);
    high_bins_Clin = max(bins_Clin);
    highest_ind_Clin = find(bins_Clin == high_bins_Clin);
%     highlight(p_in10,highest_ind_10,'NodeColor','g')
%     title({'{\bf\fontsize{13} Betweenness centrality all, highest indegree electrodes}'; 'Size of markers indicates the indegree ranking'},'FontWeight','Normal')
%     
%     subplot(1,2,2,'Position',[0.52,0.11,0.48,0.81])
%     p_in2 = plot(G_2,'Layout','force','NodeLabel',dataBase2.ch,'NodeColor','r','MarkerSize',bins_2*2,'NodeFontSize',10);
    high_bins_Prop = max(bins_Prop);
    highest_ind_Prop = find(bins_Prop == high_bins_Prop);
%     highlight(p_in2,highest_ind_2,'NodeColor','g')
%     title({'{\bf\fontsize{13} Betweenness centrality 2 stims, highest indegree electrodes}'; 'Size of markers indicates the indegree ranking'},'FontWeight','Normal')
%     
    
%     figure('Position',[1,2,1600,756])
%     subplot(1,2,1, 'Position',[0.018,0.11,0.48,0.82])
    edges_o_Clin = linspace(min(Outdegree_Clin),max(Outdegree_Clin),7);            % create 7 equally-spaced bins based on their indegree score
    bins_o_Clin = discretize(Outdegree_Clin,edges_o_Clin);
%     p_o_10 = plot(G_10,'Layout','force','NodeLabel',dataBase10.ch,'NodeColor','r','MarkerSize',bins_o_10*2,'NodeFontSize',10);
    high_bins_o_Clin = max(bins_o_Clin);
    highest_outd_Clin = find(bins_o_Clin ==high_bins_o_Clin);
%     highlight(p_o_10,highest_outd_10,'NodeColor','g')
%     title({'{\bf\fontsize{13} Betweenness centrality all stims, highest outdegree electrodes}'; 'Size of markers indicates the outdegree ranking'},'FontWeight','Normal')
%     
%     subplot(1,2,2,'Position',[0.52,0.11,0.48,0.81])
    edges_o_Prop = linspace(min(Outdegree_Prop),max(Outdegree_Prop),7);            % create 7 equally-spaced bins based on their indegree score
    bins_o_Prop = discretize(Outdegree_Prop,edges_o_Prop);
%     p_o_2 = plot(G_2,'Layout','force','NodeLabel',dataBase2.ch,'NodeColor','r','MarkerSize',bins_o_2*2,'NodeFontSize',10);
    high_bins_o_Prop = max(bins_o_Prop);
    highest_outd_Prop = find(bins_o_Prop ==high_bins_o_Prop);
%     highlight(p_o_2,highest_outd_2,'NodeColor','g')
%     title({'{\bf\fontsize{13} Betweenness centrality 2 stims, highest outdegree electrodes}'; 'Size of markers indicates the outdegree ranking'},'FontWeight','Normal')

% figure('Position',[1,2,1600,756])
% subplot(1,2,1, 'Position',[0.018,0.11,0.48,0.82])
% p10 = plot(G_10,'Layout','force','NodeLabel',dataBase10.ch,'NodeColor','r','NodeFontSize',10);
% title({'{\bf\fontsize{13}Betweenness centrality all stims}'; 'In green the electrode with a high ranking in indegree and outdegree'},'FontWeight','Normal');
% for i = 1:length(dataBase10.ch)
%     if ismember(i, highest_ind_10, 'rows') && ismember(i,highest_outd_10,'rows')               % When the electrode is highest ranked in the indegree and in outdegree
%         highlight(p10,i,'NodeColor','g','MarkerSize',10)
%     end
% end
% 
% subplot(1,2,2,'Position',[0.52,0.11,0.48,0.81])
% p2 = plot(G_2,'Layout','force','NodeLabel',dataBase2.ch,'NodeColor','r','NodeFontSize',10);
% title({'{\bf\fontsize{13}Betweenness centrality of 2 stims}'; 'In green the electrode with a high ranking in indegree and outdegree'},'FontWeight','Normal');
% for i = 1:length(dataBase2.ch)
%     if ismember(i, highest_ind_2, 'rows') && ismember(i,highest_outd_2,'rows')               % When the electrode is highest ranked in the indegree and in outdegree
%         highlight(p2,i,'NodeColor','g','MarkerSize',10)
%     end
% end

%%% Normalisation %%%
    % total number of electrodes
    stimelektot = size(dataBase_clin.ch,1);
    % total number of stimpairs
    stimptot = size(dataBase_clin.stimpnames_avg,2);
    
    % Number of times an electrode is used in a stimulation pair  
    trialelek = zeros(1,stimelektot);
    for el=1:size(elec_matClin,1)
        trialelek(el) = size(find(dataBase_clin.stimsets_avg==el),1);
    end
    
    % totaal number of possible connections
    n_outtot = zeros(1,stimelektot);
    n_intot = zeros(1,stimelektot);
    for el=1:size(elec_matClin,1)
        n_outtot(el) = trialelek(el)*(stimelektot-2);
        n_intot(el) = 2*(stimptot - trialelek(el)); 
    end
    
    % NETWERKMATEN   
    outdegreenormClin = zeros(1,stimelektot);
    indegreenormClin = zeros(1,stimelektot);
    BCnormClin = zeros(1,stimelektot); 
    
    outdegreenormProp = zeros(1,stimelektot);
    indegreenormProp = zeros(1,stimelektot);
    BCnormProp = zeros(1,stimelektot); 
    
    for el=1:size(elec_matClin,1)
        outdegreenormClin(el) = Outdegree_Clin(el)/n_outtot(el);
        indegreenormClin(el) = Indegree_Clin(el)/n_intot(el);
        BCnormClin(el) = wbc_Clin(el)/(n_intot(el)*n_outtot(el));
        
        outdegreenormProp(el) = Outdegree_Prop(el)/n_outtot(el);
        indegreenormProp(el) = Indegree_Prop(el)/n_intot(el);
        BCnormProp(el) = wbc_Prop(el)/(n_intot(el)*n_outtot(el));
    end
    
  agreement_parameter.indegreeN_Clin = indegreenormClin;
  agreement_parameter.outdegreeN_CLin = outdegreenormClin;
  agreement_parameter.BCN_all_Clin = BCnormClin;
  agreement_parameter.ERs_stimpClin = ERs_stimpClin;
  agreement_parameter.ERs_elecClin = ERs_elecClin;
  
  agreement_parameter.indegreeN_Prop = indegreenormProp;
  agreement_parameter.outdegreeN_Prop = outdegreenormProp;
  agreement_parameter.BCN_Prop = BCnormProp;
  agreement_parameter.ERs_stimpProp = ERs_stimpProp;
  agreement_parameter.ERs_elecProp = ERs_elecProp;
  
  agreement_parameter.highest_ind_Prop = highest_ind_Prop;
  agreement_parameter.highest_ind_Clin = highest_ind_Clin;
  agreement_parameter.highest_outd_Prop = highest_outd_Prop;
  agreement_parameter.highest_outd_Clin = highest_outd_Clin;

  % All variables are also saved to a excel variant
  write2excelTables(dataBase_clin, myDataPath, agreement_parameter);       % database is only needed for the electrode names and the stimnames, therefore does not matter whether database 10 or 2 is taken.
  
end


