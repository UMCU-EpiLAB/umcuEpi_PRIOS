function agreement_parameter = agreement_parameters(agreement, dataBase_clin, dataBase_prop, myDataPath)
    agreement_parameter = struct;    
    wantedAmatClin = agreement.AmatClin'; 
    wantedAmatProp = agreement.AmatProp';
    
ERs_stimpClin = sum(wantedAmatClin,2);
ERs_stimpProp = sum(wantedAmatProp,2);    

% total number of ERs evoked in an electrode
ERs_elecClin = sum(wantedAmatClin,1);
ERs_elecProp = sum(wantedAmatProp,1);
    
      
%% Work with the electrodes in the columns and rows 
% All stimulations (10)    
    elec_matClin = dataBase_clin.elec_Amat;
    elec_matClin(elec_matClin ==2) = 1;
    
    G_Clin = digraph(elec_matClin);
    wbc_Clin = centrality(G_Clin,'betweenness');        % 'Cost' edge weights specify the length of the edges.
    Indegree_Clin = centrality(G_Clin,'indegree');
    Outdegree_Clin = centrality(G_Clin,'outdegree');
    
    %edges_Clin = linspace(min(Indegree_Clin),max(Indegree_Clin),7);            % create 7 equally-spaced bins based on their indegree score (number of detected ERs)
    %bins_Clin = discretize(Indegree_Clin,edges_Clin);
    
% Only 2 stimulations, first pos and first neg direction
    elec_matProp = dataBase_prop.elec_Amat;
    elec_matProp(elec_matProp == 2) = 1;
    
    G_Prop = digraph(elec_matProp);
    wbc_Prop = centrality(G_Prop,'betweenness');        % 'Cost' edge weights specify the length of the edges.
    Indegree_Prop = centrality(G_Prop,'indegree');
    Outdegree_Prop = centrality(G_Prop,'outdegree');
    
   % edges_Prop = linspace(min(Indegree_Prop),max(Indegree_Prop),7);            % create 7 equally-spaced bins based on their indegree score
    %bins_Prop = discretize(Indegree_Prop,edges_Prop);
     

%     high_bins_Clin = max(bins_Clin);
%     highest_ind_Clin = find(bins_Clin == high_bins_Clin);
%     high_bins_Prop = max(bins_Prop);
%     highest_ind_Prop = find(bins_Prop == high_bins_Prop);
% 
%     edges_o_Clin = linspace(min(Outdegree_Clin),max(Outdegree_Clin),7);            % create 7 equally-spaced bins based on their indegree score
%     bins_o_Clin = discretize(Outdegree_Clin,edges_o_Clin);
%     high_bins_o_Clin = max(bins_o_Clin);
%     highest_outd_Clin = find(bins_o_Clin ==high_bins_o_Clin);
% 
%     edges_o_Prop = linspace(min(Outdegree_Prop),max(Outdegree_Prop),7);            % create 7 equally-spaced bins based on their indegree score
%     bins_o_Prop = discretize(Outdegree_Prop,edges_o_Prop);
%     high_bins_o_Prop = max(bins_o_Prop);
%     highest_outd_Prop = find(bins_o_Prop ==high_bins_o_Prop);

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
  agreement_parameter.outdegreeN_Clin = outdegreenormClin;
  agreement_parameter.BCN_Clin = BCnormClin;
  agreement_parameter.ERs_stimpClin = ERs_stimpClin;
  agreement_parameter.ERs_elecClin = ERs_elecClin;
  
  agreement_parameter.indegreeN_Prop = indegreenormProp;
  agreement_parameter.outdegreeN_Prop = outdegreenormProp;
  agreement_parameter.BCN_Prop = BCnormProp;
  agreement_parameter.ERs_stimpProp = ERs_stimpProp;
  agreement_parameter.ERs_elecProp = ERs_elecProp;
  
%   agreement_parameter.highest_ind_Prop = highest_ind_Prop;
%   agreement_parameter.highest_ind_Clin = highest_ind_Clin;
%   agreement_parameter.highest_outd_Prop = highest_outd_Prop;
%   agreement_parameter.highest_outd_Clin = highest_outd_Clin;

  % All variables are also saved to a excel variant
%   write2excelTables(dataBase_clin, myDataPath, agreement_parameter);       % database is only needed for the electrode names and the stimnames, therefore does not matter whether database 10 or 2 is taken.
  
end


