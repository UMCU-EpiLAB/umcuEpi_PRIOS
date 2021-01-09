function agreement_parameter = agreement_parameters(agreement, dataBase_clin, dataBase_prop)
% Determine the network characteristics/network parameters per
% electrode/stimulation pair. 
% The rewritten adjacency matrix of function rewrite_Amat.m is used to
% generate adjacency matrix with electrodes in columns and rows used for
% the in-degree, out-degree and BC per electrode.

% Pre-allocation
agreement_parameter = struct;   
% Original adjacency matrix with stimultion pairs in rows and electrodes in
% columns
wantedAmatClin = agreement.AmatClin'; 
wantedAmatProp = agreement.AmatProp';

% Determine the total number of ERs per stimulation pair (row)
ERs_stimpClin = sum(wantedAmatClin,2);
ERs_stimpProp = sum(wantedAmatProp,2);    

% Total number of ERs evoked in an electrode 
% (later the outdegree is determined which is also a measure for number of ERs evoked per electrode)
ERs_elecClin = sum(wantedAmatClin,1);
ERs_elecProp = sum(wantedAmatProp,1);
    
      
%% Work with the electrodes in the columns and rows 
% The rewritten adjacency matrix of function rewrite_Amat.m is used
% betweenness centrality is is a measure of how often an electrode is a
% bridge between other electrodes. Electrodes with a high BC are often 
% important controllers of signals because it connects multiple areas of the brain.
% In-degree is a measure for the number of responses detected on an electrode evoked by other electrodes. 
% Out-degree is a measure of the number of responses evoked by that electrode.

% Clinical-SPES  
    elec_matClin = dataBase_clin.elec_Amat;
    
    % The network parameters are later normalised, it is therefore not
    % necessary to know whether an electrode evoked ERs when present in
    % both stimulation pairs. 
    elec_matClin(elec_matClin ==2) = 1;
    
    % Digraph creates a weighted directed graph
    G_Clin = digraph(elec_matClin);
    wbc_Clin = centrality(G_Clin,'betweenness');                            
    Indegree_Clin = centrality(G_Clin,'indegree');
    Outdegree_Clin = centrality(G_Clin,'outdegree');
     
% Propofol-SPES
    elec_matProp = dataBase_prop.elec_Amat;
    elec_matProp(elec_matProp == 2) = 1;
    
    G_Prop = digraph(elec_matProp);
    wbc_Prop = centrality(G_Prop,'betweenness');        
    Indegree_Prop = centrality(G_Prop,'indegree');
    Outdegree_Prop = centrality(G_Prop,'outdegree');
    
 
%%% Normalisation %%%
% The parameters are normalised by considering the number of stimulation 
% pairs an electrode is part of in each electrode, as suggested van Blooijs (2018)

    % Total number of electrodes
    stimelektot = size(dataBase_clin.ch,1);
    % Total number of stimpairs
    stimptot = size(dataBase_clin.stimpnames_avg,2);
    
    % Number of times an electrode is part of a stimulation pair  
    trialelek = zeros(1,stimelektot);
    for el=1:size(elec_matClin,1)
        trialelek(el) = size(find(dataBase_clin.stimsets_avg==el),1);
    end
    
    % Totaal number of possible connections
    n_outtot = zeros(1,stimelektot);
    n_intot = zeros(1,stimelektot);
    for el=1:size(elec_matClin,1)
        n_outtot(el) = trialelek(el)*(stimelektot-2);
        n_intot(el) = 2*(stimptot - trialelek(el)); 
    end
    
    % Network characteristics pre-allocation  
    outdegreenormClin = zeros(1,stimelektot);
    indegreenormClin = zeros(1,stimelektot);
    BCnormClin = zeros(1,stimelektot); 
    outdegreenormProp = zeros(1,stimelektot);
    indegreenormProp = zeros(1,stimelektot);
    BCnormProp = zeros(1,stimelektot); 
    
    % Perform the normalisation
    for el=1:size(elec_matClin,1)
        outdegreenormClin(el) = Outdegree_Clin(el)/n_outtot(el);
        indegreenormClin(el) = Indegree_Clin(el)/n_intot(el);
        BCnormClin(el) = wbc_Clin(el)/(n_intot(el)*n_outtot(el));
        
        outdegreenormProp(el) = Outdegree_Prop(el)/n_outtot(el);
        indegreenormProp(el) = Indegree_Prop(el)/n_intot(el);
        BCnormProp(el) = wbc_Prop(el)/(n_intot(el)*n_outtot(el));
    end
    
  % Save the values in an struct
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
    
end


