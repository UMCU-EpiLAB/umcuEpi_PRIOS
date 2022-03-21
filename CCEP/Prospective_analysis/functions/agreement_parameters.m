function agreement_parameter = agreement_parameters(agreement, dataBase_clin, dataBase_prop)
% Determine the network characteristics/network parameters per
% electrode/stimulation pair. 
% The rewritten adjacency matrix of function rewrite_Amat.m is used to
% generate an adjacency matrix with electrodes in columns and rows used for
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
% The rewritten adjacency matrix of function rewrite_Amat.m is used.
% The betweenness centrality is a measure of how often an electrode is a
% bridge between other electrodes. Electrodes with a high BC are often 
% important controllers of signals because it connects multiple areas of the brain.
% In-degree is a measure for the number of responses detected on an electrode evoked by other electrodes. 
% Out-degree is a measure of the number of responses evoked by that electrode.

% Clinical-SPES  
elec_matClin = dataBase_clin.elec_Amat;

% The network parameters are normalised later (line 58), it is therefore not
% necessary to know whether an electrode evoked ERs when present in
% both stimulation pairs. 
elec_matClin(elec_matClin ==2) = 1;

% Digraph finds where each row is connected with a column (where a 1 is filled in)
G_Clin = digraph(elec_matClin);                     
wbc_Clin = centrality(G_Clin,'betweenness');                            
Indegree_Clin = centrality(G_Clin,'indegree');   % Equal to connections per column of elec_matClin 
Outdegree_Clin = centrality(G_Clin,'outdegree'); % Equal to connections per row of elec_matClin 
 
% Propofol-SPES
elec_matProp = dataBase_prop.elec_Amat;
elec_matProp(elec_matProp == 2) = 1;

G_Prop = digraph(elec_matProp);
wbc_Prop = centrality(G_Prop,'betweenness');        
Indegree_Prop = centrality(G_Prop,'indegree');
Outdegree_Prop = centrality(G_Prop,'outdegree');


%% Normalisation %%
% The parameters are normalised by considering the number of stimulation 
% pairs an electrode is part of in each electrode, as suggested van Blooijs (2018)

% Total number of electrodes
elec_tot = size(dataBase_clin.n1_peak_sample,1);
% Total number of stimpairs
stimp_tot = size(dataBase_clin.n1_peak_sample,2);

% Number of times an electrode is part of a stimulation pair  
times_elec_in_stimp = zeros(1,elec_tot);
for el = 1:elec_tot
    times_elec_in_stimp(el) = size(find(dataBase_clin.stimsets_avg==el),1);
end

% Totaal number of possible connections
n_outtot = zeros(1,elec_tot);
n_intot = zeros(1,elec_tot);
for el = 1:elec_tot
    % Total possible outgoing responses for electrode (el)
    % Each electrode that is part of 2 stimulation pairs (times_elec_in_stimp = 2) could, in theory,
    % evoke twice as many reponses
    n_outtot(el) = times_elec_in_stimp(el)*(elec_tot-2);        % Minus 2 because a stim-pair has 2 electrodes that cannot receive a response
    
    % Total possible incoming responses for electrode (el).
    % In theory, all stimulation pairs minus the ones that electrode (el)
    % is part of, could evoke a response on electrode (el). 
    n_intot(el) = 2*(stimp_tot - times_elec_in_stimp(el));      % Times 2 because each stim-pair had 2 electrodes and each electrode could have led to a response
end

% Pre-allocation Network characteristics  
outdegreenormClin = zeros(1,elec_tot); indegreenormClin = zeros(1,elec_tot); BCnormClin = zeros(1,elec_tot); outdegreenormProp = zeros(1,elec_tot); indegreenormProp = zeros(1,elec_tot); BCnormProp = zeros(1,elec_tot); 

% Perform the normalisation
for el = 1:elec_tot
    % Clinical SPES
    outdegreenormClin(el) = Outdegree_Clin(el)/n_outtot(el);
    indegreenormClin(el) = Indegree_Clin(el)/n_intot(el);
    BCnormClin(el) = wbc_Clin(el)/(n_intot(el)*n_outtot(el));
    
    % Propofol SPES
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


