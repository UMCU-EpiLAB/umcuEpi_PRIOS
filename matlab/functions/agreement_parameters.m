function spes = calcNetworkMeasures(spes)
% Determine the network characteristics/network parameters per
% electrode/stimulation pair. 
% The rewritten adjacency matrix of function rewrite_Amat.m is used to
% generate an adjacency matrix with electrodes in columns and rows used for
% the in-degree, out-degree and BC per electrode.

% Pre-allocation
% agreement_parameter = struct;   

% Original adjacency matrix with stimultion pairs in rows and electrodes in
% columns
% wantedAmatClin = Amat.AmatClin'; 
% wantedAmatProp = Amat.AmatProp';

% Determine the total number of ERs per stimulation pair (row)
% ERs_stimpClin = sum(wantedAmatClin,2);
% ERs_stimpProp = sum(wantedAmatProp,2);    

% Total number of ERs evoked in an electrode 
% (later the outdegree is determined which is also a measure for number of ERs evoked per electrode)
% ERs_elecClin = sum(wantedAmatClin,1);
% ERs_elecProp = sum(wantedAmatProp,1);
    
      
%% Work with the electrodes in the columns and rows 
% The rewritten adjacency matrix of function rewrite_Amat.m is used.
% The betweenness centrality is a measure of how often an electrode is a
% bridge between other electrodes. Electrodes with a high BC are often 
% important controllers of signals because it connects multiple areas of the brain.
% In-degree is a measure for the number of responses detected on an electrode evoked by other electrodes. 
% Out-degree is a measure of the number of responses evoked by that electrode.

elec_mat = spes.elec_Amat;

% The network parameters are normalised later (line 58), it is therefore not
% necessary to know whether an electrode evoked ERs when present in
% both stimulation pairs. 
elec_mat(elec_mat == 2) = 1;

% Digraph finds where each row is connected with a column (where a 1 is filled in)
G = digraph(elec_mat);                     
bc = centrality(G,'betweenness');                            
indegree = centrality(G,'indegree');   % Equal to connections per column of elec_matClin 
outdegree = centrality(G,'outdegree'); % Equal to connections per row of elec_matClin 


%% Normalisation %%
% The parameters are normalised by considering the number of stimulation 
% pairs an electrode is part of in each electrode, as suggested van Blooijs (2018)

% Total number of electrodes
elec_tot = size(spes.n1_peak_sample,1);
% Total number of stimpairs
stimp_tot = size(spes.n1_peak_sample,2);

% Number of times an electrode is part of a stimulation pair  
times_elec_in_stimp = zeros(1,elec_tot);
for el = 1:elec_tot
    times_elec_in_stimp(el) = size(find(spes.stimsets_avg==el),1);
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
outdegreenorm = zeros(1,elec_tot); indegreenorm = zeros(1,elec_tot); BCnorm = zeros(1,elec_tot); 

% Perform the normalisation
for el = 1:elec_tot

    outdegreenorm(el) = outdegree(el)/n_outtot(el);
    indegreenorm(el) = indegree(el)/n_intot(el);
    BCnorm(el) = bc(el)/(n_intot(el)*n_outtot(el));
    
end

% Save the values in an struct
spes.indegreeN = indegreenorm;
spes.outdegreeN = outdegreenorm;
spes.BCN = BCnorm;
spes.ERs_stimp = ERs_stimpClin;
spes.ERs_elec = ERs_elecClin;


end


