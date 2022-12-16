% Calculate the network measures "indegree", "outdegree" and "betweenness
% centrality" per electrode. The betweenness centrality (BC) is a measure
% of how often an electrode is a bridge between other electrodes.
% Electrodes with a high BC are often important controllers of signals,
% because it connects multiple areas of the brain. The indegree is a
% measure for the number of responses in a response channel evoked after
% stimulating other electrodes. The outdegree is a measure of the number of
% responses evoked after stimulating that electrode.
%
% NORMALISATION: We normalised the network measures to enable comparison
% between several subjects. We divided the outdegree by the possible number
% of CCEPs that can be evoked by the stimulated electrode (outgoing
% connections). We dividied the indegree by the possible number of
% stimulated electrodes that could have evoked a CCEP in the response
% electrode (ingoing connections). We divided the BC by the multiplication
% of the total number of possible outgoing and ingoing connections. 
% 
% More information can be found here: van Blooijs D, Leijten FSS, van Rijen
% PC, Meijer HGE, Huiskamp GJM. Evoked directional network characteristics
% of epileptogenic tissue derived from single pulse electrical stimulation.
% Hum Brain Mapp. 2018 Nov;39(11):4611-4622. doi: 10.1002/hbm.24309. 

% INPUT:
% - spes
%   this is a structure, containing the following fields:
% - elecAmat
%   matrix [respElec x stimElec] containing the 0's and 1's (= 1 or higher
%   when a CCEP is observed in a response electrode after stimulating
%   another electrode). 
% - cc_stimsets matrix[stim pairs x 2] containing all
%   the electrode numbers that are stimulated (C1-C2 and C2-C1 are
%   combined) - ch cell[channels x 1] containing the channel names

% OUTPUT:
% - spes
%   a structure to which the following fields are added in this function:
% - indegreeNorm
% - outdegreeNorm
% - bcNorm

function spes = calcNetworkMeasures(spes)

%% Process adjacency matrix
% make a binary matrix which does not show how often a CCEP is evoked in a
% certain response Channel after stimulating another channel, but whether a
% CCEP is evoked. 

elecAmat = spes.elecAmat;
elecAmat(elecAmat >= 1) = 1;

%% Calculate network measures
% make a graph (G) and calculate the indegree, outdegree en betweenness
% centrality (BC)

G = digraph(elecAmat);                     

indegree = centrality(G,'indegree');   
outdegree = centrality(G,'outdegree'); 
bc = centrality(G,'betweenness');                            

%% Normalisation

% Total number of electrodes
nChanTot = size(spes.ch,1);
% Total number of stimpairs
nStimPairTot = size(spes.cc_stimsets,1);

% Number of times an electrode is part of a stimulation pair  
timesChanInStimp = zeros(nChanTot,1);
for nChan = 1:nChanTot
    timesChanInStimp(nChan,:) = size(find(spes.cc_stimsets == nChan),1);
end

% Total possible outgoing connections for each stimulated electrode: Each
% electrode that is part of 2 stimulation pairs (timesChanInStimp = 2, for
% example C02 in C01-C02 and C02-C03) could, in theory, evoke twice as many
% reponses
nOutTot = timesChanInStimp*(nChanTot-2); % Minus 2 because a stim-pair has 2 electrodes in which a response cannot be evoked

% Total possible ingoing connections for each response electrode: In
% theory, all stimulation pairs minus the pairs in which the specific
% response electrode is stimulated, could evoke a response in the response
% electrode.
nInTot = 2*(nStimPairTot - timesChanInStimp); % Times 2 because each stim-pair had 2 electrodes and each electrode could have led to a response

% Apply the normalisation
outdegreeNorm = outdegree./nOutTot;
indegreeNorm = indegree./nInTot;
bcNorm = bc./(nInTot.*nOutTot);
    
%% Write variable back to struct
spes.indegreeNorm = indegreeNorm;
spes.outdegreeNorm = outdegreeNorm;
spes.bcNorm = bcNorm;

end


