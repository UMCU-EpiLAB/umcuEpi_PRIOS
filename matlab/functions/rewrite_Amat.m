% Rewrite the adjacency matrix [respElec x stim pair] to an adjacency
% matrix [respElec x stimElec]. This is necessary to enable analysis of the
% network characteristics indegree, outdegree and betweenness centrality.

% INPUT:
% - spes
%   this is a structure, containing the following fields:
%   - cc_stimsets
%   matrix[stim pairs x 2] containing all the electrode numbers that are
%   stimulated (C1-C2 and C2-C1 are combined)
%   - ch
%   cell[channels x 1] containing the channel names 
%   - n1_peak_sample_check
%   matrix[channels x stimulations] containing latency of the significant
%   n1 peaks that are detected and visually checked.

% OUTPUT:
% - elecAmat
%   matrix [respElec x stimElec] containing the 0's and 1's (= 1 when a
%   CCEP is observed in a response electrode after stimulating another
%   electrode).


function elecAmat = rewrite_Amat(spes)

%% check that n1_peak_samples = [respElec x stimsets]
if size(spes.n1_peak_sample_check,1) == size(spes.ch,1) && ...
    size(spes.n1_peak_sample_check,2) == size(spes.cc_stimsets,1)
     
    % do nothing because the size is correct
    spesAmat = ~isnan(spes.n1_peak_sample_check);

elseif size(spes.n1_peak_sample_check,2) == size(spes.ch,1) && ...
    size(spes.n1_peak_sample_check,1) == size(spes.cc_stimsets,1)
    
    spesAmat = ~isnan(spes.n1_peak_sample_check');

else
    error('The size of n1-peak_sample does not correspond with the number of channels and stimulus pairs')
end

%% 

% Pre-allocation
elecAmat = zeros(size(spesAmat,1)); % [respElec x stimElec] 

for nStim = 1:size(spesAmat,2)
    stimElec = spes.cc_stimsets(nStim,:);

    elecAmat(:,stimElec) = elecAmat(:,stimElec) + spesAmat(:,nStim);    

end

end
