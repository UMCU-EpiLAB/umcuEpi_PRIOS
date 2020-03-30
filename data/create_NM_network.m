function strout=create_NM_network(nodes, conmat, constrength)
%strout=create_NM_network(nodes, conmat, constrength)
%
% This function creates a struct containing all parameters of a network of
% neural masses
%
% Input:
%   nodes: 1xN vector of NM parameter structures.
%   conmat: NxN matrix containing the connectivity between nodes. Here conmat(i,j)
%   is the connectivity from j to i.
%   constrength: (double) global connectivity strength
%
% Output:
%   strout: a struct containing all parameters of a network of NMs
%
% This code is part of the simulation code of the manuscript 'Pathological responses to single pulse electrical stimuli in epilepsy: the role of feedforward inhibition'
% (c) 2019 Jurgen Hebbink (University of Twente, University Medical Centre Utrecht)

strout.nodes=nodes;
strout.conmat=conmat;
strout.constrength=constrength;
end

