function strout=create_NM(A,a,B,b,G,g,C,sd,alpha,beta,gamma)
% strout=create_NM(A,a,B,b,G,g,C,sd,alpha,beta,gamma)
%
% This function create a struct containing all parameters of a neural mass
%
% Input:
%   A: Synaptic gain excitatory synapses
%   a: Time scale excitatory synapses
%   B: Synaptic gain slow inhibitory synapses
%   b: Time scale slow inhibitory synapses
%   G: Synaptic gain fast inhibitory synapses
%   g: Time scale fast inhibitory synapses
%   C: Internal connectivity strength
%   sd: Standard deviation of the noise
%   alpha: Scaling constant for external input to local excitatory cells
%   beta: Scaling  constant for external input to slow inhibitory cells
%   gamma: Scaling  constant for external input to fast inhibitory cells
%
%Output:
%   strout: a struct containing the parameters for the NM
%
% This code is part of the simulation code of the manuscript 'Pathological responses to single pulse electrical stimuli in epilepsy: the role of feedforward inhibition'
% (c) 2019 Jurgen Hebbink (University of Twente, University Medical Centre Utrecht)

strout.A=A;
strout.a=a;
strout.B=B;
strout.b=b;
strout.G=G;
strout.g=g;
strout.alpha=alpha;
strout.beta=beta;
strout.gamma=gamma;

strout.Ivar=@(x) 0*x; %Time varying input (default none)

%Internal connectivity strengths
strout.C1=C;        % py->ex
strout.C2=0.8*C;    % ex->py
strout.C3=0.25*C;   % py->in
strout.C4=0.25*C;   % in->py
strout.C5=0.3*C;    % py->if
strout.C6=0.1*C;    % in->if
strout.C7=1.2*C;    % if->py
strout.C8=0.0*C;    % if->if (Not used in the manuscript, self-loop fast inhibitory population)

strout.sd=sd;

end