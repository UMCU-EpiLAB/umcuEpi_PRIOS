%% Description init_NM_network.m
% This script is used to intialize and perform the simulation of two
% coupled Wendling neural masses which receive SPES simulation.
%
% Dependencies:
%   * create_NM.m
%   * create_NM_network.m
%   * sim_NM_network.m
%
% This code is part of the simulation code of the manuscript 'Pathological responses to single pulse electrical stimuli in epilepsy: the role of feedforward inhibition'
% (c) 2019 Jurgen Hebbink (University of Twente, University Medical Centre Utrecht)
%%
tic;            % Start measuring computation time
%% Simulation settings
deltat=0.0001;  % Timestep for ODE solving
Tend=60;         % End time of the simulation
Np=5;           % Save every Np-th datapoint
%% Stimulation settings
Tin=3;          % Time before the first stimulation
Tinterstim=5;   % Time between stimulations
%Blockpulse
Tstim=0.005;    % Length of the blockpulse
Amp=1500;       % Amplitude of the blockpulse

%% Initialize neural masses

% Create parameter structs
% strout=create_NM(A,a,B,b,G,g,C,sd,alpha,beta,gamma)
NM(1)=create_NM(4.5,100,7,10,25,300,135,1,0,1,0.7);
NM(2)=create_NM(4.5,100,7,10,25,300,135,1,0,1,0.7);

%Add stimulation
NM(1).Ivar=@(t) (mod(t,Tinterstim)<Tin+Tstim).*(mod(t,Tinterstim)>=3)*(Amp);    % Stimulation at NM 1
% NM(2).Ivar=@(t)(mod(t,Tinterstim)<Tin+Tstim).*(mod(t,Tinterstim)>=3)*(Amp);   % Uncomment to add stimulation to NM2

%% Network architecture
k=20;           % Connectivity strength
cm=[0 0;1 0];   % Connectivity matrix, cm(i,j) is connection j->i.

NM=NM(:);
netw=create_NM_network(NM,cm,k);
%strout=create_NM_network(nodes, conmat, constrength)
%% perform simulation
%rng(0);                                                            % Reset random number generator if desired
[u2, x, xext, y] = sim_NM_network(Tend, deltat, Np, netw);          % Performs actual simulation
%[uout,xout,xextout,yout,yextout] = sim_NM_network(Ttot,deltat,Nout,NM_network)

%% show results
tvec=0:deltat*Np:Tend;      % Create time vector
figure(1);
plot(tvec,u2([1:4:end],:)') % Plots u_PY of Both NMM1(u2(1)) and NMM2(u2(5))
% So rows 1:4:end contain the potential of the pyramidal cells.
xlabel ('Time (s)')
ylabel ('Potential')
legend('NM1','NM2')
title('Simulated Potential of Pyramidal cells of NM1 and NM2')
%%
toc;    % Stop measuring computation time