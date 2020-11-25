
clear;
close all;

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

FI_gain = [15, 25, 40, 50, 60];                 % Fast inhibitory synaptic gain (norm = 25 mV)
SI_gain = [-7, 0, 7, 14, 28];                   % Slow inhibitory synaptic gain (norm = 7 mV)
SI_reci = [0, 5, 10, 15, 20];                   % Reciprocal of slow inhibitory time constant (norm = 10 s-1)
FI_reci = [150, 250, 300, 350 , 450];           % Reciprocal of fast inhibitory time constant (norm = 300 s-1)


%% Varying the Fast Inhibitory Gain value
% Simulating the potential of the pyramidal cells of NM1 and NM2
gcf = figure('Position',[414,688,500,374]);

for i = 1:size(FI_gain,2)

        NM(1)=create_NM(4.5,100,7,10,FI_gain(i),300,135,1,0,1,0.7);          % This is NMM 1      create_NM(A,a,B,b,G,g,C,sd,alpha,beta,gamma)
        NM(2)=create_NM(4.5,100,7,10,FI_gain(i),300,135,1,0,1,0.7);          % This is NMM 2

        %Add stimulation
        NM(1).Ivar=@(t) (mod(t,Tinterstim)<Tin+Tstim).*(mod(t,Tinterstim)>=3)*(Amp);    % Stimulation at NM 1

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

        plot(tvec,u2([1],:)','LineWidth',2)
        xlim([2.98 3.2])
        ylim([-80 10])
        xlabel ('Time (s)')
        ylabel ('Potential')
        legend('G = 15','G = 25','G = 40','G = 50','G = 60','Location','southeast')
     
            
        title(['Simulated Potential PY, FI gain'])
        hold on
        
end

%% Varying the Slow Inhibitory Gain value
% Simulating the potential of the pyramidal cells of NM1 and NM2
gcf = figure('Position',[414,688,500,374]);

for i = 1:size(SI_gain,2)

        NM(1)=create_NM(4.5,100,SI_gain(i),10,25,300,135,1,0,1,0.7);          % This is NMM 1
        NM(2)=create_NM(4.5,100,SI_gain(i),10,25,300,135,1,0,1,0.7);          % This is NMM 2

        %Add stimulation
        NM(1).Ivar=@(t) (mod(t,Tinterstim)<Tin+Tstim).*(mod(t,Tinterstim)>=3)*(Amp);    % Stimulation at NM 1

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
        
        plot(tvec,u2([1],:)','LineWidth',2)
        xlim([2.98 3.5])
        ylim([-50 8])
        xlabel ('Time (s)')
        ylabel ('Potential')
        legend('B = -7','B = 0','B = 7','B = 14','B = 28','Location','southeast')
            
        title(['Simulated Potential PY, SI gain'])
        hold on


end

