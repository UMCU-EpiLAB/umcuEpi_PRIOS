
clear;
close all;
tic;
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

FI_gain = [0, 12, 25, 50, 75, 100];
SI_gain = [-7, 0, 7, 14, 21, 28];

%% Varying the Fast Inhibitory Gain value
% Simulating the potential of the pyramidal cells of NM1 and NM2

for i = FI_gain

        NM(1)=create_NM(4.5,100,7,10,i,300,135,1,0,1,0.7);          % This is NMM 1
        NM(2)=create_NM(4.5,100,7,10,i,300,135,1,0,1,0.7);        % This is NMM 2

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
        gcf = figure();
        plot(tvec,u2([1:4:end],:)') % Plots u_PY of Both NMM1(u2(1)) and NMM2(u2(5))
        % So rows 1:4:end contain the potential of the pyramidal cells.
        xlabel ('Time (s)')
        ylabel ('Potential')
        legend('NM1','NM2')
        title(['Simulated Potential PY, G-gain = ',sprintf('%d',i)])
        %title('Simulated Potential of Pyramidal cells of NM1 and NM2')
        outlabel=sprintf('Simulated Potential PY %d.jpg',...
        i);

        path = 'C:\Users\Sifraaaaa\Documents\utwente\M3\NMM\GandB_changes\';
        saveas(gcf,[path,'/FI_gain_(G)/',outlabel],'jpg')
        saveas(gcf,[path,'/FI_gain_(G)/',outlabel])
end

%% Varying the Slow Inhibitory Gain value
% Simulating the potential of the pyramidal cells of NM1 and NM2

for i = SI_gain

        NM(1)=create_NM(4.5,100,i,10,25,300,135,1,0,1,0.7);          % This is NMM 1
        NM(2)=create_NM(4.5,100,i,10,25,300,135,1,0,1,0.7);        % This is NMM 2

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
        gcf = figure();
        plot(tvec,u2([1:4:end],:)') % Plots u_PY of Both NMM1(u2(1)) and NMM2(u2(5))
        % So rows 1:4:end contain the potential of the pyramidal cells.
        xlabel ('Time (s)')
        ylabel ('Potential')
        legend('NM1','NM2')
        title(['Simulated Potential PY, B-gain = ',sprintf('%d',i)])
        %title('Simulated Potential of Pyramidal cells of NM1 and NM2')
        outlabel=sprintf('Simulated Potential PY %d.jpg',...
        i);

        path = 'C:\Users\Sifraaaaa\Documents\utwente\M3\NMM\GandB_changes\';
        saveas(gcf,[path,'/SI_gain_(B)/',outlabel],'jpg')

end

toc;
