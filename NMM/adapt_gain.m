B_gain = [0, 12, 25, 50, 75, 100];

for i = B_gain
    for j = 1:length(B_gain)
        name(j) = B_gain(j);
        NM(j)=create_NM(4.5,100,7,10,i,300,135,1,0,1,0.7); 





        %Add stimulation
        NM(j).Ivar=@(t) (mod(t,Tinterstim)<Tin+Tstim).*(mod(t,Tinterstim)>=3)*(Amp);    % Stimulation at NM 1
        %NM(2).Ivar=@(t)(mod(t,Tinterstim)<Tin+Tstim).*(mod(t,Tinterstim)>=3)*(Amp);   % Uncomment to add stimulation to NM2

        %% Network architecture
        k=20;           % Connectivity strength
        cm=[0 0;1 0];   % Connectivity matrix, cm(i,j) is connection j->i.

        NM(j)=NM(:);
        netw(j)=create_NM_network(NM(j),cm,k);
        %strout=create_NM_network(nodes, conmat, constrength)
        %% perform simulation
        %rng(0);                                                            % Reset random number generator if desired
        [u2, x, xext, y] = sim_NM_network(Tend, deltat, Np, netw(j));          % Performs actual simulation
        %[uout,xout,xextout,yout,yextout] = sim_NM_network(Ttot,deltat,Nout,NM_network)

        %% show results
        tvec=0:deltat*Np:Tend;      % Create time vector
        figure(j);
        plot(tvec,u2([1:4:end],:)') % Plots u_PY of Both NMM1(u2(1)) and NMM2(u2(5))
        % So rows 1:4:end contain the potential of the pyramidal cells.
        xlabel ('Time (s)')
        ylabel ('Potential')
        legend('NM1','NM2')
        %title('Simulated Potential of Pyramidal cells of NM1 and NM2')
        title(sprintf('Simulated Potential PY of NM1 Bgain= %s', 25))
    end
end


