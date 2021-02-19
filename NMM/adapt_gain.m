
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

SI_gain = [repmat(7.6, [1,10]), repmat(7.9, [1,10])];                     % Slow inhibitory synaptic gain (norm = 7 mV)
SI_reci = [repmat(4.6, [1,10]), repmat(4.6, [1,10])]   ;                  % Reciprocal of slow inhibitory time constant (norm = 10 s-1)

FI_gain = [repmat(25, [1,10]), repmat(16.3, [1,10])]    ;               % Fast inhibitory synaptic gain (norm = 25 mV)
FI_reci = [repmat(300, [1,10]),repmat(180, [1,10])];                               % Reciprocal of fast inhibitory time constant (norm = 300 s-1)8


%% Varying the Fast Inhibitory Gain value
% Simulating the potential of the pyramidal cells of NM1 and NM2
gcf = figure('Position',[272,203,1220,819]);
j = 1;


for i = 1:size(SI_gain,2)   
    
%     for j = 1:size(SI_reci,2) 
        
        NM(1)=create_NM(4.5,100, SI_gain(i),SI_reci(i),  FI_gain(i),FI_reci(i) , 135,1,0,1,0.7);          % This is NMM 1      create_NM(A,a,B,b,G,g,C,sd,alpha,beta,gamma)
        NM(2)=create_NM(4.5,100,SI_gain(i),SI_reci(i),FI_gain(i),FI_reci(i),135,1,0,1,0.7);          % This is NMM 2

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

        fs = 120001/60;
        ts = 1/fs;
        tstart = round(2.98 * fs);
        tend = round(3.4 *fs);
%         
        if i <11
            p(i) = plot(tvec(tstart:tend),u2(1,tstart:tend)','LineWidth',1,'Color',[0.6350, 0.0780, 0.1840]); 
        else %i>=11
            p(i) = plot(tvec(tstart:tend),u2(1,tstart:tend)','LineWidth',1,'Color',[0, 0.4470, 0.7410]); 
        end
        
%         xlim([2.98 3.08])
        ylim([-50 10])
        xlabel ('Time (ms)','FontSize',14,'FontWeight','bold')
        ylabel ('Potential mV','FontSize',14,'FontWeight','bold')
        
        ax = gca;                      
%         ax.XTickLabel = [ax.XTick(1) :0.01: ax.XTick(end)] - 3;
        ax.FontSize = 12;
        hold on
%        

        %% Determine the N1
        [yN, xN] = findpeaks(-u2(1,tstart:tend)','MinPeakHeight',-10,'MinPeakDistance',1);
        
        if ~isempty(yN)
            
            if numel(yN) > 1            % Select the N1 with the higest amplitude (yN1)
                Loc_max = find(yN == max(yN));
                xN = xN(Loc_max);
                yN = -max(yN);
                
                xN = (tstart*ts)+(xN*ts);
                yN1(i,j) = yN;
                xN1(i,j) = xN;
            else
                xN = (tstart*ts)+(xN*ts);
                yN1(i,j) = yN;
                xN1(i,j) = xN;
            end
                
        else
            yN1(i,j) = NaN;
            xN1(i,j) = NaN;
            
        end
        
        % Determine the P1
        [yP, xP] = findpeaks(u2(1,tstart:tend)','MinPeakHeight',-20);
          
        if ~isempty(yP)

            if numel(yP) >1         % When there are more than 1 peaks
                
                % Convert the samples to seconds to compare the P1 with the
                % N1 found above
                xP1_after_n1 = (tstart*ts)+(xP*ts)';

                % When there is a P1 after the N11
                if any(xP1_after_n1 > xN1(i,j))
                
                    % Determine the P1 after the N1
                    loc_P1 = min(find(xP1_after_n1 > xN1(i,j)));
                    xP1(i,j) = xP1_after_n1(loc_P1);
                    yP1(i,j) = yP(loc_P1); 


                    % for the P0, find the highest peak before the N1.
                    % Loc_max = find(yP == max(yP)); 
                else 
                     xP1(i,j) = NaN;                yP1(i,j) = NaN;            
                end   
            else
                xP = (tstart*ts)+(xP*ts);
                xP1(i,j) = xP;
                yP1(i,j) = yP;
                
            end
        else
            xP1(i,j) = NaN;                yP1(i,j) = NaN;            
        end
        
        
        
        %%% Determine the N2
        %N2 must be later than the P1         
         [yN, xN] = findpeaks(-u2(1,tstart:tend)','MinPeakHeight',-2,'MinPeakDistance',1);
        
        if ~isempty(yN)

            if numel(yN) >1         % When there are more than 1 peaks
                % The N2 must be later than the P1
                % convert sample to seconds, to compare with the found P1
                % above
                xN2_after_P1 = (tstart*ts)+(xN*ts)';
                
                % When there is a P1 after the N11
                if any(xN2_after_P1 > xP1(i,j))
                
                    % Determine the biggest yP1 after the N1
                    
                    loc_N2 = find(xN2_after_P1 > xP1(i,j));         % find the points that are after the N1
                    yN2(i,j) =  -max(yN(loc_N2));                    % Find the largest (negative) ampliude of the P1
                                                                    
                    loc_yN2 = find(ismember(-yN , yN2(i,j)));        % Determine the xN2
                    xN2(i,j) = xN2_after_P1(loc_yN2);


                    % for the P0, find the highest peak before the N1.
                    % Loc_max = find(yP == max(yP)); 
                    
                else 
                     xN2(i,j) = NaN;                yN2(i,j) = NaN;            

                end
                
            
            else
                xN = (tstart*ts)+(xN*ts);
                xN2(i,j) = xN;
                yN2(i,j) = yN;
                
            end
 
        else
            xN2(i,j) = NaN;                yN2(i,j) = NaN;            
        end
         
%     end
        

end


%
% N1
xN1C = mean(xN1(1:10));
yN1C = mean(yN1(1:10));

xN1P = mean(xN1(11:20));
yN1P = mean(yN1(11:20));

%P1
xP1C = mean(xP1(1:10));
yP1C = mean(yP1(1:10));

xP1P = mean(xP1(11:20));
yP1P = mean(yP1(11:20));

%N2
xN2C = mean(xN2(1:10));
yN2C = mean(yN2(1:10));

xN2P = mean(xN2(11:20));
yN2P = mean(yN2(11:20));


hold on
p(size(SI_gain,2) +1) = plot(xN1C,yN1C,'o','MarkerSize',9,'MarkerEdgeColor',[0.93,0.69,0.13],'MarkerFaceColor',[0.93,0.69,0.13]);
p(size(SI_gain,2) +2) = plot(xN1P,yN1P,'o','MarkerSize',9,'MarkerEdgeColor',[0.93,0.89,0.13],'MarkerFaceColor',[0.93,0.89,0.13]);

p(size(SI_gain,2) +3) = plot(xP1C,yP1C,'o','MarkerSize',9,'MarkerEdgeColor',[0.99,0.00,0.99],'MarkerFaceColor',[0.99,0.00,0.99]);
p(size(SI_gain,2) +4) = plot(xP1P,yP1P,'o','MarkerSize',9,'MarkerEdgeColor',[0.99,0.80,0.99],'MarkerFaceColor',[0.99,0.80,0.99]);

p(size(SI_gain,2) +5) = plot(xN2C,yN2C,'o','MarkerSize',9,'MarkerEdgeColor',[0.47,0.57,0.19],'MarkerFaceColor',[0.47,0.57,0.19]);
p(size(SI_gain,2) +6) = plot(xN2P,yN2P,'o','MarkerSize',9,'MarkerEdgeColor',[0.47,0.7,0.29],'MarkerFaceColor',[0.47,0.87,0.29]);



legendInfo{1} = [sprintf('Clin: B = %1.1f, b = %1.1f, G = %1.1f, g = %1.0f',SI_gain(1), SI_reci(1),FI_gain(1),FI_reci(1))]    ;   
legendInfo{2} = [sprintf('Prop: B = %1.1f, b = %1.1f, G = %1.1f, g = %1.0f',SI_gain(11), SI_reci(11),FI_gain(11),FI_reci(11))]    ;   

legendInfo{3} = sprintf('N1-clin: Lat = %1.0f ms, Amp = %1.1f mV',(xN1C-3)*1000,yN1C);
legendInfo{4} = sprintf('N1-prop: Lat = %1.0f ms, Amp = %1.1f mV',(xN1P-3)*1000,yN1P);

legendInfo{5} = sprintf('P1-clin: Lat = %1.0f ms, Amp = % 1.1f mV',(xP1C-3)*1000,yP1C);
legendInfo{6} = sprintf('P1-prop: Lat = %1.0f ms, Amp = %1.1f mV',(xP1P-3)*1000,yP1P);

legendInfo{7} = sprintf('N2-clin: Lat = %1.0f ms, Amp = %1.1f mV',(xN2C-3)*1000,yN2C);
legendInfo{8} = sprintf('N2-prop: Lat = %1.0f ms, Amp = %1.1f mV',(xN2P-3)*1000,yN2P);



legend(p([1,11, 21:end]),legendInfo,'Location','Southeast')
xlim([2.99,3.4])
% ax.XTickLabel = [-10, 0,10,20,30,40,50,60,70,80,90];
















