
clear;
ccep_allPat.mode = 'NMM';
myDataPath = setLocalDataPath(ccep_allPat);

%% Simulation settings
deltat=0.0001;  % Timestep for ODE solving
Tend=60;         % End time of the simulation
Np=1;           % Save every Np-th datapoint

%% Stimulation settings
Tin=3;          % Time before the first stimulation
Tinterstim=5;   % Time between stimulations
%Blockpulse
Tstim=0.005;    % Length of the blockpulse
Amp=1500;       % Amplitude of the blockpulse

SI_gain = repmat(7.0, [1,31]);                     % Slow inhibitory synaptic gain (norm = 7 mV)
SI_reci = repmat(4.6, [1,31]);                  % Reciprocal of slow inhibitory time constant (norm = 10 s-1)

FI_gain = [10 :0.5: 25 ] ;  %17;     repmat(25, [1,11]);     % Fast inhibitory synaptic gain (norm = 25 mV)
FI_reci = [145:5: 295] ; % 300;    %repmat(300, [1,11])];                               % Reciprocal of fast inhibitory time constant (norm = 300 s-1)

EX_gain = repmat(4.5, [1,31]);       %[3:0.3:5];                                 % Excitatory synaptic gain (norm = 4.5 mV)
EX_reci = repmat(100, [1,31]);       %[70 :5:100];                      %repmat(100, [1,31]) [60 :5:110]]);   % Excitatory time constant (norm = 100 s-1)

%% Varying the Fast Inhibitory Gain value
% Simulating the potential of the pyramidal cells of NM1 and NM2
gcf = figure('Position',[272,489,930,533]);
j = 1;
n = size(EX_gain,2);
c = flipud(jet(n));

for i = 1:size(EX_gain,2)   
    
     for j = 1:size(EX_reci,2) 
        
        NM(1)=create_NM(4.5,100,7,4.0,11,175,135,1,0,1,0.7);         %EX_gain(i),EX_reci(i),SI_gain(i),SI_reci(i),FI_gain(i),FI_reci(j)  % This is NMM 1      create_NM(A,a,B,b,G,g,C,sd,alpha,beta,gamma)
        NM(2)=create_NM(4.5,100,7,10,25,300,135,1,0,1,0.7);          % This is NMM 2

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

        fs = 600001/60;
        ts = 1/fs;
        tstart = round(2.98 * fs);
        tend = round(3.4 *fs);
%         
%         if i > 10
            p(i,j) = plot(tvec(tstart:tend),u2(1,tstart:tend)','LineWidth',3,'Color','b'); 
            hold on
%         else
%             p(i) = plot(tvec(tstart:tend),u2(1,tstart:tend)','LineWidth',1,'Color',[0, 0.4470, 0.7410]); 
%         end
        
%         xlim([2.98 3.08])
        ylim([-50 10])
        xlabel ('Time (ms)','FontSize',14,'FontWeight','bold')
        ylabel ('Potential (mV)','FontSize',14,'FontWeight','bold')
        title('Four times more noise')
        ax = gca;                      
%         ax.XTickLabel = [ax.XTick(1) :0.01: ax.XTick(end)] - 3;
        ax.FontSize = 12;
        hold on
        
%         legendCell = cellstr(num2str(EX_reci', 'a=%-d'));
%         legend(legendCell)
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
         
    end
        

end



% Save settings and outcome
NMM.settings.SI_gain = SI_gain;
NMM.settings.SI_timeC = SI_reci;
NMM.settings.FI_gain = FI_gain;
NMM.settings.FI_timeC = FI_reci;
NMM.settings.EX_gain = EX_gain;
NMM.settings.EX_timeC = EX_reci;

% Save in ms
xN1_plot = (xN1-3)*1000;
xP1_plot = (xP1-3)*1000;
xN2_plot = (xN2-3)*1000;

NMM.outcome.N1_lat = xN1_plot;
NMM.outcome.P1_lat = xP1_plot;
NMM.outcome.N2_lat = xN2_plot;
  

fileName = ['NMM_setting_and_outcome_',datestr(now, 'dd-mmm-yy, HH:MM'),'.mat'];
save([myDataPath.save_data_loc,fileName], 'NMM');

%% Plot the results of adaptations to the NMM in the colourplot
% The gain is on the y-axis, the variations of the time-constant are
% displayed with various colors. 
% x-axis is the latency (ms)

if exist('NMM','var')      % after the script above is ran
   colourPlot_latencies(NMM, myDataPath)

    
else % load save data before running colourPlot_latencies
    
    files = dir(fullfile(myDataPath.load_data_loc));
    Loc_NMMFile = find(contains({files(:).name},'NMM'));
    load(fullfile(files(Loc_NMMFile).folder,files(Loc_NMMFile).name));
    
    colourPlot_latencies(NMM, myDataPath)

end


%% Plot the contour plot of the latencies
% Plot the results of the NMM in a contour plot. The time constant on the
% x-axis, the gain on the y-axis. A logarithmic scale is used to indicate
% the deviation from the desirec value (i.e. the value found during the in
% vivo studies).

if exist('NMM','var')      % after the script above is ran
    
    contourPlot_latencies(NMM, myDataPath)
    
else % load save data before running colourPlot_latencies
    
    files = dir(fullfile(myDataPath.load_data_loc));
    Loc_NMMFile = find(contains({files(:).name},'NMM'));
    load(fullfile(files(Loc_NMMFile).folder,files(Loc_NMMFile).name));
    
    contourPlot_latencies(NMM, myDataPath)

end

