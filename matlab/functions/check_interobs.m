
function dataBase = check_interobs(dataBase, myDataPath)

close all
T = readtable(fullfile(myDataPath.CCEP_interObVar,"Different_ratings_PRIOS04_REC2Stim05.xlsx"));

tt = dataBase.tt;

n1_peak_amplitude = dataBase.ccep.n1_peak_amplitude;
n1_peak_sample = dataBase.ccep.n1_peak_sample;

% Preallocation
n1_peak_amplitude_check_check = NaN(size(n1_peak_amplitude));
n1_peak_sample_check_check = NaN(size(n1_peak_sample));

% Find stimpair column number
for stimp = 1:size(T,1)
    Col_stimp = find(T{stimp,3:4} == dataBase.cc_stimsets_avg(:,1:2));
    T{stimp,5} = Col_stimp(1);
end

for Event = 1:size(T,1)                           % For each stimulation pair
    
    stimp = T{Event,5};
    chan =T{Event,1}; 
                                       
            H=figure(1);
            H.Units = 'normalized';
            H.Position = [0.13 0.31 0.77 0.7];
            this_plot = squeeze(dataBase.cc_epoch_sorted_select_avg(chan,stimp,:,:));                   
            this_plot_avg = squeeze(dataBase.cc_epoch_sorted_avg(chan,stimp,:));
               
            subplot(1,2,1)
            plot(tt,this_plot,':r','linewidth',1);
            hold on
            plot(tt,this_plot_avg,'k','linewidth',2);
            plot(tt(n1_peak_sample(chan,stimp)),this_plot_avg(n1_peak_sample(chan,stimp)),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3)
            hold off
            xlim([-0.2 0.5]); ylim([-1000 1000]); xlabel('time(s)'); ylabel('potential (\muV)');
            title(sprintf('Electrode %s, stimulating %s',dataBase.ch{chan},dataBase.stimpnames_avg{stimp}))
            
            
            subplot(1,2,2)
            plot(tt,this_plot,':r','linewidth',1);
            hold on
            plot(tt,this_plot_avg,'k','linewidth',2);
            plot(tt(n1_peak_sample(chan,stimp)),this_plot_avg(n1_peak_sample(chan,stimp)),'o','MarkerEdgeColor','b','MarkerFaceColor','b')
            hold off
            xlim([-0.02 0.09]); ylim([-400 400]); title('Zoomed average signal'); xlabel('Time (s)'); ylabel('Potential (\muV)')
            
            % Create patch to indicate the 9 ms interval
            patch([0 0.009 0.009 0],[-400 -400 1000 1000],[0.6,0.2,0.2], 'EdgeAlpha',0)
            alpha(0.1)

            currkey = 0;
            fprintf('N1 [y/n], if incorrect N1, select correct N1 and press enter \n')
            
            % select new N1 or categorize as good N1 or no N1
            % When incorrect N1 is selected, click on correct N1, a blue
            % stip will occure, then press enter! The new coordinates will
            % show in n1_peak_amplitude and sample.
            while ~strcmp(currkey,{'y','n',char(13)})
                cp =[];
                w = waitforbuttonpress; % 0 = mouse, other = key
                if w == 0
                    % draw correct N1
                    cp = get(gca,'CurrentPoint');
                                                              
                    % find sample number closest to the selected point
                    % Try to select the N1-peak as correctly as possible!
                    [~,sampnum] = min(abs(tt-cp(1,1)));
                    
                    locsamp = sampnum;
                    
                    hold on
                    plot(tt(locsamp),this_plot_avg(locsamp),'bo','MarkerFaceColor','r','MarkerSize',4); drawnow;
                    hold off
                    
                    n1_peak_sample_check_check(chan,stimp) = locsamp ;
                    n1_peak_amplitude_check_check(chan,stimp) = this_plot_avg(locsamp) ;
                    
                elseif w == 1
                    currkey = get(gcf,'CurrentCharacter');
                    
                    if strcmp(currkey,'y') && isempty(cp)
                        n1_peak_amplitude_check_check(chan,stimp) = n1_peak_amplitude(chan,stimp) ;
                        n1_peak_sample_check_check(chan,stimp) = n1_peak_sample(chan,stimp) ;
                    elseif strcmp(currkey,'n')
                        n1_peak_amplitude_check_check(chan,stimp) = NaN ;
                        n1_peak_sample_check_check(chan,stimp) = NaN ;
                    end
                end
            end
            
end

dataBase.ccep.n1_peak_amplitude_check_check = n1_peak_amplitude_check_check;
dataBase.ccep.n1_peak_sample_check_check = n1_peak_sample_check_check;

%% Visualise the difference of latency of observer and detector
% Find the number of samples after the stimulation artefact
% 9 ms is 19 samples
% Fs = 2048 
Latency_interOb = dataBase.ccep.n1_peak_sample_check_check-(2*2048);            % To samples after stimulation artefact
Num_interOb = find(~isnan(Latency_interOb));
Latency_interOb(isnan(Latency_interOb)) = [];
numel(find(Latency_interOb<19))

% Latencies of detector detected ERs
Latency_detector = dataBase.ccep.n1_peak_sample-(2*2048);
Latency_detector = Latency_detector(Num_interOb)'; %#ok<FNDSB>

figure()
boxplot([Latency_interOb' Latency_detector'],'Notch','on', ...
        'Labels',{'Checked Latencies N1','Detector latencies N1'})
ylabel('Number of samples after stimulation artefact')

% Save
outlabel=('N1_latency_interObs.jpg');
path = fullfile(myDataPath.CCEP_interObVar,'N1_latency_interOb/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'jpg')

end


