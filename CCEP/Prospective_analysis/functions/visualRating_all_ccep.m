
function dataBase = visualRating_all_ccep(dataBase)


% INSTRUCTIONS
% select new point: select new point > enter or y
% correct: y
% incorrect: n or enter

close all

tt = dataBase.tt;

 n1_peak_amplitude = dataBase.ccep.n1_peak_amplitude;
 n1_peak_sample = dataBase.ccep.n1_peak_sample;
% 
% Preallocation
n1_peak_amplitude_all_check = NaN(size(n1_peak_amplitude));
n1_peak_sample_all_check = NaN(size(n1_peak_sample));

for stimp = 1:size(dataBase.cc_epoch_sorted_avg,2)
        
    for chan =1 :size(dataBase.cc_epoch_sorted_avg,1)
        
%         if ~isnan(dataBase.ccep.n1_peak_sample(chan,stimp))
            % figure with left the epoch, and right zoomed in
            H=figure(1);
            H.Units = 'normalized';
            H.Position = [0.13 0.31 0.77 0.7];
            this_plot = squeeze(dataBase.cc_epoch_sorted_select_avg(chan,stimp,:,:));                                
            this_plot(:,tt>-0.001 & tt<0.01) = NaN;            
            
            this_plot_avg = squeeze(dataBase.cc_epoch_sorted_avg(chan,stimp,:));
            this_plot_avg(tt>-0.001 & tt<0.01) = NaN;                      
               
            subplot(1,2,1)
            plot(tt,this_plot,':r','linewidth',1);
            hold on
            plot(tt,this_plot_avg,'k','linewidth',2);
            
            % If a ER was detected with the detector, then plot a blue dot
            % at the N1 peak
            if ~isnan(dataBase.ccep.n1_peak_sample(chan,stimp))
                plot(tt(n1_peak_sample(chan,stimp)),this_plot_avg(n1_peak_sample(chan,stimp)),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3)
            end
            
            hold off
            xlim([-2 2])
            ylim([-2000 2000])
            xlabel('time(s)')
            ylabel('amplitude(uV)')
            title(sprintf('Electrode %s, stimulating %s',dataBase.ch{chan},dataBase.stimpnames_avg{stimp}))
            
            subplot(1,2,2)
            plot(tt,this_plot,':r','linewidth',1);
            hold on
            plot(tt,this_plot_avg,'k','linewidth',2);
            
            % If a ER was detected with the detector, then plot a blue dot
            % at the N1 peak
            if ~isnan(dataBase.ccep.n1_peak_sample(chan,stimp))
                plot(tt(n1_peak_sample(chan,stimp)),this_plot_avg(n1_peak_sample(chan,stimp)),'o','MarkerEdgeColor','b','MarkerFaceColor','b')
            end
            
            hold off
            xlim([-0.2 0.5])
            ylim([-750 750])
            title('Zoomed average signal')
            xlabel('Time (s)')
            ylabel('Voltage (uV)')
            
            currkey = 0;
            fprintf('When a blue dot is already shown, press y/n. No ER, press n, When ER, select the correct N1 and press enter. \n')
            
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
                    [~,sampnum] = min(abs(tt-cp(1,1)));
                    
                    % find nearby peak
                    [~,locs] = findpeaks(-1*this_plot_avg(sampnum-50:sampnum+50),...
                        'NPeaks',1,'SortStr','descend');
                    
                    % find x-position of nearby peak
                    locsamp = sampnum-50+locs-1;
                    
                    hold on
                    plot(tt(locsamp),this_plot_avg(locsamp),'bo','MarkerFaceColor','b','MarkerSize',4); drawnow;
                    hold off
                    
                    n1_peak_sample_all_check(chan,stimp) = locsamp ;
                    n1_peak_amplitude_all_check(chan,stimp) = this_plot_avg(locsamp) ;
                    
                elseif w == 1
                    currkey = get(gcf,'CurrentCharacter');
                    
                    if strcmp(currkey,'y') && isempty(cp)
                        n1_peak_amplitude_all_check(chan,stimp) = n1_peak_amplitude(chan,stimp) ;
                        n1_peak_sample_all_check(chan,stimp) = n1_peak_sample(chan,stimp) ;
                    elseif strcmp(currkey,'n')
                        n1_peak_amplitude_all_check(chan,stimp) = NaN ;
                        n1_peak_sample_all_check(chan,stimp) = NaN ;
                    end
                end
            end
    end
end

dataBase.ccep.n1_peak_amplitude_all_check = n1_peak_amplitude_all_check;
dataBase.ccep.n1_peak_sample_all_check = n1_peak_sample_all_check;
