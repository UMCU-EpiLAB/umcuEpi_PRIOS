function dataBase = visualRating_N2(dataBase, ccep)

% INSTRUCTIONS
% select new point: select new point > enter or y
% correct: y
% incorrect: n or enter

close all

tt = dataBase.tt;

N1_amplitude = ccep.n1_peak_amplitude_check;
N1_sample = ccep.n1_peak_sample_check;

% Preallocation
n2_amplitude = NaN(size(N1_amplitude));
n2_sample = NaN(size(N1_sample));

for stimp = 1:size(dataBase.cc_epoch_sorted_avg,2)
        
    for chan =1 :size(dataBase.cc_epoch_sorted_avg,1)
        
        if ~isnan(N1_amplitude(chan,stimp))              
           
            % figure with N1 in both clinical and propofol SPES
            H=figure(1);
            H.Units = 'normalized';
            H.Position = [0.13 0.31 0.77 0.7];
            this_plot = squeeze(dataBase.cc_epoch_sorted_select_avg(chan,stimp,:,:));           
             
            plot_avg = squeeze(dataBase.cc_epoch_sorted_avg(chan,stimp,:));
            plot_avg(tt>-0.001 & tt<0.01) = NaN;                      
               
            
            plot(tt,this_plot,':r','linewidth',1);
            hold on
            plot(tt,plot_avg,'k','linewidth',2);
            hold off
            xlim([-0.5 1])
            ylim([-2000 2000])
            xlabel('time(s)')
            ylabel('amplitude(uV)')
            title(sprintf('Electrode %s, stimulating %s',dataBase.ch{chan},dataBase.stimpnames_avg{stimp}))
            
           
            currkey = 0;
            fprintf('Select N2 and press enter, no N2 press n \n')
            
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
                    [~,locs] = findpeaks(-1*plot_avg(sampnum-50:sampnum+50),...
                        'NPeaks',1,'SortStr','descend');
                    
                    % find x-position of nearby peak
                    locsamp = sampnum-50+locs-1;
                    
                    hold on
                    plot(tt(locsamp),plot_avg(locsamp),'bo','MarkerFaceColor','b','MarkerSize',4); drawnow;
                    hold off
                    
                    n2_sample(chan,stimp) = locsamp ;
                    n2_amplitude(chan,stimp) = plot_avg(locsamp) ;
                    
                elseif w == 1
                    currkey = get(gcf,'CurrentCharacter');
                    
                    if strcmp(currkey,'y') && isempty(cp)
                        n2_amplitude(chan,stimp) = n1_peak_amplitude(chan,stimp) ;
                        n2_sample(chan,stimp) = n1_peak_sample(chan,stimp) ;
                    elseif strcmp(currkey,'n')
                        n2_amplitude(chan,stimp) = NaN ;
                        n2_sample(chan,stimp) = NaN ;
                    end
                end
            end
            
        end
        
    end
end

dataBase.ccep.n2_amplitude = n2_amplitude;
dataBase.ccep.n2_sample = n2_sample;