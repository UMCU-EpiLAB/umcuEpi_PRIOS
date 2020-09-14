function dataBase = plot_all_signals(dataBase)
%
% Sifra Blok, 2020, UMC Utrecht & University of Twente
% Dorien van Blooijs, 2020, UMC Utrecht


check_allSig_amplitude = NaN(size(dataBase.ch,1),size(dataBase.stimpnames_avg,2));
check_allSig_sample = NaN(size(dataBase.ch,1),size(dataBase.stimpnames_avg,2));
                                      
                    
tt = dataBase.tt;    
   for stimp = 1:size(dataBase.cc_epoch_sorted_avg,2)            % for all stimulation pairs
         Stimpnm = dataBase.stimpnames_avg{stimp};    
        for elec = 1:size(dataBase.cc_epoch_sorted_avg,1)            % for all electrodes
             elecnm = dataBase.ch{elec};
                ccep_plot = squeeze(dataBase.cc_epoch_sorted_avg(elec,stimp,:));
                ccep_plot(tt>-0.010 & tt<0.010) = NaN;

                figure()
                plot(tt,  ccep_plot);   
                str = sprintf('Stimulation pair %s on %s', Stimpnm, elecnm);
                title(str)
                xlim([-0.2 0.5])
                ylim([-750 750])
                ylabel('Average per electrodes (mV)')
                xlabel('time (s)') 
                set(gca,'YTick',1,'YTickLabel',elecnm) ;
                
                currkey = 0;
                fprintf('When no N1, press N, otherwise click on N1, press enter \n')
            
            while ~strcmp(currkey,{'n',char(13)})
                cp =[];
                w = waitforbuttonpress;
                                    
                if w == 0
                    % click on  N1
                    cp = get(gca,'CurrentPoint');
                   
                    % find sample number closest to the selected point
                    [~,sampnum] = min(abs(tt-cp(1,1)));
                    
                     % find nearby peak
                    [~,locs] = findpeaks(-1*ccep_plot(sampnum-50:sampnum+50),...
                        'NPeaks',1,'SortStr','descend');
                    
                    % find x-position of nearby peak
                    locsamp = sampnum-50+locs-1;
                    
                    hold on
                    plot(tt(locsamp),ccep_plot(locsamp),'bo','MarkerFaceColor','b','MarkerSize',4); drawnow;
                    hold off
                    
                    check_allSig_sample(elec,stimp) = locsamp;
                    check_allSig_amplitude(elec,stimp) = ccep_plot(locsamp) ;
                    
               elseif w == 1                % when no ER is detected ('n')
                   check_allSig_amplitude(elec,stimp) = NaN ;
                   check_allSig_sample(elec,stimp) = NaN ;
                end
                                      
            end
        end            
    end
    
    close all
   
   dataBase.ccep.check_allSig_amplitude = check_allSig_amplitude;
   dataBase.ccep.check_allSig_sample = check_allSig_sample;
end

          
