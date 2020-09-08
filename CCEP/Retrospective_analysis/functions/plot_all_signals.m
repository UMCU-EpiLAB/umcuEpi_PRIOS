function dataBase = plot_all_signals(dataBase)
%
% Sifra Blok, 2020, UMC Utrecht & University of Twente
% Dorien van Blooijs, 2020, UMC Utrecht

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
                fprintf('N1 [y/n], if incorrect N1, select correct N1 and press enter \n')
            
            while ~strcmp(currkey,{'y','n',char(13)})
                cp =[];
                w = waitforbuttonpress;
                                    
                if w == 1
                    currkey = get(gcf,'CurrentCharacter');
                    
                    ER_check_amplitude = zeros(size(dataBase.ch,1),size(dataBase.stimpnames_avg,2));
                    ER_check_sample = zeros(size(dataBase.ch,1),size(dataBase.stimpnames_avg,2));
                    %%% dit klopt nog niet, hij slaat het niet op...
                    if strcmp(currkey,'y') && isempty(cp)
                        ER_check_amplitude(elec,stimp) = 1 ;
                        ER_check_sample(elec,stimp) = 1 ;
                    elseif strcmp(currkey,'n')
                        ER_check_amplitude(elec,stimp) = 0 ;
                        ER_check_sample(elec,stimp) = 0 ;
                    end
                end
            end
            
        end
    
    close all
   end 
   
   dataBase.ccep.ER_check_amplitude = ER_check_amplitude;
   dataBase.ccep.ER_check_sample = ER_check_sample;
end

          